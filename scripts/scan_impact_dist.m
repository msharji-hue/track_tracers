function OUT = scan_impact_dist(root, varargin)
%SCAN_IMPACT_DIST  Diagnostic scan of impactDistPx. READ-ONLY.
%
%   Writes nothing, overwrites nothing, re-runs no kinematics.
%
%   Uses kin.rodBedDist_px, which the pipeline already computes and saves, so
%   nothing here re-derives geometry from tracks.
%
%   SIGN CONVENTION (from the saved calib: lineA=-36, lineB=0, lineC=144)
%   The bed line is vertical at x = 4 px and the normalised signed distance is
%   (4 - x) px. The projectile falls along +x, so rodBedDist_px DECREASES
%   through the drop. impactDistPx = -370 means contact is declared when the
%   reference marker reaches x = 374.
%
%   THE THREE QUESTIONS
%   Q1  Is the trigger value inside the recorded range at all?
%       If rodBedDist_px never crosses impactDistPx, the anchor pins to an
%       endpoint and both v0 and d are truncated.
%   Q2  Where does the trigger actually land relative to the record?
%       Few frames before impact_index = the trace starts mid-event.
%   Q3  What does the data say impactDistPx SHOULD be?
%       Peak speed marks contact (free fall accelerates, the bed decelerates)
%       and that criterion does not use impactDistPx. Reading rodBedDist_px at
%       the peak-speed frame gives an EMPIRICAL impactDistPx per trial.
%
%   VALIDATION
%   GB/full is known-good (v0/v_freefall = 1.01-1.05 in its mid-range). If the
%   peak-speed method does not recover roughly the saved value for GB/full,
%   the method is not working here and the rest of the scan is void.
%
%   USAGE
%       OUT = scan_impact_dist('D:\ME_GRANULAB\JerboaImpact');
%
%   OPTIONS
%       'Conditions'  default ["GB/full","CHIN/as_poured","CHIN/dense"]
%       'SmoothN'     frames for the speed smoother (default 5)
%       'MinPre'      pre-impact frames below which a trace is flagged (default 10)

opt.Conditions = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
opt.SmoothN    = 5;
opt.MinPre     = 10;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

fprintf('\n=== scan_impact_dist  (READ-ONLY) ===\n');

KN = dir(fullfile(root,'03_RESULTS','**','*_kin.mat')); KN = KN(~[KN.isdir]);
CS = dir(fullfile(root,'03_RESULTS','**','*_kin_scalars.csv')); CS = CS(~[CS.isdir]);
if isempty(KN), error('scan_impact_dist:noKin','No *_kin.mat under %s', root); end

scalarTags = string(erase({CS.name},'_kin_scalars.csv'))';
kinTags    = string(erase({KN.name},'_kin.mat'))';
orphan     = ~ismember(kinTags, scalarTags);
fprintf('_kin.mat: %d   _kin_scalars.csv: %d\n', numel(KN), numel(CS));
if any(orphan)
    fprintf('%d _kin.mat have NO matching scalars CSV (excluded from this scan):\n', sum(orphan));
    fprintf('   %s\n', kinTags(orphan));
end

% ------------------------------------------------------------- 1. per trial
rows = {};
for i = 1:numel(KN)
    if orphan(i), continue; end
    S = load(fullfile(KN(i).folder, KN(i).name));
    if ~isfield(S,'kin') || ~isfield(S,'meta'), continue; end
    kin = S.kin; meta = S.meta; cal = S.calib;

    cond = string(meta.material) + "/" + string(meta.container);
    if ~ismember(cond, opt.Conditions), continue; end
    if ~isfield(kin,'rodBedDist_px'), continue; end

    rb  = kin.rodBedDist_px(:);
    v   = kin.v(:);
    idp = cal.impactDistPx;
    iImp = kin.impact_index;

    ok = isfinite(rb) & isfinite(v);
    if nnz(ok) < 20, continue; end
    k = find(ok);

    % Q1 -- is the trigger value inside the recorded distance range?
    rbLo = min(rb(ok)); rbHi = max(rb(ok));
    inRange = (idp >= rbLo) && (idp <= rbHi);

    % Q3 -- empirical trigger from the peak-speed frame
    sp  = movmean(abs(v(k)), opt.SmoothN);
    [~, kPk] = max(sp);
    fPk = k(kPk);
    empIDP = rb(fPk);

    rows{end+1} = table(string(meta.trialTag), cond, ...
        local_h(meta), idp, ...
        rb(k(1)), rbLo, rbHi, inRange, ...
        iImp, sum(k < iImp), numel(k), ...
        fPk, empIDP, empIDP - idp, ...
        kin.v0_cm_s, kin.d_final_cm, kin.t_stop_s, kin.stopFrame, ...
        'VariableNames', {'trialTag','condition','h_true_mm','impactDistPx_used', ...
            'rb_first_px','rb_min_px','rb_max_px','trigger_in_range', ...
            'impact_index','nPreImpact','nTracked', ...
            'frame_vpeak','empirical_IDP_px','delta_px', ...
            'v0_cm_s','d_final_cm','t_stop_s','stopFrame'}); %#ok<AGROW>
end
if isempty(rows), error('scan_impact_dist:noMatch','No trials matched.'); end
OUT.trials = vertcat(rows{:});

% free-fall cross-check on v0
OUT.trials.v0_ff     = sqrt(2*980*OUT.trials.h_true_mm/10);
OUT.trials.v0_ratio  = abs(OUT.trials.v0_cm_s) ./ OUT.trials.v0_ff;
OUT.trials.midEvent  = OUT.trials.nPreImpact < opt.MinPre;

% ------------------------------------------------------------- 2. summary
fprintf('\n--- Q1/Q2: is the trigger inside the record, and where? ---\n');
fprintf('  %-16s %4s  %-14s %-16s %s\n','condition','n','trigger in range','median nPreImpact','midEvent');
srows = {};
for c = unique(OUT.trials.condition,'stable')'
    m = OUT.trials.condition == c;
    fprintf('  %-16s %4d  %5d / %-6d   %6.0f            %d\n', c, sum(m), ...
        sum(OUT.trials.trigger_in_range(m)), sum(m), ...
        median(OUT.trials.nPreImpact(m)), sum(OUT.trials.midEvent(m)));
end

fprintf('\n--- Q3: EMPIRICAL impactDistPx from the peak-speed frame ---\n');
fprintf('  %-16s %4s  %-32s %-18s %s\n','condition','n','empirical IDP (px)','delta vs used','v0/v_ff');
for c = unique(OUT.trials.condition,'stable')'
    m = OUT.trials.condition == c;
    e = OUT.trials.empirical_IDP_px(m);
    fprintf('  %-16s %4d  median %7.1f  IQR [%7.1f,%7.1f]  %+8.1f px      %.2f\n', ...
        c, sum(m), median(e), prctile(e,25), prctile(e,75), ...
        median(OUT.trials.delta_px(m)), median(OUT.trials.v0_ratio(m)));
    srows{end+1} = table(c, sum(m), median(e), prctile(e,25), prctile(e,75), std(e), ...
        median(OUT.trials.delta_px(m)), median(OUT.trials.v0_ratio(m)), ...
        sum(OUT.trials.midEvent(m)), sum(~OUT.trials.trigger_in_range(m)), ...
        'VariableNames',{'condition','n','emp_IDP_median','emp_IDP_q25','emp_IDP_q75', ...
            'emp_IDP_sd','delta_median_px','v0ratio_median','n_midEvent','n_outOfRange'}); %#ok<AGROW>
end
OUT.summary = vertcat(srows{:});

iG = OUT.summary.condition == "GB/full";
if any(iG)
    fprintf(['\n  VALIDATION: GB/full empirical = %.0f px vs used %.0f px (delta %+.0f).\n' ...
             '  If that delta is small, the method works and the other rows are\n' ...
             '  trustworthy. If it is large, stop -- the method is not valid here.\n'], ...
        OUT.summary.emp_IDP_median(iG), median(OUT.trials.impactDistPx_used), ...
        OUT.summary.delta_median_px(iG));
end

% ------------------------------------------------- 3. drift with drop height
fprintf('\n--- Does the empirical trigger drift with drop height? ---\n');
fprintf('  (a single constant is only defensible if it does not)\n');
for c = unique(OUT.trials.condition,'stable')'
    m = OUT.trials.condition == c;
    h = OUT.trials.h_true_mm(m); e = OUT.trials.empirical_IDP_px(m);
    if numel(unique(h)) < 3, continue; end
    r = corr(h, e);
    fprintf('  %-16s corr(h, empirical IDP) = %+.2f  %s\n', c, r, ...
        local_tern(abs(r) > 0.5, '<-- drifts; no single value fits', ''));
end

% ------------------------------------------------------------- 4. figures
col = containers.Map({'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {[0 .45 .74],[.30 .75 .93],[.85 .33 .10],[.64 .08 .18]});
figure('Color','w','Position',[40 60 1500 440]);
tl = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
for c = unique(OUT.trials.condition,'stable')'
    m = OUT.trials.condition == c; cc = col(char(c));
    plot(ax1, OUT.trials.h_true_mm(m), OUT.trials.empirical_IDP_px(m), 'o', ...
        'MarkerSize',5,'MarkerFaceColor',cc,'MarkerEdgeColor','none','DisplayName',c);
end
yline(ax1, median(OUT.trials.impactDistPx_used), 'k--', 'used');
yline(ax1, -290, 'k:', 'old override -290');
xlabel(ax1,'true drop height (mm)'); ylabel(ax1,'empirical impactDistPx (px)');
title(ax1,'What the data say the trigger should be');
legend(ax1,'Location','best','Box','off','FontSize',8);

ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
for c = unique(OUT.trials.condition,'stable')'
    m = OUT.trials.condition == c; cc = col(char(c));
    plot(ax2, OUT.trials.h_true_mm(m), max(OUT.trials.nPreImpact(m),0.5), 'o', ...
        'MarkerSize',5,'MarkerFaceColor',cc,'MarkerEdgeColor','none','DisplayName',c);
end
yline(ax2, opt.MinPre, 'r--', 'starts mid-event');
set(ax2,'YScale','log');
xlabel(ax2,'true drop height (mm)'); ylabel(ax2,'frames before impact\_index');
title(ax2,'Is the impact inside the record?');

ax3 = nexttile(tl); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
for c = unique(OUT.trials.condition,'stable')'
    m = OUT.trials.condition == c; cc = col(char(c));
    plot(ax3, OUT.trials.delta_px(m), OUT.trials.v0_ratio(m), 'o', ...
        'MarkerSize',5,'MarkerFaceColor',cc,'MarkerEdgeColor','none','DisplayName',c);
end
yline(ax3, 1, 'k--'); xline(ax3, 0, 'k--');
xlabel(ax3,'empirical - used  (px)'); ylabel(ax3,'v_0 / v_{freefall}');
title(ax3,'Does the trigger error explain the v_0 error?');
title(tl,'impactDistPx diagnostic -- read-only, nothing rewritten','FontWeight','bold');

fprintf('\nNothing written. OUT.trials and OUT.summary returned.\n\n');
end

% ------------------------------------------------------------------ helpers
function h = local_h(meta)
if isfield(meta,'dropHeight_true_mm') && isfinite(meta.dropHeight_true_mm)
    h = meta.dropHeight_true_mm;
else
    h = meta.dropHeight_mm;
end
end
function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
