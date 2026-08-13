function [FIT, DATA] = depth_velocity_scaling(root, varargin)
%DEPTH_VELOCITY_SCALING  L0 baseline: penetration depth vs impact velocity.
%
%   Tests the fixed-geometry reduction of Newhall & Durian (2003) Eq.(4).
%   For a single foot the mass m and area A are constant, so
%
%       d = K * v0^(2/3)        (equivalently  d = K * (v0^2)^(1/3) )
%
%   This script tests ONLY that. No force-law fitting, no k or C, no A(z),
%   no area normalisation, no pooling across substrates.
%
%   Geometry: default / middle jerboa foot, side-toe splay +/-7.92 deg.
%   One geometry only, so this is NOT a shape-dependence comparison; it is
%   the baseline the other two feet will later be measured against.
%
%   INPUT SET
%   Kinematics exist on disk only for the final kept trials (earlier outputs
%   were deleted in the pipeline rewrite), so the glob over *_kin_scalars.csv
%   IS the kept set. No exclusion table is read. Expected counts are asserted
%   instead, per condition and in total.
%
%   CONVENTIONS (per pipeline changelog)
%     v0_cm_s              measured impact velocity; the reported predictor
%     dropHeight_true_mm   physical height; use for height-based grouping
%     dropHeight_mm        labelled value; IDENTIFIER ONLY (GB/shallow labels
%                          are reversed on disk)
%     d_final_cm           depth at t_stop, i.e. maximum penetration
%     depth cutoff         4.00 cm (visible model height)
%     shallow bed          2.50 cm (GB/shallow container floor)
%
%   USAGE
%       root = 'D:\ME_GRANULAB\JerboaImpact';
%       [FIT, DATA] = depth_velocity_scaling(root);
%
%   OPTIONS (name/value)
%       'OutDir'        output folder (default <root>/03_RESULTS/_batch_logs)
%       'ExpectKept'    expected total (default 269)
%       'DepthCut'      model-height cutoff, cm (default 4.00, drawn only)
%       'ShallowBed'    GB/shallow bed depth, cm (default 2.50, drawn only)
%
%   Base MATLAB only; no toolbox required.

% ------------------------------------------------------------------ options
if nargin < 1 || isempty(root)
    error('depth_velocity_scaling:noRoot', 'Supply the JerboaImpact data root.');
end

opt.OutDir     = fullfile(root, '03_RESULTS', '_batch_logs');
opt.ExpectKept = 269;
opt.DepthCut   = 4.00;
opt.ShallowBed = 2.50;

for i = 1:2:numel(varargin)
    f = varargin{i};
    if ~isfield(opt, f)
        error('depth_velocity_scaling:badOption', 'Unknown option: %s', f);
    end
    opt.(f) = varargin{i+1};
end

FOOT_LABEL = '\pm7.92\circ foot';
G_CM_S2    = 980;                    % matches get_calibration

% expected kept per condition (pipeline changelog, final counts)
expCond = struct('GB_full',154, 'GB_shallow',27, ...
                 'CHIN_as_poured',40, 'CHIN_dense',48);

fprintf('\n=== depth_velocity_scaling | %s ===\n', FOOT_LABEL);

% ------------------------------------------------------- 1. gather scalars
D = dir(fullfile(root, '03_RESULTS', '**', '*_kin_scalars.csv'));
D = D(~[D.isdir]);
if isempty(D)
    error('depth_velocity_scaling:noFiles', ...
          'No *_kin_scalars.csv found under %s', fullfile(root,'03_RESULTS'));
end

n = numel(D);
trialTag = strings(n,1);   condition = strings(n,1);
hLab_mm  = nan(n,1);       hTrue_mm  = nan(n,1);   hCorr = false(n,1);
phi      = nan(n,1);       v0_cm_s   = nan(n,1);
d_final  = nan(n,1);       t_stop    = nan(n,1);
impactFr = nan(n,1);       stopFr    = nan(n,1);

for i = 1:n
    f = fullfile(D(i).folder, D(i).name);
    try
        T = readtable(f);
    catch ME
        warning('depth_velocity_scaling:readFail', ...
                'Skipping %s (%s)', D(i).name, ME.message);
        continue
    end
    if height(T) < 1, continue; end

    trialTag(i)  = string(erase(D(i).name, '_kin_scalars.csv'));
    condition(i) = local_getstr(T, 'condition');
    hLab_mm(i)   = local_getnum(T, 'dropHeight_mm');
    hTrue_mm(i)  = local_getnum(T, 'dropHeight_true_mm');
    hCorr(i)     = local_getnum(T, 'heightCorrected') == 1;
    phi(i)       = local_getnum(T, 'phi');
    v0_cm_s(i)   = local_getnum(T, 'v0_cm_s');
    d_final(i)   = local_getnum(T, 'd_final_cm');
    t_stop(i)    = local_getnum(T, 't_stop_s');
    impactFr(i)  = local_getnum(T, 'impact_frame');
    stopFr(i)    = local_getnum(T, 'stop_frame');
end

valid = trialTag ~= "";
fprintf('kin_scalars files found : %d  (expected %d)\n', sum(valid), opt.ExpectKept);
if sum(valid) ~= opt.ExpectKept
    warning('depth_velocity_scaling:countMismatch', ...
            ['Found %d trials, expected %d. Kinematics should exist only for ' ...
             'the final kept set.'], sum(valid), opt.ExpectKept);
end

% dropHeight_true_mm is required: GB/shallow labels are reversed on disk.
if ~any(isfinite(hTrue_mm))
    error('depth_velocity_scaling:noTrueHeight', ...
          ['dropHeight_true_mm not found in any _kin_scalars.csv. These files ' ...
           'predate the height correction; re-run Stage B before analysing.']);
end

% ------------------------------------------------------- 2. assemble table
depth_cm = d_final;
usable   = valid & isfinite(v0_cm_s) & isfinite(depth_cm) & ...
           v0_cm_s > 0 & depth_cm > 0;
nDrop    = sum(valid) - sum(usable);
if nDrop > 0
    fprintf('dropped (NaN or <=0)    : %d\n', nDrop);
end

DATA = table(trialTag(usable), condition(usable), hLab_mm(usable), ...
             hTrue_mm(usable), hCorr(usable), phi(usable), ...
             abs(v0_cm_s(usable)), depth_cm(usable), t_stop(usable), ...
             impactFr(usable), stopFr(usable), ...
    'VariableNames', {'trialTag','condition','dropHeight_mm', ...
                      'dropHeight_true_mm','heightCorrected','phi', ...
                      'v0_cm_s','d_final_cm','t_stop_s', ...
                      'impact_frame','stop_frame'});

% Free-fall cross-check, built on the TRUE height. Reported only.
DATA.v0_ff_true = sqrt(2 * G_CM_S2 * DATA.dropHeight_true_mm / 10);
DATA.v0_ratio   = DATA.v0_cm_s ./ DATA.v0_ff_true;

fprintf('usable for fitting      : %d\n', height(DATA));
fprintf('height-corrected rows   : %d  (GB/shallow)\n\n', sum(DATA.heightCorrected));

% ------------------------------------------------------- 3. per-condition fit
condOrder = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
present   = unique(DATA.condition, 'stable');
conds     = [condOrder(ismember(condOrder, present)), ...
             reshape(present(~ismember(present, condOrder)), 1, [])];
nc = numel(conds);
rows = cell(nc,1);

for c = 1:nc
    m = DATA.condition == conds(c);
    v = DATA.v0_cm_s(m);   d = DATA.d_final_cm(m);   r = DATA.v0_ratio(m);

    A = local_ols(log10(v), log10(d));      % free exponent
    B = local_ols(v.^(2/3), d);             % exponent fixed at 2/3

    rows{c} = table(conds(c), numel(v), local_expect(expCond, conds(c)), ...
        min(v), max(v), min(d), max(d), ...
        A.b(2), A.ci(2,1), A.ci(2,2), ...
        (A.ci(2,1) <= 2/3) && (A.ci(2,2) >= 2/3), ...
        A.b(1), A.R2, ...
        B.b(2), B.b(1), B.ci(1,1), B.ci(1,2), B.R2, B.rmse, ...
        median(r), min(r), max(r), sum(r < 0.7 | r > 1.3), ...
        'VariableNames', {'condition','n','n_expected', ...
            'v0_min','v0_max','d_min','d_max', ...
            'exponent','exp_lo95','exp_hi95','exp_includes_2_3', ...
            'log10K','R2_loglog', ...
            'slope_a','intercept_d0','d0_lo95','d0_hi95','R2_linear','rmse_cm', ...
            'v0ratio_med','v0ratio_min','v0ratio_max','n_v0ratio_off'});
end
FIT = vertcat(rows{:});

% ------------------------------------------------------------- 4. console
fprintf('--- counts ---\n');
for c = 1:nc
    fprintf('  %-16s n = %3d   expected %3d %s\n', FIT.condition(c), ...
        FIT.n(c), FIT.n_expected(c), ...
        local_tern(FIT.n(c)==FIT.n_expected(c), '', '  <-- MISMATCH'));
end

fprintf('\n--- Fit A: free exponent, d ~ v0^n   (target n = 2/3 = 0.667) ---\n');
for c = 1:nc
    fprintf('  %-16s n = %6.3f  [%6.3f, %6.3f]  R2 = %.3f %s\n', ...
        FIT.condition(c), FIT.exponent(c), FIT.exp_lo95(c), FIT.exp_hi95(c), ...
        FIT.R2_loglog(c), local_tern(FIT.exp_includes_2_3(c), ' <- 2/3 in CI', ''));
end

fprintf('\n--- Fit B: d = a*v0^(2/3) + d0 ---\n');
for c = 1:nc
    fprintf('  %-16s a = %7.4f   d0 = %6.3f cm [%6.3f, %6.3f]   R2 = %.3f   rmse = %.3f cm\n', ...
        FIT.condition(c), FIT.slope_a(c), FIT.intercept_d0(c), ...
        FIT.d0_lo95(c), FIT.d0_hi95(c), FIT.R2_linear(c), FIT.rmse_cm(c));
end

fprintf('\n--- v0 vs free-fall from TRUE height (cross-check, nothing excluded) ---\n');
for c = 1:nc
    fprintf('  %-16s median %.2f   range [%.2f, %.2f]   outside 0.7-1.3: %3d/%3d\n', ...
        FIT.condition(c), FIT.v0ratio_med(c), FIT.v0ratio_min(c), ...
        FIT.v0ratio_max(c), FIT.n_v0ratio_off(c), FIT.n(c));
end
fprintf('\n');

% -------------------------------------------------------------- 5. figure
col = containers.Map( ...
    {'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {[0.00 0.45 0.74],[0.30 0.75 0.93],[0.85 0.33 0.10],[0.64 0.08 0.18]});
mk  = containers.Map( ...
    {'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {'o','s','^','d'});

fig = figure('Color','w','Position',[80 80 1180 480]);
tl  = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

% -- Panel A: log-log, free exponent
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log','YScale','log');
hA = gobjects(nc,1); lblA = strings(nc,1);
for c = 1:nc
    key = char(conds(c));  m = DATA.condition == conds(c);
    v = DATA.v0_cm_s(m);   d = DATA.d_final_cm(m);
    cc = local_mapget(col, key, [0.4 0.4 0.4]);
    hA(c) = plot(ax1, v, d, local_mapget(mk,key,'o'), 'LineStyle','none', ...
        'MarkerSize',5, 'MarkerFaceColor',cc, 'MarkerEdgeColor','none');
    vv = linspace(min(v), max(v), 50);
    plot(ax1, vv, 10^FIT.log10K(c) * vv.^FIT.exponent(c), '-', ...
        'Color',cc, 'LineWidth',1.6);
    lblA(c) = sprintf('%s  (n = %.3f)', conds(c), FIT.exponent(c));
end
xg = exp(mean(log(DATA.v0_cm_s)));  yg = exp(mean(log(DATA.d_final_cm)));
xr = [min(DATA.v0_cm_s) max(DATA.v0_cm_s)];
hRef = plot(ax1, xr, yg*(xr/xg).^(2/3), 'k--', 'LineWidth',1.4);
xlabel(ax1,'impact velocity  v_0  (cm s^{-1})');
ylabel(ax1,'penetration depth  d\_final  (cm)');
title(ax1,'A.  free exponent:  d \propto v_0^{\itn}');
legend(ax1,[hA;hRef],[cellstr(lblA);{'slope 2/3 guide'}], ...
    'Location','southeast','Box','off','FontSize',8);

% -- Panel B: d vs v0^(2/3), linear axes
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
hB = gobjects(nc,1); lblB = strings(nc,1);
for c = 1:nc
    key = char(conds(c));  m = DATA.condition == conds(c);
    x = DATA.v0_cm_s(m).^(2/3);  d = DATA.d_final_cm(m);
    cc = local_mapget(col, key, [0.4 0.4 0.4]);
    hB(c) = plot(ax2, x, d, local_mapget(mk,key,'o'), 'LineStyle','none', ...
        'MarkerSize',5, 'MarkerFaceColor',cc, 'MarkerEdgeColor','none');
    xx = linspace(min(x), max(x), 50);
    plot(ax2, xx, FIT.slope_a(c)*xx + FIT.intercept_d0(c), '-', ...
        'Color',cc, 'LineWidth',1.6);
    lblB(c) = sprintf('%s  (d_0 = %.2f cm, R^2 = %.2f)', ...
        conds(c), FIT.intercept_d0(c), FIT.R2_linear(c));
end
yline(ax2, opt.DepthCut, '--', sprintf('model height %.2f cm', opt.DepthCut), ...
    'Color',[0.45 0.45 0.45], 'LabelHorizontalAlignment','left', 'FontSize',8);
if any(DATA.condition == "GB/shallow")
    yline(ax2, opt.ShallowBed, ':', sprintf('shallow bed %.2f cm', opt.ShallowBed), ...
        'Color',local_mapget(col,'GB/shallow',[0.3 0.75 0.93]), ...
        'LabelHorizontalAlignment','right', 'FontSize',8);
end
xlabel(ax2,'v_0^{2/3}  ((cm s^{-1})^{2/3})');
ylabel(ax2,'penetration depth  d\_final  (cm)');
title(ax2,'B.  exponent fixed:  d = a\cdotv_0^{2/3} + d_0');
legend(ax2,hB,cellstr(lblB),'Location','northwest','Box','off','FontSize',8);

title(tl, sprintf('Penetration depth vs impact velocity  |  %s  |  %d trials  |  depth = d\\_final', ...
    FOOT_LABEL, height(DATA)), 'FontWeight','bold');

% ---------------------------------------------------------------- 6. save
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));

pngPath  = fullfile(opt.OutDir, sprintf('depth_velocity_scaling_%s.png', stamp));
figPath  = fullfile(opt.OutDir, sprintf('depth_velocity_scaling_%s.fig', stamp));
fitPath  = fullfile(opt.OutDir, sprintf('depth_velocity_fits_%s.csv',   stamp));
dataPath = fullfile(opt.OutDir, sprintf('depth_velocity_data_%s.csv',   stamp));

exportgraphics(fig, pngPath, 'Resolution', 200);
savefig(fig, figPath);
writetable(FIT,  fitPath);
writetable(DATA, dataPath);

fprintf('wrote:\n  %s\n  %s\n  %s\n  %s\n\n', pngPath, figPath, fitPath, dataPath);

end % ============================================================== main end


% ------------------------------------------------------------------ helpers

function S = local_ols(x, y)
%LOCAL_OLS  y = b(1) + b(2)*x by ordinary least squares, with 95% CIs.
x = x(:); y = y(:);
n = numel(x);  X = [ones(n,1), x];
b = X \ y;  res = y - X*b;  dof = n - 2;
if dof < 1
    S.b = b; S.ci = nan(2,2); S.R2 = NaN; S.rmse = NaN; S.n = n;
    return
end
s2 = (res'*res)/dof;
se = sqrt(diag(s2 * ((X'*X) \ eye(2))));
tc = local_tcrit95(dof);
S.b    = b;
S.ci   = [b - tc*se, b + tc*se];
S.R2   = 1 - (res'*res)/sum((y - mean(y)).^2);
S.rmse = sqrt(s2);
S.n    = n;
end

function t = local_tcrit95(dof)
try, t = tinv(0.975, dof); catch, t = 1.96; end
end

function v = local_getnum(T, name)
v = NaN;
if ismember(name, T.Properties.VariableNames)
    x = T.(name)(1);
    if isnumeric(x) || islogical(x), v = double(x);
    else, v = str2double(string(x));
    end
end
end

function s = local_getstr(T, name)
s = "";
if ismember(name, T.Properties.VariableNames)
    s = strtrim(string(T.(name)(1)));
end
end

function e = local_expect(expCond, condName)
%LOCAL_EXPECT  Expected kept count for a condition string; NaN if unknown.
f = char(strrep(strrep(condName, '/', '_'), '-', '_'));
if isfield(expCond, f), e = expCond.(f); else, e = NaN; end
end

function v = local_mapget(M, key, default)
if isKey(M, key), v = M(key); else, v = default; end
end

function s = local_tern(cond, a, b)
if cond, s = a; else, s = b; end
end
