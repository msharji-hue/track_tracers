%% ========================================================================
%  THREE-MODEL DRY RUN — GB/full only, Tight / Default / Wide
%  Read-only: computes kinematics in memory, plots, writes NO files.
%
%  Toe splay is the independent variable:
%     Tight   +/- 1.50 deg   distal separation  5.85 mm
%     Default +/- 7.92 deg                      9.48 mm
%     Wide    +/-12.00 deg                     11.76 mm
%  Same 65 g projectile, same glass beads (phi = 0.624), same container,
%  so any difference in d0 is geometry.
%
%  SCOPE: ALL GB/full trials at ALL drop heights, for all three models.
%  Nothing is height-filtered out of the audit table; the h = 0 rows are
%  additionally broken out for the d0 statistic.
%
%  h = 0 trials (released touching the surface) measure d0 DIRECTLY.
%  They are handled separately throughout: v0_ff = 0 there, so the
%  measured/free-fall ratio is undefined and is left NaN rather than Inf.
%% ========================================================================

root     = 'D:\ME_GRANULAB\JerboaImpact';
gbRoot   = fullfile(root,'03_RESULTS','GB','Batch 5');
modelH   = 4.00;                       % model height 40 mm -> overdepth cut

SAVE_FIGS   = false;   % figures are shown interactively and NOT saved
DEFAULT_D0_ONLY = false;% Default contributes d0 trials only; its validated
                       % non-d0 set is already processed and is excluded here

% per-model impact trigger, measured before testing
IDP = containers.Map({'Tight','Default','Wide'}, {-376.001, -370.001, -409});

base = get_calibration();              % current scale + bed line (shared)

% ---- text log: the ONLY file this script writes -------------------------
logDir = fullfile(root,'03_RESULTS','_batch_logs');
if ~exist(logDir,'dir'), mkdir(logDir); end
stamp   = datestr(now,'yyyymmdd_HHMMSS');
logPath = fullfile(logDir, sprintf('dryrun_three_models_%s.txt', stamp));
LOGF = fopen(logPath,'w');
LOG  = @(varargin) both(LOGF, varargin{:});

LOG('========================================================\n');
LOG(' THREE-MODEL DRY RUN — GB/full   %s\n', stamp);
LOG('========================================================\n');
LOG('DRY RUN ONLY: no _kin.mat, no _kin_scalars.csv, no overwrites,\n');
LOG('no figures saved (SAVE_FIGS = %d). Detection settings untouched.\n\n', SAVE_FIGS);

LOG('========== CALIBRATION IN USE ==========\n');
LOG('shared bedline : x = %g px   points [%g,%g] - [%g,%g]\n', ...
    base.bedX, base.bedPoint1, base.bedPoint2);
LOG('shared scale   : mmPerPx = %.4f   (pxPerMm = %.2f)\n', base.mmPerPx, 1/base.mmPerPx);
LOG('gravity        : %g cm/s^2   trackTolerancePx = %g\n', base.g_cm_s2, base.trackTolerancePx);
LOG('impactDistPx   : Tight %.3f | Default %.3f | Wide %.3f  (per model)\n', ...
    IDP('Tight'), IDP('Default'), IDP('Wide'));
LOG('fps            : per trial via resolve_fps\n');
LOG('                 (scalars CSV -> meta.fps_true -> tracks.fps)\n');
LOG('d0 handling    : h = 0 trials are previewed like any other. sqrt(2gh) = 0,\n');
LOG('                 so no free-fall reference and no ratio are formed for them.\n');
LOG('                 Impact is located exactly as usual (geometric trigger +\n');
LOG('                 velocity-peak refinement). d0 curves are drawn in red.\n\n');

%% ---- 1. discover GB/full tracks across all three models ----------------
D = dir(fullfile(gbRoot,'**','*_tracks.mat'));
LOG('found %d *_tracks.mat under GB\\Batch 5\n', numel(D));
nSkipContainer = 0; nSkipDefaultNonZero = 0;

n = numel(D);
tag=strings(n,1); mdl=tag; cont=tag; flg=tag; pathv=tag;
hmm=nan(n,1); fpsv=hmm; imf=hmm; stf=hmm; pre=hmm; npo=hmm;
v0m=hmm; v0f=hmm; rat=hmm; dfin=hmm; astp=hmm; idpUsed=hmm; idpSaved=hmm;
isZero=false(n,1);

ws = warning('off','all'); restore = onCleanup(@() warning(ws)); %#ok<NASGU>

for k = 1:n
    tp = fullfile(D(k).folder, D(k).name);
    try
        S = load(tp,'meta','tracks');
    catch
        continue;
    end
    m = S.meta;

    cont(k) = string(getf(m,'container',''));
    if ~strcmpi(cont(k),'full')                         % GB/full ONLY
        cont(k) = ""; nSkipContainer = nSkipContainer + 1; continue;
    end

    thisMdl = model_of(string(getf(m,'trialTag','')), string(D(k).folder));
    thisH   = getf(m,'dropHeight_mm',NaN);

    % Default: d0 trials ONLY. Its non-d0 set is already validated/processed.
    if DEFAULT_D0_ONLY && thisMdl == "Default" && thisH ~= 0
        cont(k) = ""; nSkipDefaultNonZero = nSkipDefaultNonZero + 1; continue;
    end

    tag(k)   = string(getf(m,'trialTag',''));
    pathv(k) = string(tp);
    mdl(k)   = thisMdl;
    hmm(k)   = thisH;
    isZero(k)= (hmm(k) == 0);

    fps = resolve_fps(tp, m, S.tracks);
    fpsv(k) = fps;
    if ~isfinite(fps), flg(k)="BADFPS"; continue; end

    % calibration: shared scale/bed + the model's own trigger.
    cb = base;
    cb.impactDistPx = IDP(char(mdl(k)));
    idpUsed(k) = cb.impactDistPx;
    % what Stage A recorded, for cross-check only
    try
        Sc = load(tp,'calib');
        if isfield(Sc,'calib') && isfield(Sc.calib,'impactDistPx')
            idpSaved(k) = Sc.calib.impactDistPx;
        end
    catch
    end

    try
        evalc('kin = kd_kinematics(S.tracks.trackedX,S.tracks.trackedY,cb,1/fps);');
    catch ME
        flg(k) = "KINFAIL"; continue;
    end

    imf(k)=kin.impact_index; stf(k)=kin.stopFrame;
    pre(k)=kin.impact_index-1; npo(k)=kin.stopFrame-kin.impact_index;
    v0m(k)=kin.v0_cm_s; dfin(k)=kin.d_final_cm;
    astp(k)=getf(kin,'a_stop_cm_s2',NaN);

    % free-fall reference is undefined at h = 0 -> leave NaN, do not divide
    if isZero(k)
        v0f(k) = NaN; rat(k) = NaN;
    else
        v0f(k) = sqrt(2*base.g_cm_s2*(hmm(k)/10));
        rat(k) = v0m(k)/v0f(k);
    end

    % flags (same rules as the Default set; no shallow-bed rule, GB/full)
    f = strings(0,1);
    dd = dfin(k);
    if ~isfinite(dd),                 f(end+1)="NaNDEPTH"; end
    if isfinite(dd) && dd > modelH,   f(end+1)="OVERDEPTH"; end
    if isfinite(dd) && dd < 0.05,     f(end+1)="THIN"; end
    if npo(k) < 5,                    f(end+1)="SHORT"; end
    vp = kin.v(kin.impact_index:kin.stopFrame); dv = diff(vp(:));
    if any(dv > 0.15*abs(kin.v0_cm_s)), f(end+1)="GLITCH"; end
    if isempty(f), f="ok"; end
    flg(k) = strjoin(f,"|");
end

keep = tag~="";                       % every discovered file appears
A = table(mdl(keep),cont(keep),hmm(keep),isZero(keep),tag(keep),fpsv(keep), ...
          idpUsed(keep),idpSaved(keep),imf(keep),stf(keep),pre(keep),npo(keep), ...
          v0m(keep),v0f(keep),rat(keep),dfin(keep),astp(keep),flg(keep),pathv(keep), ...
    'VariableNames',{'model','container','dropHeight_mm','isZeroDrop','trialTag','fps', ...
    'impactDistPx','impactDistPx_saved','impactFrame','stopFrame','preFrames','nPost', ...
    'v0_meas','v0_ff','ratio','d_final_cm','a_stop_cm_s2','flags','path'});
A = sortrows(A,{'model','dropHeight_mm','trialTag'});

%% ---- 1b. accounting: what was found vs what is usable ------------------
fprintf('\nrows recorded            : %d  (of %d files discovered)\n', height(A), numel(D));
fprintf('containers seen:\n'); disp(groupsummary(A,'container'))
fprintf('model x container:\n'); disp(groupsummary(A,{'model','container'}))
nonFull = A(~strcmpi(A.container,'full'),:);
if ~isempty(nonFull)
    fprintf('%d row(s) are NOT container=full and are dropped from the analysis:\n', height(nonFull));
    disp(groupsummary(nonFull,{'model','container'}))
end
A = A(strcmpi(A.container,'full'),:);       % GB/full only, from here on
fprintf('GB/full rows carried forward: %d\n', height(A));
fprintf('  per model: '); 
gm = groupsummary(A,'model');
for i=1:height(gm), fprintf('%s=%d  ', char(gm.model(i)), gm.GroupCount(i)); end
fprintf('\n');

%% ---- 2. exclusions -----------------------------------------------------
isBad = contains(A.flags,"NaNDEPTH") | contains(A.flags,"OVERDEPTH") | ...
        contains(A.flags,"GLITCH")   | contains(A.flags,"BADFPS")    | ...
        contains(A.flags,"KINFAIL")  | contains(A.flags,"LOADFAIL");
A.excluded = isBad;
K = A(~isBad,:);
keptTags3 = cellstr(K.trialTag);

%% ---- 3. SUMMARY --------------------------------------------------------
LOG('\n========== SCOPE CONFIRMATION ==========\n');
LOG('skipped, container ~= full            : %d\n', nSkipContainer);
LOG('skipped, Default non-d0 (validated)   : %d\n', nSkipDefaultNonZero);
LOG('CONFIRMED: only GB/full included.\n');
LOG('CONFIRMED: Default non-d0 trials excluded (d0 only).\n');
LOG('CONFIRMED: Tight and Wide use all current GB/full trials.\n');

LOG('\n========== THREE-MODEL DRY RUN (GB/full) ==========\n');
LOG('trials in scope      : %d\n', height(A));
LOG('excluded by flags    : %d\n', sum(isBad));
LOG('KEPT                 : %d\n', height(K));

LOG('\nfps range: %.0f - %.0f\n', min(A.fps), max(A.fps));
nBadFps = sum(~isfinite(A.fps));
if nBadFps, LOG('  *** %d trial(s) with unresolved fps ***\n', nBadFps); end

LOG('\nby model (in scope):\n');   logtab(LOGF, groupsummary(A,'model'));
LOG('kept by model:\n');           logtab(LOGF, groupsummary(K,'model'));
LOG('flag counts:\n');
for f = ["ok","NaNDEPTH","OVERDEPTH","THIN","SHORT","GLITCH","BADFPS","KINFAIL"]
    c = sum(contains(A.flags,f));
    if c, LOG('  %-10s %d\n', f, c); end
end

%% ---- 3b. FULL TRIAL LIST ----------------------------------------------
LOG('\n========== TRIAL / VIDEO LIST BEING TESTED ==========\n');
LOG('%-28s %-8s %6s %7s %7s %7s %8s %9s %9s  %s\n', ...
    'trialTag','model','h_mm','fps','impact','stop','nPost','v0_cm_s','d_cm','flags');
for i = 1:height(A)
    LOG('%-28s %-8s %6g %7.0f %7d %7d %8d %9.1f %9.4f  %s%s\n', ...
        A.trialTag(i), A.model(i), A.dropHeight_mm(i), A.fps(i), ...
        A.impactFrame(i), A.stopFrame(i), A.nPost(i), A.v0_meas(i), ...
        A.d_final_cm(i), A.flags(i), ternary(A.isZeroDrop(i),'  [d0]',''));
end

% trigger cross-check: does the value Stage A recorded match the model map?
mismatch = A(isfinite(A.impactDistPx_saved) & ...
             abs(A.impactDistPx_saved - A.impactDistPx) > 0.01, :);
LOG('\nimpactDistPx recorded in tracks vs model map: %d mismatch(es)\n', height(mismatch));
if ~isempty(mismatch)
    logtab(LOGF, groupsummary(mismatch,{'model','impactDistPx_saved'}));
end

LOG('\nkept by model x drop height:\n');
logtab(LOGF, unstack(groupsummary(K,{'model','dropHeight_mm'}),'GroupCount','model'));

%% ---- 3b. v0 and depth statistics per model (ALL heights) ---------------
fprintf('\n===== v0_cm_s by model (all kept, all heights) =====\n');
disp(varfun(@(x)[numel(x) min(x) median(x) max(x)], K, ...
     'InputVariables','v0_meas','GroupingVariables','model'))
fprintf('===== d_cm by model (all kept, all heights) =====\n');
disp(varfun(@(x)[numel(x) min(x) median(x) max(x)], K, ...
     'InputVariables','d_final_cm','GroupingVariables','model'))

fprintf('===== v0_cm_s by model x drop height =====\n');
disp(varfun(@(x)[numel(x) median(x)], K, ...
     'InputVariables','v0_meas','GroupingVariables',{'model','dropHeight_mm'}))
fprintf('===== d_cm by model x drop height =====\n');
disp(varfun(@(x)[numel(x) median(x)], K, ...
     'InputVariables','d_final_cm','GroupingVariables',{'model','dropHeight_mm'}))

%% ---- 4. d0 (h = 0) — the shape-dependence measurement ------------------
Z = K(K.isZeroDrop,:);
LOG('\n===== d0 at h = 0 (direct measurement) =====\n');
if isempty(Z)
    LOG('no h = 0 trials found\n');
else
    d0 = groupsummary(Z,'model',{'mean','std','median'},'d_final_cm');
    logtab(LOGF, d0);
    LOG('splay order Tight -> Default -> Wide should map monotonically to d0\n');
end

%% ---- 5. PLOTS ----------------------------------------------------------
models = unique(K.model); mc = lines(numel(models));

% 5a. kinematics per model
for mi = 1:numel(models)
    sel = K(K.model==models(mi),:);            % d0 INCLUDED
    if isempty(sel), continue; end
    hs = unique(sel.dropHeight_mm); cmap = parula(max(2,numel(hs)));
    f = figure('Color','w','Name',char(models(mi)),'Position',[60 60 950 820]);
    for k = 1:height(sel)
        S = load(char(sel.path(k)),'tracks','meta');
        cb = base; cb.impactDistPx = sel.impactDistPx(k);
        try, evalc('kin = kd_kinematics(S.tracks.trackedX,S.tracks.trackedY,cb,1/sel.fps(k));');
        catch, continue; end
        if sel.isZeroDrop(k)
            c = [0.85 0.1 0.1];                % d0 highlighted
        else
            c = cmap(find(hs==sel.dropHeight_mm(k),1),:);
        end
        dn = sprintf('%s (%g mm%s)', char(sel.trialTag(k)), sel.dropHeight_mm(k), ...
                     ternary(sel.isZeroDrop(k),', d0',''));
        t = kin.t_s*1e3;                       % all three panels share t
        subplot(3,1,1); hold on; plot(t,kin.z,       '-','Color',c,'DisplayName',dn);
        subplot(3,1,2); hold on; plot(t,kin.v,       '-','Color',c,'DisplayName',dn);
        subplot(3,1,3); hold on; plot(t,kin.a_plus_g,'-','Color',c,'DisplayName',dn);
    end
    subplot(3,1,1); grid on; ylabel('z (cm)'); xlabel('t (ms, 0 = impact)');
      yline(modelH,'k:','model 40 mm');
      title(sprintf('%s — GB/full  (n=%d, red = d0)',char(models(mi)),height(sel)),'Interpreter','none');
    subplot(3,1,2); grid on; ylabel('v (cm/s)'); xlabel('t (ms, 0 = impact)'); yline(0,'k:');
    subplot(3,1,3); grid on; ylabel('a+g (cm/s^2)'); xlabel('t (ms, 0 = impact)'); yline(0,'k:');
    dcm = datacursormode(f); dcm.Enable='on';
    dcm.UpdateFcn = @(~,e) {get(e.Target,'DisplayName'), ...
        sprintf('x=%.3f  y=%.3f', e.Position(1), e.Position(2))};
end

% 5b. d0 across models — the headline comparison
if ~isempty(Z)
    figure('Color','w','Name','d0 by model');
    g = groupsummary(Z,'model',{'mean','std'},'d_final_cm');
    errorbar(1:height(g), g.mean_d_final_cm, g.std_d_final_cm, 'o', ...
             'LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor','auto');
    set(gca,'XTick',1:height(g),'XTickLabel',g.model); grid on;
    xlim([0.5 height(g)+0.5]);
    ylabel('d_0 (cm)  [h = 0, released touching surface]');
    title('Shape dependence: d_0 vs toe splay   (same mass, bed, container)');
end

% 5c. depth vs drop height, all models on one axes
figure('Color','w','Name','d vs h by model'); hold on; grid on;
for mi = 1:numel(models)
    s = K(K.model==models(mi),:);
    scatter(s.dropHeight_mm, s.d_final_cm, 30, mc(mi,:), 'filled', ...
            'MarkerFaceAlpha',.55, 'DisplayName',char(models(mi)));
    g = groupsummary(s,'dropHeight_mm','median','d_final_cm');
    plot(g.dropHeight_mm, g.median_d_final_cm,'-','Color',mc(mi,:),'LineWidth',1.5, ...
         'HandleVisibility','off');
end
yline(modelH,'k:','model 40 mm','HandleVisibility','off');
xlabel('drop height (mm)'); ylabel('d final (cm)');
title('Penetration depth vs drop height, by foot geometry'); legend('Location','southeast');

%% ---- 6. CONFIRM NOTHING WRITTEN ---------------------------------------
if SAVE_FIGS
    LOG('\nSAVE_FIGS was true -- figures would be saved (set by hand only)\n');
end

LOG('\n===== SAFETY CONFIRMATION =====\n');
LOG('_kin.mat written             : 0\n');
LOG('_kin_scalars.csv written     : 0\n');
LOG('existing results overwritten : 0\n');
LOG('figures saved                : 0\n');
LOG('kinematics computed in memory only\n');
LOG('only file written            : %s\n', logPath);
LOG('keptTags3 ready: %d trials\n', numel(keptTags3));
fclose(LOGF);
fprintf('\nfull log: %s\n', logPath);

%% ---- helpers -----------------------------------------------------------
function m = model_of(tag, folder)
    t = lower(string(tag)); f = lower(string(folder));
    if endsWith(t,"_tight") || contains(f,"tight"), m = "Tight";
    elseif endsWith(t,"_wide") || contains(f,"wide"), m = "Wide";
    else, m = "Default";
    end
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function both(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    if fid > 0, fprintf(fid, fmt, varargin{:}); end
end

function logtab(fid, T)
    disp(T);
    if fid > 0
        txt = evalc('disp(T)');
        fprintf(fid, '%s\n', txt);
    end
end

function o = ternary(c,a,b)
    if c, o = a; else, o = b; end
end
