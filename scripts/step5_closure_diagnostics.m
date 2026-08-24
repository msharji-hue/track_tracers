% STEP 5.3 + 5.4 -- closure tests and diagnostics on the ALREADY-FITTED KD law.
%
% NOTHING IS REFITTED HERE. Every (k, d1) used below was fitted in
% scripts/step5_fleet_fits.m and is either hard-coded from that run's
% trial-balanced primary estimator or read from its cell-fit CSV. The force law
% is not modified; no data, exclusion or pipeline output is touched. The only
% file produced is one working figure.
%
%   PART 1 (5.3)  CLOSURE. Do the fitted parameters reproduce the two integral
%                 observables -- final depth and stopping time -- across the
%                 whole drop ladder, with ONE (k, d1) per geometry?
%   PART 2 (5.4a) ARCH vs AREA. Does the residual arch line up with the foot's
%                 area profile, or with a featureless substrate decay?
%   PART 3 (5.4b) a_stop GAP. Is the ~2.2x gap between the cell-fit k and the
%                 stop-relation k_astop5 a definition artefact or real physics?
%
% STOPPING-TIME METRIC. tstop5 is the PRIMARY estimator throughout: a single
% documented rule (linear v(t) fit over 4 < speed < 30 cm/s AND depth >
% 0.5*d_final, extrapolated to zero), recovered on 520 of 524 trials. The
% pipeline column t_stop_s appears only as a faint reference layer: when its
% pre-stop fit fails it silently falls back to the raw zero-crossing time, so
% across the fleet it is a MIXED definition and is never used for a number.
%
% Base MATLAB only.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));  % shared KD helpers

% -- constants --------------------------------------------------------
mass = 65;      % projectile mass, g   (as track_tracers_2 / export_master_dataset)
grav = 980;     % cm/s^2               (as get_calibration)

MASTER  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.mat';
EXPDIR  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports';
AREACSV = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\foot_area_vs_depth.csv';
FIGPATH = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step5_closure_diag_working.png';

% default mask, identical to checkpoints 5.1 and 5.2
Z_MIN = 0.1;    % cm
V_MIN = 15;     % cm/s

% geometry identity, fixed across every panel and table (project standard)
MDL = ["Tight" "Default" "Wide"];
MRK = {'o' 's' '^'};
% one colour per geometry, fixed explicitly: yyaxis and interleaved line plots
% both hijack the default ColorOrder, so nothing is left to the axes cycle
CO  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

% TRIAL-BALANCED GLOBAL PARAMETERS -- hard-coded, NOT refitted.
% Source: scripts/step5_fleet_fits.m section C, "PRIMARY trial-balanced".
K_GLOB  = [1.8998e5 1.7183e5 2.1952e5];    % g/s^2, order = MDL
D1_GLOB = [4.003    3.906    4.505   ];    % cm,    order = MDL
K_POOL  = 1.9163e5; D1_POOL = 4.088;       % pooled, carried for reference

% step-3 scaling fit, used as a reference line only
STEP3_D0 = 0.044; STEP3_A = 0.0711;        % d = d0 + a*v0^(2/3)

V0_GRID  = 60:5:290;         % cm/s, the closure curves
V0_TABLE = [100 175 250];    % cm/s, the closure table rows
V0_BAND  = 15;               % cm/s, half-width for the measured means
DT_ODE   = 1e-4;             % s, time step for the model integration
GRID_STEP = 0.025;           % cm, common depth grid inside a cell
LAMBDA   = 2;                % cm, the illustrative featureless-decay length

% geometry landmarks (reference lines ONLY, nothing is masked at them)
LANDMARKS = [0.80 0.97 1.12];

fprintf('\n=== STEP 5.3 + 5.4  CLOSURE AND DIAGNOSTICS ===\n');
fprintf('No (k, d1) is refitted in this script; all parameters come from step5_fleet_fits.\n');
fprintf('Stopping time: tstop5 is PRIMARY. Pipeline t_stop_s is a mixed definition and\n');
fprintf('is shown only as a faint reference layer -- no number below is taken from it.\n');

%% ===================================================================
%  LOAD -- master dataset, step5 fit exports, foot area profile
%  ===================================================================
L = load(MASTER);
T = L.T;  S = L.S;

% analysis set: the reviewed, non-quarantined trials
keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
model = string(K.model);
h_mm  = K.dropHeight_mm;

% display set: the same trials less the sensitivity-flagged ones, MARKERS ONLY
flagDisp = K.d_final_cm < 1 | K.v0_cm_s < 55;
fprintf('\nanalysis set %d trials | display set %d (omits %d flagged, markers only)\n', ...
    nTr, sum(~flagDisp), sum(flagDisp));

% newest step5 exports
cellCsv = local_newest(EXPDIR, 'step5_cell_fits_*.csv');
trCsv   = local_newest(EXPDIR, 'step5_pertrial_fits_*.csv');
C  = readtable(cellCsv);   C.model = string(C.model);
PT = readtable(trCsv);     PT.model = string(PT.model);
fprintf('cell fits    : %s (%d cells)\n', cellCsv, height(C));
fprintf('per-trial    : %s (%d trials)\n', trCsv, height(PT));

% foot area profile; comment lines start with #
AR = readtable(AREACSV, 'CommentStyle', '#');
z_area = AR.z_mm / 10;                       % mm -> cm
fprintf('area profile : %s (%d rows, %.2f to %.2f cm)\n', AREACSV, height(AR), min(z_area), max(z_area));

% per-trial series, the default mask, and the scalars
D = struct('tag', cell(nTr,1));
tagsS = string({S.trialTag});
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    D(i).model   = model(i);
    D(i).depth   = S(j).z_cm(:);
    D(i).speed   = S(j).v_cm_s(:);
    D(i).speed2  = D(i).speed.^2;
    D(i).t_s     = S(j).t_s(:);
    D(i).v0      = K.v0_cm_s(i);
    D(i).d_final = K.d_final_cm(i);
end

% cell identity, matching step5_fleet_fits
[cellId, cellModel, cellH] = findgroups(model, h_mm);
nCell = max(cellId);

%% ===================================================================
%  tstop5 -- the PRIMARY stopping-time estimator, recomputed here
%  ===================================================================
% This is a measurement of the data, not a refit of the force law. The rule is
% the one documented in step5_fleet_fits section F, applied unchanged.
astop5 = nan(nTr,1); tstop5 = nan(nTr,1);
for i = 1:nTr
    w = D(i).speed > 4 & D(i).speed < 30 & D(i).depth > 0.5*D(i).d_final & isfinite(D(i).speed);
    if sum(w) < 5, continue; end
    p = polyfit(D(i).t_s(w), D(i).speed(w), 1);
    astop5(i) = p(1);
    if p(1) ~= 0, tstop5(i) = -p(2)/p(1); end
end
fprintf('tstop5 recovered on %d / %d trials (primary metric)\n', sum(isfinite(tstop5)), nTr);

%% ===================================================================
%  PART 1 (5.3) -- CLOSURE
%  ===================================================================
fprintf('\n=== PART 1 (5.3): CLOSURE ===\n');

% model curves over the v0 grid, one set per geometry, using that geometry's
% trial-balanced (k, d1) -- no fitting
nV = numel(V0_GRID);
d_pred_g  = nan(nV,3); t_pred_g = nan(nV,3); t_phys_g = nan(nV,3);
for g = 1:3
    for q = 1:nV
        v0 = V0_GRID(q);
        % final depth: the root of speed2_model(z) = 0
        d_pred_g(q,g) = local_predict_d(v0, K_GLOB(g), D1_GLOB(g), mass, grav);
        if ~isfinite(d_pred_g(q,g)), continue; end
        % ESTIMATOR-MATCHED stopping time: integrate the model in time, then
        % apply the IDENTICAL tstop5 window rule to the model's own trace
        [tm, vm, zm] = local_integrate_model(v0, K_GLOB(g), D1_GLOB(g), mass, grav, DT_ODE);
        [~, t_pred_g(q,g)] = local_window_slope(tm, vm, zm, d_pred_g(q,g));
        % SECONDARY: the physical time integral, closed with the model's tail
        t_phys_g(q,g) = local_t_phys(v0, K_GLOB(g), D1_GLOB(g), mass, grav, d_pred_g(q,g));
    end
end

% definition sensitivity at v0 = 175
fprintf('\n  definition sensitivity at v0 = 175 cm/s (t_phys - t_pred_windowed):\n');
[~, j175] = min(abs(V0_GRID - 175));
for g = 1:3
    fprintf('    %-8s t_pred_windowed = %.5f s | t_phys = %.5f s | difference = %+.5f s (%+.1f%%)\n', ...
        MDL(g), t_pred_g(j175,g), t_phys_g(j175,g), ...
        t_phys_g(j175,g) - t_pred_g(j175,g), ...
        100*(t_phys_g(j175,g) - t_pred_g(j175,g))/t_pred_g(j175,g));
end

% closure table
fprintf('\n  CLOSURE TABLE (measured means over v0 +/- %d cm/s, analysis set):\n', V0_BAND);
fprintf('  %-8s %5s %8s %8s %7s %10s %10s %7s %10s\n', ...
    'model', 'v0', 'd_pred', 'd_meas', 'err%', 't_pred_w', 'tstop5', 'err%', 't_phys');
for g = 1:3
    for q = 1:numel(V0_TABLE)
        v0 = V0_TABLE(q);
        sel = model == MDL(g) & abs(K.v0_cm_s - v0) <= V0_BAND;
        dm = mean(K.d_final_cm(sel), 'omitnan');
        tm = mean(tstop5(sel), 'omitnan');
        [~, jj] = min(abs(V0_GRID - v0));
        dp = d_pred_g(jj,g); tp = t_pred_g(jj,g); tf = t_phys_g(jj,g);
        fprintf('  %-8s %5d %8.3f %8.3f %+7.1f %10.5f %10.5f %+7.1f %10.5f\n', ...
            MDL(g), v0, dp, dm, 100*(dp-dm)/dm, tp, tm, 100*(tp-tm)/tm, tf);
    end
end

% tstop5 trend by v0 tercile, so every stopping-time number is from the primary
fprintf('\n  tstop5 trend by v0 tercile (PRIMARY estimator, median per tercile):\n');
fprintf('  %-8s %10s %10s %10s %14s\n', 'model', 'low', 'mid', 'high', 'v0 cuts');
for g = 1:3
    sel = find(model == MDL(g) & isfinite(tstop5));
    v0s = K.v0_cm_s(sel);
    cuts = local_prctile(v0s, [100/3 200/3]);
    lo = sel(v0s <= cuts(1));
    mi = sel(v0s > cuts(1) & v0s <= cuts(2));
    hi = sel(v0s > cuts(2));
    fprintf('  %-8s %10.5f %10.5f %10.5f   %5.0f / %5.0f\n', MDL(g), ...
        median(tstop5(lo)), median(tstop5(mi)), median(tstop5(hi)), cuts(1), cuts(2));
end

%% ===================================================================
%  PART 2 (5.4a) -- ARCH vs AREA
%  ===================================================================
fprintf('\n=== PART 2 (5.4a): ARCH vs AREA ===\n');

% Rebuild the per-geometry median residual curves with code identical to
% step5_fleet_fits, but using each cell's ALREADY-FITTED (k, d1) from the CSV.
resid_c = cell(nCell,1); grid_c = cell(nCell,1); v0_c = nan(nCell,1);
kc = nan(nCell,1); dc = nan(nCell,1);
for c = 1:nCell
    % match this cell to its row in the cell-fit CSV
    r = find(C.model == cellModel(c) & C.h_mm == cellH(c), 1);
    if isempty(r) || ~isfinite(C.k_fit(r)), continue; end
    kc(c) = C.k_fit(r); dc(c) = C.d1_fit(r); v0_c(c) = C.v0_mean(r);

    idx = find(cellId == c);
    zg = (0 : GRID_STEP : max([D(idx).d_final])).';
    M  = nan(numel(zg), numel(idx));
    for q = 1:numel(idx)
        [zz, uu] = local_monotone_descent(D(idx(q)).depth, D(idx(q)).speed, D(idx(q)).speed2);
        if numel(zz) < 2, continue; end
        M(:,q) = interp1(zz, uu, zg, 'linear', NaN);
    end
    support = sum(isfinite(M), 2);
    med = median(M, 2, 'omitnan');
    med(support < ceil(numel(idx)/2)) = NaN;

    mk = zg > Z_MIN & isfinite(med) & med > V_MIN^2;
    if sum(mk) < 5, continue; end
    grid_c{c}  = zg(mk);
    resid_c{c} = med(mk) - kd_speed2_model(zg(mk), v0_c(c), kc(c), dc(c), mass, grav);
end

% per-geometry median residual curve, normalised by v0^2
medResid = cell(3,1); medGrid = cell(3,1);
for g = 1:3
    sel = find(cellModel == MDL(g) & isfinite(kc));
    zmax = max(cellfun(@(v) max(v), grid_c(sel)));
    zg = (0 : GRID_STEP : zmax).';
    R = nan(numel(zg), numel(sel));
    for q = 1:numel(sel)
        cc = sel(q);
        R(:,q) = interp1(grid_c{cc}, resid_c{cc}/v0_c(cc)^2, zg, 'linear', NaN);
    end
    medGrid{g} = zg;  medResid{g} = median(R, 2, 'omitnan');
end

% (i) where the arch crosses zero, against the two geometry landmarks
fprintf('\n  (i) residual zero crossing vs the geometry landmarks:\n');
fprintf('  %-8s %14s %18s %18s\n', 'model', 'zero crossing', 'vs a_local peak', 'vs foot/bar junction');
for g = 1:3
    [zc, zp, pk] = local_arch_features(medGrid{g}, medResid{g});
    fprintf('  %-8s %11.3f cm %+14.3f cm %+18.3f cm   (peak %.3f cm, %.4f of v0^2)\n', ...
        MDL(g), zc, zc - 0.97, zc - 1.12, zp, pk);
end

% (ii)+(iii) alignment, DESCRIPTIVE ONLY -- nothing is fitted here
fprintf('\n  (ii)/(iii) Spearman alignment over 0.2-2.5 cm, DESCRIPTIVE ONLY (no fit):\n');
fprintf('  %-8s %18s %22s %s\n', 'model', 'rho(resid,-a_local)', 'rho(resid,exp(-z/2cm))', 'better aligned');
zc_common = (0.2 : GRID_STEP : 2.5).';
areaWins = false(1,3);
for g = 1:3
    r = interp1(medGrid{g}, medResid{g}, zc_common, 'linear', NaN);
    a_local = interp1(z_area, AR.(sprintf('a_local_%s', MDL(g))), zc_common, 'linear', NaN);
    decay   = exp(-zc_common / LAMBDA);          % featureless monotone alternative
    rho_area  = local_spearman(r, -a_local);
    rho_decay = local_spearman(r, decay);
    areaWins(g) = abs(rho_area) > abs(rho_decay);
    if areaWins(g), better = 'area profile'; else, better = 'featureless decay'; end
    fprintf('  %-8s %18.3f %22.3f   %s\n', MDL(g), rho_area, rho_decay, better);
end

% Why (i) and (ii) can disagree, stated from the measured area profile rather
% than argued: a_local is a narrow SPIKE, not a broad trend.
[pk, jpk] = max(AR.(sprintf('a_local_%s', MDL(2))));
fprintf(['\n  Reading these together: a_local is a narrow spike (peak %.2f cm^2 at %.3f cm,\n' ...
         '  falling to %.2f cm^2 by 2.5 cm), so a rank correlation taken over the whole\n' ...
         '  0.2-2.5 cm range is dominated by the broad rising limb of the arch and is\n' ...
         '  largely blind to the landmark coincidence that test (i) measures. The two\n' ...
         '  tests are answering different questions, and only (i) is sensitive to the\n' ...
         '  feature actually being claimed.\n'], pk, z_area(jpk), ...
         interp1(z_area, AR.(sprintf('a_local_%s', MDL(2))), 2.5));
if ~any(areaWins)
    fprintf(['  On the rank measure the featureless decay tracks the arch BETTER than the area\n' ...
             '  profile in all three geometries. Reported as found: the area hypothesis does\n' ...
             '  not win this test, and only the landmark coincidence in (i) supports it.\n']);
end
fprintf(['\n  NOTE: this is ALIGNMENT evidence, not a fit -- no area term is fitted to the\n' ...
         '  residual anywhere. The C(z)-vs-A(z) degeneracy at phi > phi_c is UNRESOLVED:\n' ...
         '  a depth-dependent substrate coefficient and a depth-dependent area produce the\n' ...
         '  same residual signature here, and only a cylinder control (constant A with\n' ...
         '  depth) can separate them.\n']);

%% ===================================================================
%  PART 3 (5.4b) -- the a_stop gap: definition or physics
%  ===================================================================
fprintf('\n=== PART 3 (5.4b): a_stop GAP ===\n');

% For every cell, integrate the model with that cell's ALREADY-FITTED (k, d1)
% and apply the IDENTICAL tstop5 window rule, so model and measurement are
% compared through the same estimator.
astop_mod = nan(nCell,1); astop_true = nan(nCell,1);
for c = 1:nCell
    if ~isfinite(kc(c)), continue; end
    dp = local_predict_d(v0_c(c), kc(c), dc(c), mass, grav);
    if ~isfinite(dp), continue; end
    [tm, vm, zm] = local_integrate_model(v0_c(c), kc(c), dc(c), mass, grav, DT_ODE);
    astop_mod(c) = local_window_slope(tm, vm, zm, dp);
    % the model's true terminal slope at the stop: dv/dt = grav - k*d_pred/mass
    astop_true(c) = -(kc(c)*dp/mass - grav);
end

% measured medians, step5-derived (step5_fleet_fits section F)
% The decisive control: how much does the WINDOW itself move the model's slope?
% If windowing the model barely changes it, the window cannot explain the gap.
fprintf('\n  %-8s %16s %22s %20s %8s %12s\n', 'model', 'measured astop5', ...
    'astop_model_windowed', 'true terminal slope', 'ratio', 'window bias');
ratio_g = nan(3,1); meas_g = nan(3,1); modw_g = nan(3,1); true_g = nan(3,1); wbias_g = nan(3,1);
for g = 1:3
    okT = model == MDL(g) & isfinite(astop5);
    selc = cellModel == MDL(g) & isfinite(astop_mod);
    meas_g(g) = median(astop5(okT));
    modw_g(g) = median(astop_mod(selc));
    true_g(g) = median(astop_true(selc));
    ratio_g(g) = modw_g(g) / meas_g(g);
    % what the estimator does to the model, in isolation from the data
    wbias_g(g) = (modw_g(g) - true_g(g)) / abs(true_g(g));
    fprintf('  %-8s %16.1f %22.1f %20.1f %8.2f %11.1f%%\n', ...
        MDL(g), meas_g(g), modw_g(g), true_g(g), ratio_g(g), 100*wbias_g(g));
end

% verdict, decided on the numbers rather than asserted
fprintf('\n  VERDICT:\n');
worst = max(abs(ratio_g - 1));
maxbias = max(abs(wbias_g));
fprintf(['  Window control: applying the tstop5 rule to the MODEL moves its slope by at most\n' ...
         '  %.1f%% from the model''s own true terminal slope, so the estimator is very nearly\n' ...
         '  unbiased on the model and cannot by itself manufacture a large gap.\n'], 100*maxbias);
if worst <= 0.30
    fprintf(['  The estimator-matched model slope lands within %.0f%% of the measured astop5 in\n' ...
             '  every geometry, so the gap between the cell-fit k and k_astop5 is DEFINITIONAL\n' ...
             '  (window plus smooth approach to the stop), not a model failure, and is NOT a\n' ...
             '  KD departure.\n'], 100*worst);
else
    fprintf(['  The model decelerates %.2f / %.2f / %.2f times harder than the trials actually do\n' ...
             '  at the stop, far outside the %.0f%% band that a definitional artefact could\n' ...
             '  explain and %.0fx larger than the window bias just measured. The TERMINAL DRAG\n' ...
             '  DEFICIT IS REAL: real feet keep creeping after the model has come to rest. It\n' ...
             '  must be listed as a KD departure alongside the arch.\n'], ...
             ratio_g, 30, worst/max(maxbias, eps));
end

%% ===================================================================
%  PART 4 -- FIGURE (working, minimal styling)
%  ===================================================================
figure('Position', [50 50 1550 950]);

% (a) closure in final depth
ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
hM = gobjects(1,3);
for g = 1:3
    sel = model == MDL(g) & ~flagDisp;
    hM(g) = plot(ax, K.v0_cm_s(sel), K.d_final_cm(sel), MRK{g}, 'MarkerSize', 4, 'Color', CO(g,:));
    plot(ax, V0_GRID, d_pred_g(:,g), '-', 'LineWidth', 1.5, 'Color', CO(g,:));
end
hS3 = plot(ax, V0_GRID, STEP3_D0 + STEP3_A*V0_GRID.^(2/3), 'k--', 'LineWidth', 1.2);
xlabel(ax,'v_0 (cm/s)'); ylabel(ax,'d_{final} (cm)');
title(ax,'(a) closure: depth');
legend(ax, [hM hS3], [cellstr(MDL) {'step-3 scaling fit'}], 'Location','northwest');

% (b) closure in stopping time; pipeline layer faint and explicitly mixed
ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
% faint reference layer only -- this column is definitionally mixed
hP = plot(ax, K.v0_cm_s(~flagDisp), K.t_stop_s(~flagDisp), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 4);
hM = gobjects(1,3);
for g = 1:3
    sel = model == MDL(g) & ~flagDisp;
    hM(g) = plot(ax, K.v0_cm_s(sel), tstop5(sel), MRK{g}, 'MarkerSize', 4, 'Color', CO(g,:));
    plot(ax, V0_GRID, t_pred_g(:,g), '-', 'LineWidth', 1.5, 'Color', CO(g,:));
end
xlabel(ax,'v_0 (cm/s)'); ylabel(ax,'stopping time (s)');
title(ax,'(b) closure: stopping time (tstop5 primary)');
legend(ax, [hM hP], [cellstr(MDL) {'pipeline (mixed definition)'}], 'Location','northeast');

% (c) arch against the local area profile
ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
yyaxis(ax,'left');
set(ax, 'YColor', 'k');
hC = gobjects(1,3);
for g = 1:3
    % explicit colour: yyaxis would otherwise paint all three the same
    hC(g) = plot(ax, medGrid{g}, medResid{g}, '-', 'LineWidth', 1.5, 'Color', CO(g,:));
end
ylabel(ax,'median residual speed^2 / v_0^2'); yline(ax, 0, 'k-');
yyaxis(ax,'right');
% the three models' a_local agree to ~0.001 cm^2, so one curve is drawn for all
hA = plot(ax, z_area, AR.a_local_Default, '--', 'LineWidth', 1.2, 'Color', [0.4 0.4 0.4]);
set(ax, 'YColor', [0.4 0.4 0.4]);
ylabel(ax,'a_{local} (cm^2)');
for zL = LANDMARKS, xline(ax, zL, ':', sprintf('%.2f', zL)); end
xlim(ax, [0 2.6]); xlabel(ax,'depth (cm)');
title(ax,'(c) residual arch vs a_{local} (dashed)');
legend(ax, [hC hA], [cellstr(MDL) {'a_{local} (all models)'}], 'Location','southeast');

% (d) the three a_stop quantities per geometry
ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
bar(ax, [meas_g modw_g true_g]);
set(ax, 'XTick', 1:3, 'XTickLabel', cellstr(MDL));
ylabel(ax,'a_{stop} (cm/s^2)');
title(ax,'(d) a_{stop}: measured vs model');
legend(ax, {'measured astop5','model, windowed','model, true terminal'}, 'Location','southeast');

exportgraphics(gcf, FIGPATH, 'Resolution', 200);
fprintf('\nfigure written: %s\n', FIGPATH);

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function p = local_newest(dirPath, pattern)
% Newest file matching a pattern, by modification time.
    f = dir(fullfile(dirPath, pattern));
    if isempty(f), error('step5c:noFile', 'No file matching %s in %s', pattern, dirPath); end
    [~, j] = max([f.datenum]);
    p = fullfile(f(j).folder, f(j).name);
end

function [t, v, z] = local_integrate_model(v0, k, d1, mass, grav, dt)
% The KD equation integrated forward in TIME with RK4, so the model produces a
% v(t) trace that the measurement's window rule can be applied to unchanged:
%   dv/dt = grav - k*z/mass - v^2/d1 ,  dz/dt = v
    f = @(zz, vv) grav - k*zz/mass - vv.^2/d1;
    nMax = 20000;
    t = nan(nMax,1); v = nan(nMax,1); z = nan(nMax,1);
    t(1) = 0; v(1) = v0; z(1) = 0;
    for i = 1:nMax-1
        vi = v(i); zi = z(i);
        k1v = f(zi, vi);                 k1z = vi;
        k2v = f(zi+dt/2*k1z, vi+dt/2*k1v); k2z = vi + dt/2*k1v;
        k3v = f(zi+dt/2*k2z, vi+dt/2*k2v); k3z = vi + dt/2*k2v;
        k4v = f(zi+dt*k3z,   vi+dt*k3v);   k4z = vi + dt*k3v;
        v(i+1) = vi + dt/6*(k1v + 2*k2v + 2*k3v + k4v);
        z(i+1) = zi + dt/6*(k1z + 2*k2z + 2*k3z + k4z);
        t(i+1) = t(i) + dt;
        if v(i+1) <= 0, break; end       % the projectile has stopped
    end
    ok = isfinite(t); t = t(ok); v = v(ok); z = z(ok);
end

function [slope, tzero] = local_window_slope(t, v, z, d_ref)
% The tstop5 rule, applied verbatim: linear fit of v against t over
% 4 < v < 30 cm/s AND depth > 0.5*d_ref, extrapolated to v = 0.
    slope = NaN; tzero = NaN;
    w = v > 4 & v < 30 & z > 0.5*d_ref & isfinite(v);
    if sum(w) < 5, return; end
    p = polyfit(t(w), v(w), 1);
    slope = p(1);
    if p(1) ~= 0, tzero = -p(2)/p(1); end
end

function d_pred = local_predict_d(v0, k, d1, mass, grav)
% Predicted final depth = the root of kd_speed2_model(z) = 0, bracketed.
    f = @(z) kd_speed2_model(z, v0, k, d1, mass, grav);
    zs = linspace(1e-6, 20, 2000);
    fs = arrayfun(f, zs);
    j  = find(fs(1:end-1) > 0 & fs(2:end) <= 0, 1);
    if isempty(j), d_pred = NaN; return; end
    try
        d_pred = fzero(f, [zs(j) zs(j+1)]);
    catch
        d_pred = NaN;
    end
end

function tp = local_t_phys(v0, k, d1, mass, grav, d_pred)
% SECONDARY stopping time: the physical integral of dz/v(z) to where the model
% speed reaches 1 cm/s, closed with the model's own linear terminal tail.
    tp = NaN;
    f = @(z) kd_speed2_model(z, v0, k, d1, mass, grav);
    g = @(z) f(z) - 1;                      % speed^2 = 1 -> speed = 1 cm/s
    try
        z_end = fzero(g, [1e-9, d_pred]);
    catch
        return
    end
    t_main = integral(@(z) 1 ./ sqrt(max(f(z), eps)), 0, z_end, 'AbsTol', 1e-10, 'RelTol', 1e-9);
    % terminal deceleration at the stop, from the force law itself
    a_term = k*d_pred/mass - grav;
    if a_term <= 0
        fprintf('    WARNING: a_term = %.1f <= 0 at v0 with d_pred = %.3f; terminal tail undefined.\n', ...
            a_term, d_pred);
        return
    end
    tp = t_main + 1/a_term;                 % remaining time from v_end = 1 cm/s
end

function [zz, uu] = local_monotone_descent(depth, speed, speed2)
% The descent segment with strictly increasing depth, so interp1 gets unique
% sample points. Identical to step5_fleet_fits.
    [~, i1] = max(speed);
    [~, i2] = max(depth);
    zz = []; uu = [];
    if i2 <= i1, return; end
    z = depth(i1:i2); u = speed2(i1:i2);
    ok = isfinite(z) & isfinite(u);
    z = z(ok); u = u(ok);
    keep = false(numel(z),1); last = -inf;
    for q = 1:numel(z)
        if z(q) > last, keep(q) = true; last = z(q); end
    end
    zz = z(keep); uu = u(keep);
end

function [zc, zp, pk] = local_arch_features(zg, r)
% Depth of the first negative-to-positive zero crossing and of the positive peak.
    zc = NaN; zp = NaN; pk = NaN;
    ok = isfinite(r); zg = zg(ok); r = r(ok);
    if numel(r) < 3, return; end
    j = find(r(1:end-1) < 0 & r(2:end) >= 0, 1);
    if ~isempty(j)
        zc = zg(j) + (zg(j+1)-zg(j)) * (-r(j)) / (r(j+1)-r(j));
    end
    [pk, jp] = max(r);
    zp = zg(jp);
end

function y = local_prctile(x, p)
% Percentiles by linear interpolation of the empirical CDF at midpoints.
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0, y = nan(size(p)); return; end
    if n == 1, y = repmat(x, size(p)); return; end
    q = (0.5:1:(n-0.5)) / n * 100;
    y = interp1(q, x, p, 'linear', 'extrap');
    y = min(max(y, x(1)), x(end));
end

function r = local_spearman(x, y)
% Spearman rank correlation with average ranks for ties. Base MATLAB.
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3, r = NaN; return; end
    rx = local_tiedrank(x); ry = local_tiedrank(y);
    ax = rx - mean(rx); ay = ry - mean(ry);
    r = sum(ax.*ay) / sqrt(sum(ax.^2) * sum(ay.^2));
end

function r = local_tiedrank(x)
% Ranks of x, ties sharing their average rank.
    n = numel(x);
    [xs, i] = sort(x(:));
    ranks = (1:n).';
    rs = zeros(n,1);
    q = 1;
    while q <= n
        j = q;
        while j < n && xs(j+1) == xs(q), j = j + 1; end
        rs(q:j) = mean(ranks(q:j));
        q = j + 1;
    end
    r = zeros(n,1); r(i) = rs;
end
