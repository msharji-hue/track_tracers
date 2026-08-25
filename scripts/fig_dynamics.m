% FIG_DYNAMICS -- the four final dynamic-study figures.
%
% Presentation plus the prescribed decomposition arithmetic only. NOTHING is
% refitted: (k, c) come from step5_fleet_fits and the per-cell CSV, the force
% law is untouched, and no data, exclusion or pipeline output is modified. The
% only files written are the eight figure files named below.
%
%   DYN-A  fig_tstop                 stopping time vs v0, per geometry
%   DYN-B  fig_scaling_models        depth scaling with three model curves
%   DYN-C  fig_force_decomposition   UFL Fig-3a/3b analogues on corrected a+g
%   DYN-D  fig_parameter_ladder      cell parameters and the model ladder
%
% AUTHORITATIVE INPUT. master_trials_20260824_221720 -- the export in which the
% a_plus_g -2g convention bug is FIXED (a+g = g - a, rest state +g). Older
% stamps must not be read; see PR "Fix -2g a_plus_g convention in master export".
%
% PUBLIC NOTATION. The drag parameter is written c = m/d1 everywhere a reader
% can see it. d1 appears only inside the solver, never in figure text.
%
% Base MATLAB only.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));  % shared KD helpers

% -- constants --------------------------------------------------------
mass = 65;      % g
grav = 980;     % cm/s^2

MASTER = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260824_221720.mat';
EXPDIR = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports';
OUTDIR = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures';

% shared M0, trial-balanced (step5_fleet_fits section C). NOT refitted here.
k0   = 1.9163e5;      % dyn/cm
c0   = 15.898;        % g/cm   -- the public drag parameter
d1_0 = mass / c0;     % cm     -- solver-internal only, never printed on a figure

% step-3 published coefficients, reproduced below and asserted to 3 s.f.
STEP3_D0 = 0.044; STEP3_A = 0.0711; STEP3_N = 0.606;

% Geometry identity -- copied from scripts/fig_scaling.m lines 101-106 so the
% two figures cannot drift. Verified: MARK = o/s/^ and COL rows exactly as
% fig_scaling defines them for Tight/Default/Wide.
MDL  = ["Tight" "Default" "Wide"];
MARK = ["o" "s" "^"];
COL  = [0.8500 0.3250 0.0980;      % Tight
        0      0.4470 0.7410;      % Default
        0.9290 0.6940 0.1250];     % Wide
MSZ  = 4;
KDCOL = [0.4940 0.1840 0.5560];    % KD curve: distinct from every geometry

V_MIN = 15;                        % cm/s, the frame mask used throughout DYN-C
ZC    = [0.3 0.5 0.7 0.9];         % DYN-C row-1 bin centres
ZHALF = 0.1;

fprintf('\n=== fig_dynamics ===\n');
fprintf('  input : %s\n', MASTER);
fprintf('  outdir: %s\n', OUTDIR);

%% ===================================================================
%  LOAD + SANITY ASSERTIONS
%  ===================================================================
L = load(MASTER);
T = L.T;  S = L.S;
keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
model = string(K.model);
tagsS = string({S.trialTag});

% display set: markers only. Fits and printed statistics keep all 524.
flagDisp = K.d_final_cm < 1 | K.v0_cm_s < 55;
fprintf('\n  analysis set %d trials | display set %d (omits %d flagged, markers only)\n', ...
    nTr, sum(~flagDisp), sum(flagDisp));

% ---- a_plus_g convention assertion -------------------------------------
% The prescribed check is the median a+g over the last 20 finite frames against
% +980. That band cannot be met by construction: kd_kinematics NaNs a+g outside
% [impact_index, stopFrame], so the retained record STOPS at the stop and never
% contains the rest state. The last finite frames are still decelerating, where
% a+g = g + |a| > g. The quantity is reported, and the assertion that actually
% discriminates the convention is enforced instead: a+g must equal g - a on
% every frame, and must sit nearer +g than -g (the old formula gave -g at rest).
j1 = find(tagsS == string(K.trialTag(1)), 1);
ap1 = S(j1).a_plus_g_cm_s2(:);  f1 = ap1(isfinite(ap1));
m20 = median(f1(max(1,end-19):end));
fprintf('\n  a+g sanity, trial %s:\n', K.trialTag(1));
fprintf('    median over last 20 finite frames = %.1f cm/s^2 (prescribed target +980 +/- 100)\n', m20);
if abs(m20 - grav) >= 100
    fprintf(['    NOTE: outside the +/-100 band, and that is EXPECTED, not a failure --\n' ...
             '    a+g is NaN past stopFrame so the record never reaches rest. Enforcing the\n' ...
             '    discriminating checks instead (identity with g - a, and sign).\n']);
end
nExact = 0;
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    ap = S(j).a_plus_g_cm_s2(:); a = S(j).a_cm_s2(:);
    o = isfinite(ap) & isfinite(a);
    if any(o) && max(abs(ap(o) - (grav - a(o)))) < 1e-9, nExact = nExact + 1; end
end
fprintf('    a+g == g - a exactly on every finite frame: %d / %d trials\n', nExact, nTr);
if nExact ~= nTr
    error('figdyn:apgConvention', ...
        ['a+g does not equal g - a on %d of %d trials. This export carries the stale ' ...
         '(-a - g) convention; read the corrected stamp.'], nTr-nExact, nTr);
end
if ~(m20 > 0 && abs(m20-grav) < abs(m20+grav))
    error('figdyn:apgSign', ...
        ['a+g is nearer -g than +g (median %.1f): the stale convention is present.'], m20);
end
fprintf('    PASS: convention is g - a with rest state +g.\n');

% ---- step-3 coefficient reproduction ------------------------------------
x23 = K.v0_cm_s.^(2/3);  y = K.d_final_cm;  N = numel(y);
c23 = [ones(N,1) x23] \ y;
rssN = @(n) local_rss_pow(n, K.v0_cm_s, y);
nfree = fminbnd(rssN, 0.3, 1.2);
cfree = [ones(N,1) K.v0_cm_s.^nfree] \ y;
fprintf('\n  step-3 reproduction (all %d trials):\n', N);
fprintf('    fixed 2/3 : d0 = %.4f (published %.3f) | a = %.5f (published %.4f)\n', ...
    c23(1), STEP3_D0, c23(2), STEP3_A);
fprintf('    free n    : n  = %.4f (published %.3f)\n', nfree, STEP3_N);
local_assert3sf('step-3 d0', c23(1), STEP3_D0);
local_assert3sf('step-3 a',  c23(2), STEP3_A);
local_assert3sf('step-3 n',  nfree,  STEP3_N);
fprintf('    PASS: all three reproduce to 3 significant figures.\n');

% ---- per-cell fits (already fitted; read, never refitted) ---------------
cellCsv = local_newest(EXPDIR, 'step5_cell_fits_*.csv');
C = readtable(cellCsv);  C.model = string(C.model);
fprintf('\n  cell fits: %s (%d cells)\n', cellCsv, height(C));

%% ===================================================================
%  DYN-A  fig_tstop
%  ===================================================================
% CANONICAL stopping time = the pipeline t_stop_s. One definition for the whole
% manuscript, identical to fig_kinematics and to Katsuragi & Durian (2007):
% linear extrapolation of v(t) to zero over the pre-stop window above twice the
% maximum rebound speed. The Step-5 window estimators are sensitivity
% diagnostics and appear nowhere in this figure.
fprintf('\n=== DYN-A: fig_tstop ===\n');

% -- DEFINITION QA -------------------------------------------------------
% src/kd_kinematics.m find_stop: the rebound-gated segment is built backwards
% from the first zero crossing while v > rebFactor*rebound. The slope is fitted
% only when that segment has >= 2 frames, and fit_slope itself returns NaN below
% 3 points; when no slope is produced the function FALLS BACK to t_stop =
% t(cross), the raw first-zero-crossing frame time. a_stop is assigned in that
% same branch and nowhere else, so a_stop = NaN is a RECORDED marker that the
% fallback fired -- no reprocessing is needed to count it.
fallbackFrac = nan(1,3);
for g = 1:3
    s = model == MDL(g);
    fallbackFrac(g) = mean(~isfinite(K.a_stop_cm_s2(s)));
end
fprintf('  fallback rule: extrapolation segment too short -> t_stop = t(first zero crossing)\n');
fprintf('  branch recorded? kin.method is a single constant string for all %d trials and does\n', nTr);
fprintf('  NOT record it; a_stop = NaN does, since a_stop is set only in the extrapolation branch.\n');
for g = 1:3
    fprintf('    %-8s fallback fired on %5.1f%% of trials (%d of %d)\n', MDL(g), ...
        100*fallbackFrac(g), sum(~isfinite(K.a_stop_cm_s2(model==MDL(g)))), sum(model==MDL(g)));
end
fprintf('    overall  %5.1f%%\n', 100*mean(~isfinite(K.a_stop_cm_s2)));

figA = figure('Color','w','Units','inches','Position',[1 1 6.75 2.85], 'PaperPositionMode','auto');
axA = gobjects(1,3);
for g = 1:3
    axA(g) = subplot(1,3,g); hold on; grid on; box on;
    sd = model == MDL(g) & ~flagDisp;
    % per-trial cloud, faint
    scatter(K.v0_cm_s(sd), K.t_stop_s(sd), 14, 'Marker', char(MARK(g)), ...
        'MarkerEdgeColor', COL(g,:), 'MarkerEdgeAlpha', 0.30);
    % per-height medians, full weight
    sa = model == MDL(g);
    hh = unique(K.dropHeight_mm(sa));
    mv = nan(numel(hh),1); mt = nan(numel(hh),1);
    for q = 1:numel(hh)
        s2 = sa & K.dropHeight_mm == hh(q);
        mv(q) = median(K.v0_cm_s(s2), 'omitnan');
        mt(q) = median(K.t_stop_s(s2), 'omitnan');
    end
    plot(mv, mt, char(MARK(g)), 'MarkerSize', 6, 'MarkerFaceColor', COL(g,:), ...
        'MarkerEdgeColor', COL(g,:), 'LineStyle','none');
    xlabel('v_0 (cm/s)');
    if g == 1, ylabel('stopping time (s)'); end
    title(MDL(g), 'FontWeight','normal');
end
linkaxes(axA, 'y');
ylim(axA(1), [0 max(K.t_stop_s(~flagDisp))*1.05]);
pA = local_save(figA, OUTDIR, 'fig_tstop');

% CAPTION-NOTES (DYN-A)
%   Stopping time is the pipeline t_stop_s throughout: KD 2007's linear
%   extrapolation of v(t) to zero over the pre-stop window above twice the
%   maximum rebound speed (UFL Supplemental Fig. 1), the same definition
%   fig_kinematics uses.
%   Fallback: when the rebound-gated segment cannot form (fewer than the frames
%   the slope fit needs), find_stop falls back to t_stop = t(first zero
%   crossing). The branch is not recorded directly by kin.method, but a_stop is
%   assigned only in the extrapolation branch, so a_stop = NaN counts it: the
%   fallback fired on 98.9 / 98.3 / 98.9 % of Tight / Default / Wide trials.
%   The definition hierarchy should be stated in Methods.
%   Sensitivity: window-based re-estimates (Step-5 diagnostics) shift low-v0
%   values by up to ~30%, so the MAGNITUDE of the declining trend is
%   estimator-sensitive; its direction and the high-v0 flattening are not.

%% ===================================================================
%  DYN-B  fig_scaling_models
%  ===================================================================
fprintf('\n=== DYN-B: fig_scaling_models ===\n');
vg = linspace(min(K.v0_cm_s), max(K.v0_cm_s), 300).';
dKD = arrayfun(@(v) local_predict_d(v, k0, d1_0, mass, grav), vg);

figB = figure('Color','w','Units','inches','Position',[1 1 3.375 2.85], 'PaperPositionMode','auto');
hold on; grid on; box on;
hG = gobjects(1,3);
for g = 1:3
    sd = model == MDL(g) & ~flagDisp;
    hG(g) = plot(K.v0_cm_s(sd), K.d_final_cm(sd), char(MARK(g)), 'MarkerSize', MSZ, ...
        'MarkerFaceColor', COL(g,:), 'MarkerEdgeColor', COL(g,:), 'LineStyle','none');
end
h1 = plot(vg, c23(1) + c23(2)*vg.^(2/3), 'k-',  'LineWidth', 1.4);
h2 = plot(vg, cfree(1) + cfree(2)*vg.^nfree, 'k--', 'LineWidth', 1.2);
h3 = plot(vg, dKD, '-', 'Color', KDCOL, 'LineWidth', 1.6);
xlabel('v_0 (cm/s)'); ylabel('d (cm)');
legend([hG h1 h2 h3], [cellstr(MDL), ...
    {sprintf('fixed 2/3 fit'), sprintf('free exponent n = %.3f', nfree), ...
     'KD force law (no free parameters)'}], ...
    'Location','northwest', 'FontSize', 7, 'Box','off');
% headroom so the legend text clears the curves it sits over
ylim([0.5 4.5]);
pB = local_save(figB, OUTDIR, 'fig_scaling_models');

%% ===================================================================
%  DYN-C  fig_force_decomposition
%  ===================================================================
fprintf('\n=== DYN-C: fig_force_decomposition ===\n');

figC = figure('Color','w','Units','inches','Position',[1 1 6.75 5.2], 'PaperPositionMode','auto');

% ---- Row 1: the independent check on c ---------------------------------
cBins = nan(1,3); pPar = nan(1,3);
for g = 1:3
    ax = subplot(2,3,g); hold on; grid on; box on;
    % one point per trial per bin
    v2p = []; agp = []; binp = [];
    idx = find(model == MDL(g));
    for q = 1:numel(idx)
        i = idx(q);
        j = find(tagsS == string(K.trialTag(i)), 1);
        z = S(j).z_cm(:); v = S(j).v_cm_s(:); ag = S(j).a_plus_g_cm_s2(:);
        for b = 1:4
            m = abs(z - ZC(b)) <= ZHALF & v > V_MIN & isfinite(ag);
            if ~any(m), continue; end
            v2p(end+1,1) = mean(v(m).^2);   %#ok<AGROW>
            agp(end+1,1) = mean(ag(m));     %#ok<AGROW>
            binp(end+1,1) = b;              %#ok<AGROW>
        end
    end
    % common-slope joint fit: bin indicators + one shared slope
    Ind = double(binp == 1:4);
    Xc  = [Ind v2p];
    bc  = Xc \ agp;
    rssC = sum((agp - Xc*bc).^2);
    cBins(g) = mass * bc(5);
    % separate slopes, for the parallel-lines F test
    Xf  = [Ind Ind.*v2p];
    bf  = Xf \ agp;
    rssF = sum((agp - Xf*bf).^2);
    nP = numel(agp);
    Fst = ((rssC - rssF)/(8-5)) / (rssF/(nP-8));
    pPar(g) = local_fpval(Fst, 3, nP-8);
    % four shades of the geometry colour, palest shallowest
    wsh = [0.55 0.37 0.18 0];
    hB = gobjects(1,4);
    for b = 1:4
        sh = COL(g,:) + (1-COL(g,:))*wsh(b);
        m = binp == b;
        hB(b) = plot(v2p(m), agp(m), char(MARK(g)), 'MarkerSize', 3, ...
            'MarkerFaceColor', sh, 'MarkerEdgeColor', sh, 'LineStyle','none');
        xr = linspace(min(v2p(m)), max(v2p(m)), 20);
        plot(xr, bc(b) + bc(5)*xr, '-', 'Color', sh, 'LineWidth', 1.0);
    end
    xlabel('v^2 (cm^2/s^2)');
    if g == 1, ylabel('a+g (cm/s^2)'); end
    title(MDL(g), 'FontWeight','normal');
    text(0.04, 0.95, sprintf('c_{bins} = %.1f g/cm\nparallel p = %.3f', cBins(g), pPar(g)), ...
        'Units','normalized', 'VerticalAlignment','top', 'FontSize', 7);
    if g == 3
        legend(hB, arrayfun(@(b) sprintf('z = %.1f cm', ZC(b)), 1:4, 'UniformOutput', false), ...
            'Location','southeast', 'FontSize', 6, 'Box','off');
    end
end

% ---- Row 2: consistency display ----------------------------------------
ax2 = subplot(2,1,2); hold on; grid on; box on;
zEdges = 0.1:0.1:3.2;
hC = gobjects(1,3);
localSlope = nan(1,3);
for g = 1:3
    zz = []; yy = [];
    idx = find(model == MDL(g));
    for q = 1:numel(idx)
        j = find(tagsS == string(K.trialTag(idx(q))), 1);
        z = S(j).z_cm(:); v = S(j).v_cm_s(:); ag = S(j).a_plus_g_cm_s2(:);
        m = v > V_MIN & isfinite(ag);
        zz = [zz; z(m)]; yy = [yy; ag(m) - v(m).^2*c0/mass];   %#ok<AGROW>
    end
    zm = nan(numel(zEdges)-1,1); ym = nan(numel(zEdges)-1,1);
    for b = 1:numel(zEdges)-1
        m = zz > zEdges(b) & zz <= zEdges(b+1);
        if sum(m) >= 200                       % retain only well-populated bins
            zm(b) = 0.5*(zEdges(b)+zEdges(b+1));
            ym(b) = median(yy(m));
        end
    end
    hC(g) = plot(zm, ym, '-', 'Color', COL(g,:), 'LineWidth', 1.4);
    % local slope over 1.3-2.2 cm, the region between the two marks
    w = isfinite(zm) & zm >= 1.3 & zm <= 2.2;
    pl = [ones(sum(w),1) zm(w)] \ ym(w);
    localSlope(g) = pl(2);
end
zl = linspace(0.1, 3.2, 50);
hK = plot(zl, k0*zl/mass, 'k-', 'LineWidth', 1.2);
xline(1.0, ':', 'Color', [0.6 0.6 0.6]);       % light marks, reference only
xline(2.4, ':', 'Color', [0.6 0.6 0.6]);
xlabel('z (cm)'); ylabel('(a+g) - c_0 v^2/m  (cm/s^2)');
legend([hC hK], [cellstr(MDL), {'trajectory-fit k (not refit)'}], ...
    'Location','northwest', 'FontSize', 7, 'Box','off');
pC = local_save(figC, OUTDIR, 'fig_force_decomposition');

% ---- for the record -----------------------------------------------------
fprintf('\n  c_bins vs c0 = %.3f g/cm:\n', c0);
for g = 1:3
    fprintf('    %-8s c_bins = %6.2f g/cm  (%+.1f%% vs c0) | parallel-lines p = %.4f\n', ...
        MDL(g), cBins(g), 100*(cBins(g)-c0)/c0, pPar(g));
end
fprintf('  row-2 local slope over 1.3-2.2 cm vs k0 = %.4e:\n', k0);
for g = 1:3
    fprintf('    %-8s slope = %8.1f cm/s^2 per cm -> m*slope = %.4e (%+.1f%% vs k0)\n', ...
        MDL(g), localSlope(g), mass*localSlope(g), 100*(mass*localSlope(g)-k0)/k0);
end

% CAPTION-NOTES (DYN-C)
%   Row 1 is the independent check on c, NOT a measurement of k: per-bin
%   intercepts are noise-limited and are not presented as independent k
%   measurements (within-trial corr(z, v^2) = -0.994; design condition number
%   12566). Only the common slope, reported as c_bins = m*s, is read off it.
%   Row 2 is a consistency display against the trajectory-fit k, not a refit.
%   The shallow region is estimator-limited -- the two-stage differentiation
%   biases it, and a synthetic null built from the shared parameters reproduces
%   ~88% of the intercept excess -- and the deep rolloff reflects near-stop
%   dominance. a+g is the corrected g - a throughout.

%% ===================================================================
%  DYN-D  fig_parameter_ladder
%  ===================================================================
fprintf('\n=== DYN-D: fig_parameter_ladder ===\n');
okC = isfinite(C.k_fit) & isfinite(C.d1_fit);
cCell = mass ./ C.d1_fit;          % per cell first, medianed after

figD = figure('Color','w','Units','inches','Position',[1 1 6.75 2.85], 'PaperPositionMode','auto');

% (a) cell k by geometry
subplot(1,3,1); hold on; grid on; box on;
for g = 1:3
    s = okC & C.model == MDL(g);
    xj = g + 0.16*(rand(sum(s),1)-0.5);
    % plotted in units of 1e5 so no axis exponent label can collide with the title
    plot(xj, C.k_fit(s)/1e5, char(MARK(g)), 'MarkerSize', 3.5, ...
        'MarkerFaceColor', COL(g,:), 'MarkerEdgeColor', COL(g,:), 'LineStyle','none');
    plot(g+[-0.28 0.28], median(C.k_fit(s))/1e5*[1 1], 'k-', 'LineWidth', 1.6);
end
set(gca,'XTick',1:3,'XTickLabel',cellstr(MDL)); xlim([0.5 3.5]);
ylabel('cell k (10^5 dyn/cm)'); title('(a) k by geometry','FontWeight','normal');

% (b) cell c by geometry -- c = m/d1 per cell, then medianed
subplot(1,3,2); hold on; grid on; box on;
for g = 1:3
    s = okC & C.model == MDL(g);
    xj = g + 0.16*(rand(sum(s),1)-0.5);
    plot(xj, cCell(s), char(MARK(g)), 'MarkerSize', 3.5, ...
        'MarkerFaceColor', COL(g,:), 'MarkerEdgeColor', COL(g,:), 'LineStyle','none');
    plot(g+[-0.28 0.28], median(cCell(s))*[1 1], 'k-', 'LineWidth', 1.6);
end
set(gca,'XTick',1:3,'XTickLabel',cellstr(MDL)); xlim([0.5 3.5]);
ylabel('cell c = m/d_1 (g/cm)'); title('(b) c by geometry','FontWeight','normal');

% (c) the ladder: dBIC relative to M0, with the selection fractions
LAD   = ["M0" "M1" "M2" "M3"];
LADP  = [2 4 4 6];
wRMSE = [5000.9 4995.4 4998.2 4978.7];
BICv  = [8938.7 8950.1 8950.7 8959.1];
SELp  = [97.4 0.2 2.0 0.4];
dBIC  = BICv - BICv(1);
subplot(1,3,3); hold on; grid on; box on;
bh = bar(1:4, dBIC, 0.6, 'EdgeColor', 'k');
% every row of CData must be set, or the unassigned bars fall back to the
% default colour and the M0 highlight silently does nothing
bh.FaceColor = 'flat';
bh.CData = repmat([0.75 0.75 0.75], 4, 1);
bh.CData(1,:) = KDCOL;                       % M0, the selected model
% M0's bar has zero height by construction (dBIC against itself), so the
% selected model is marked on its label rather than by bar colour.
for q = 1:4
    if q == 1
        text(q, dBIC(q) + 1.0, sprintf('%.1f%%', SELp(q)), 'HorizontalAlignment','center', ...
            'FontSize', 8, 'FontWeight','bold', 'Color', KDCOL);
    else
        text(q, dBIC(q) + 1.0, sprintf('%.1f%%', SELp(q)), 'HorizontalAlignment','center', ...
            'FontSize', 7);
    end
end
% single-line labels: an embedded newline is split across consecutive ticks
set(gca,'XTick',1:4,'XTickLabel', arrayfun(@(q) sprintf('%s (%d)', LAD(q), LADP(q)), ...
    1:4, 'UniformOutput', false), 'XTickLabelRotation', 20);
ylabel('\DeltaBIC vs M0'); ylim([-2 max(dBIC)+5]);
title('(c) model ladder','FontWeight','normal');
text(0.5, 0.92, 'M0 selected in 97.4% of resamples', 'Units','normalized', ...
    'HorizontalAlignment','center', 'FontSize', 8, 'FontWeight','bold');
pD = local_save(figD, OUTDIR, 'fig_parameter_ladder');

fprintf('\n  ladder (given, not refitted): wRMSE %s | BIC %s | selection %s\n', ...
    mat2str(wRMSE), mat2str(BICv), mat2str(SELp));
fprintf('  cell medians: k = %s | c = %s\n', ...
    mat2str(arrayfun(@(g) median(C.k_fit(okC & C.model==MDL(g))), 1:3, 'UniformOutput', true), 5), ...
    mat2str(arrayfun(@(g) median(cCell(okC & C.model==MDL(g))), 1:3, 'UniformOutput', true), 4));

%% ===================================================================
%  OUTPUT PATHS
%  ===================================================================
fprintf('\n=== WRITTEN ===\n');
for p = [pA pB pC pD], fprintf('  %s\n', p); end

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function paths = local_save(fig, outDir, stem)
% Vector PDF + 300 dpi PNG, matching fig_scaling. Refuses to overwrite.
    pdfPath = fullfile(outDir, [stem '.pdf']);
    pngPath = fullfile(outDir, [stem '.png']);
    if isfile(pdfPath) || isfile(pngPath)
        error('figdyn:wouldOverwrite', ...
            'Refusing to overwrite an existing figure: %s(.pdf/.png)', fullfile(outDir, stem));
    end
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 300);
    paths = [string(pdfPath) string(pngPath)];
end

function r = local_rss_pow(n, v0, y)
% RSS of d = d0 + a*v0^n at fixed n, coefficients concentrated out.
    X = [ones(numel(y),1) v0.^n];
    r = sum((y - X*(X\y)).^2);
end

function local_assert3sf(name, got, want)
% Assert agreement with the published coefficient at 3-significant-figure
% precision, judged as a RELATIVE difference below 5e-3. An absolute
% 3-s.f. test is not usable here because the published values are quoted at
% mixed precision (d0 = 0.044 carries only two significant figures), so the
% reproduced 0.04414 would fail a literal digit comparison while rounding to
% exactly the published 0.044.
    rel = abs(got - want) / abs(want);
    fprintf('      %-12s reproduced %.6g vs published %.6g  (relative %.2e)\n', name, got, want, rel);
    if rel > 5e-3
        error('figdyn:assert3sf', ...
            '%s: reproduced %.6g but published %.6g (relative %.2e exceeds 5e-3)', ...
            name, got, want, rel);
    end
end

function d_pred = local_predict_d(v0, k, d1, mass, grav)
% Predicted final depth = the root of kd_speed2_model(z) = 0.
    f = @(z) kd_speed2_model(z, v0, k, d1, mass, grav);
    zs = linspace(1e-6, 20, 2000);
    fs = arrayfun(f, zs);
    j  = find(fs(1:end-1) > 0 & fs(2:end) <= 0, 1);
    if isempty(j), d_pred = NaN; return; end
    try, d_pred = fzero(f, [zs(j) zs(j+1)]); catch, d_pred = NaN; end
end

function p = local_fpval(F, df1, df2)
% Upper-tail F probability via the regularised incomplete beta (base MATLAB):
%   P(X > F) = I_{df2/(df2+df1 F)}(df2/2, df1/2)
    if ~isfinite(F) || F <= 0, p = 1; return; end
    p = betainc(df2/(df2 + df1*F), df2/2, df1/2);
end

function p = local_newest(dirPath, pattern)
% Newest file matching a pattern, by modification time.
    f = dir(fullfile(dirPath, pattern));
    if isempty(f), error('figdyn:noFile', 'No file matching %s in %s', pattern, dirPath); end
    [~, j] = max([f.datenum]);
    p = fullfile(f(j).folder, f(j).name);
end
