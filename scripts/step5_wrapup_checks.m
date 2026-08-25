% STEP 5 WRAP-UP -- two slim checks that decide how the manuscript phrases the
% two open items. Nothing is refitted: every (k, d1) is read from the existing
% step5_cell_fits CSV, the force law is untouched, and no data, exclusion or
% pipeline output is modified. The only file written is one working figure.
%
%   CHECK 1  ARCH SCALING. Is the residual arch DEPTH-locked (its zero crossing
%            sits at a fixed depth in the bed) or TRAJECTORY-locked (it sits at
%            a fixed fraction of that cell's final depth)? This decides one
%            clause in the departure list and guards against over-reading the
%            a_local coincidence found in 5.4a.
%
%   CHECK 2  TSTOP WINDOW SENSITIVITY. tstop5 uses an ABSOLUTE speed window
%            (4-30 cm/s). At low v0 that window is a large fraction of the
%            trajectory -- at v0 = 70 it opens at 0.43*v0, deep in the curved
%            part of v(t) -- so the linear extrapolation truncates the soft tail
%            hardest exactly where the decline lives. Re-running with a RELATIVE
%            window documents that with numbers instead of an argument.
%
% Base MATLAB only.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));  % shared KD helpers

% -- constants --------------------------------------------------------
mass = 65;      % projectile mass, g
grav = 980;     % cm/s^2

MASTER  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.mat';
EXPDIR  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports';
FIGPATH = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step5_wrapup_working.png';

% default mask, identical to checkpoints 5.1-5.4
Z_MIN = 0.1;    % cm
V_MIN = 15;     % cm/s

MDL = ["Tight" "Default" "Wide"];
MRK = {'o' 's' '^'};
CO  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

GRID_STEP = 0.025;   % cm, common depth grid inside a cell
NBOOT = 1000;        % cluster-bootstrap samples

% check-2 windows
ABS_LO = 4;  ABS_HI = 30;        % cm/s, the standard tstop5 window
REL_LO = 0.05; REL_HI = 0.25;    % fractions of v0, the relative alternative

fprintf('\n=== STEP 5 WRAP-UP CHECKS ===\n');
fprintf('No refitting: (k, d1) are read from the existing cell-fit export.\n');

%% ===================================================================
%  LOAD
%  ===================================================================
L = load(MASTER);
T = L.T;  S = L.S;

keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
model = string(K.model);
h_mm  = K.dropHeight_mm;

cellCsv = local_newest(EXPDIR, 'step5_cell_fits_*.csv');
C = readtable(cellCsv);  C.model = string(C.model);
fprintf('cell fits : %s (%d cells)\n', cellCsv, height(C));

% per-trial series and scalars
D = struct('model', cell(nTr,1));
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

[cellId, cellModel, cellH] = findgroups(model, h_mm);
nCell = max(cellId);

%% ===================================================================
%  CHECK 1 -- ARCH SCALING
%  ===================================================================
fprintf('\n=== 1. ARCH SCALING (all %d cells) ===\n', nCell);

z_cross = nan(nCell,1); d_final_mean = nan(nCell,1); v0_mean = nan(nCell,1);
for c = 1:nCell
    % this cell's already-fitted parameters
    r = find(C.model == cellModel(c) & C.h_mm == cellH(c), 1);
    if isempty(r) || ~isfinite(C.k_fit(r)), continue; end
    kc = C.k_fit(r); dc = C.d1_fit(r); v0_mean(c) = C.v0_mean(r);

    idx = find(cellId == c);
    d_final_mean(c) = mean([D(idx).d_final]);

    % ensemble-median speed2 curve, built exactly as in step5_fleet_fits
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
    resid = med(mk) - kd_speed2_model(zg(mk), v0_mean(c), kc, dc, mass, grav);

    % the arch's upward zero crossing, taken AFTER the shallow negative lobe
    z_cross(c) = local_zero_cross_after_lobe(zg(mk), resid);
end

ok = isfinite(z_cross) & isfinite(d_final_mean);
ratio = z_cross ./ d_final_mean;      % trajectory-locked would hold this fixed
fprintf('  z_cross recovered on %d / %d cells\n', sum(ok), nCell);

fprintf('\n  %-8s %6s %10s %10s %14s %14s\n', ...
    'model', 'cells', 'med z_cr', 'IQR z_cr', 'med z_cr/d_f', 'IQR z_cr/d_f');
for s = 1:4
    if s < 4, sel = ok & cellModel == MDL(s); else, sel = ok; end
    nm = local_tern(s < 4, char(MDL(min(s,3))), 'Pooled');
    fprintf('  %-8s %6d %10.3f %10.3f %14.3f %14.3f\n', nm, sum(sel), ...
        median(z_cross(sel)), local_iqr(z_cross(sel)), ...
        median(ratio(sel)),   local_iqr(ratio(sel)));
end

% slope of z_cross on d_final_mean, with a cluster bootstrap over cells
fprintf('\n  slope of z_cross on d_final_mean (cluster bootstrap, B = %d):\n', NBOOT);
rng(1);
for s = 1:4
    if s < 4, sel = find(ok & cellModel == MDL(s)); else, sel = find(ok); end
    nm = local_tern(s < 4, char(MDL(min(s,3))), 'Pooled');
    X = [ones(numel(sel),1) d_final_mean(sel)];
    b = X \ z_cross(sel);
    bs = nan(NBOOT,1);
    for q = 1:NBOOT
        p = sel(randi(numel(sel), numel(sel), 1));      % resample cells
        Xb = [ones(numel(p),1) d_final_mean(p)];
        cb = Xb \ z_cross(p);
        bs(q) = cb(2);
    end
    ci = local_prctile(bs, [2.5 97.5]);
    fprintf('    %-8s slope = %+.4f [%+.4f %+.4f] cm per cm\n', nm, b(2), ci(1), ci(2));
    if s == 4, slope_pool = b(2); ci_pool = ci; end
end

% decision, on the two spreads plus the pooled slope CI
iqr_abs = local_iqr(z_cross(ok));
rel_abs = iqr_abs / median(z_cross(ok));               % relative spread, depth-locked view
rel_rat = local_iqr(ratio(ok)) / median(ratio(ok));    % relative spread, trajectory-locked view
fprintf('\n  spreads: z_cross IQR = %.3f cm (relative %.3f) | z_cross/d_final relative %.3f\n', ...
    iqr_abs, rel_abs, rel_rat);
if ci_pool(1) <= 0 && ci_pool(2) >= 0 && iqr_abs < 0.2
    fprintf(['  DECISION: DEPTH-LOCKED. The pooled slope CI [%+.4f %+.4f] includes zero and the\n' ...
             '  crossing sits within an IQR of %.3f cm of a fixed depth, so the arch is tied to\n' ...
             '  position in the bed rather than to how far that trial travelled.\n'], ci_pool, iqr_abs);
elseif rel_rat < rel_abs
    fprintf(['  DECISION: TRAJECTORY-LOCKED. z_cross/d_final is relatively tighter (%.3f vs %.3f),\n' ...
             '  so the crossing tracks a fixed fraction of the trial''s own final depth.\n'], rel_rat, rel_abs);
else
    fprintf(['  DECISION: MIXED. The pooled slope is %+.4f [%+.4f %+.4f] cm per cm and neither\n' ...
             '  view is clearly tighter (depth %.3f vs trajectory %.3f), so the arch is not\n' ...
             '  cleanly locked to either depth or trajectory.\n'], slope_pool, ci_pool, rel_abs, rel_rat);
end

%% ===================================================================
%  CHECK 2 -- TSTOP ESTIMATOR SENSITIVITY
%  ===================================================================
fprintf('\n=== 2. TSTOP WINDOW SENSITIVITY (6 cells) ===\n');
fprintf('  Geometries: Tight and Wide (the ends of the ladder).\n');
fprintf('  (a) absolute window %g-%g cm/s   (b) relative window %.2f*v0-%.2f*v0\n', ...
    ABS_LO, ABS_HI, REL_LO, REL_HI);

pickG = [1 3];                       % Tight and Wide
pick = [];                           % the six chosen cells
for g = pickG
    sel = find(cellModel == MDL(g) & isfinite(v0_mean));
    [~, ord] = sort(v0_mean(sel));
    pick = [pick; sel(ord(1)); sel(ord(round(numel(ord)/2))); sel(ord(end))];  %#ok<AGROW>
end

fprintf('\n  %-8s %8s %10s %12s %12s %10s %14s\n', ...
    'model', 'v0_mean', 'trials', 'tstop abs', 'tstop rel', 'ratio', 'window (rel)');
tab = nan(6,3);
for q = 1:numel(pick)
    c = pick(q);
    idx = find(cellId == c);
    ta = nan(numel(idx),1); tr = nan(numel(idx),1);
    for i = 1:numel(idx)
        d = D(idx(i));
        % (a) the standard absolute window
        ta(i) = local_tstop(d, ABS_LO, ABS_HI);
        % (b) the same rule with the window scaled to this trial's own v0
        tr(i) = local_tstop(d, REL_LO*d.v0, REL_HI*d.v0);
    end
    tab(q,:) = [v0_mean(c), median(ta,'omitnan'), median(tr,'omitnan')];
    fprintf('  %-8s %8.1f %10d %12.5f %12.5f %10.3f %6.1f-%.1f cm/s\n', ...
        cellModel(c), v0_mean(c), numel(idx), tab(q,2), tab(q,3), tab(q,3)/tab(q,2), ...
        REL_LO*v0_mean(c), REL_HI*v0_mean(c));
end

% does the decline reappear once the window scales with v0?
fprintf('\n');
for gi = 1:2
    rows = (gi-1)*3 + (1:3);
    dAbs = 100*(tab(rows(3),2) - tab(rows(1),2)) / tab(rows(1),2);
    dRel = 100*(tab(rows(3),3) - tab(rows(1),3)) / tab(rows(1),3);
    fprintf('  %-8s low->high v0 change: absolute window %+.1f%% | relative window %+.1f%%\n', ...
        MDL(pickG(gi)), dAbs, dRel);
end
dAbsAll = mean([100*(tab(3,2)-tab(1,2))/tab(1,2), 100*(tab(6,2)-tab(4,2))/tab(4,2)]);
dRelAll = mean([100*(tab(3,3)-tab(1,3))/tab(1,3), 100*(tab(6,3)-tab(4,3))/tab(4,3)]);
if dRelAll < dAbsAll - 5
    fprintf(['\n  The declining trend REAPPEARS under the relative window (mean low->high change\n' ...
             '  %+.1f%% relative vs %+.1f%% absolute), confirming that the flat tstop5 trend is\n' ...
             '  window-induced at low v0 rather than a property of the trajectories.\n'], dRelAll, dAbsAll);
else
    fprintf(['\n  The decline does NOT reappear under the relative window (mean low->high change\n' ...
             '  %+.1f%% relative vs %+.1f%% absolute), so the flat tstop5 trend is not explained\n' ...
             '  by the absolute window alone.\n'], dRelAll, dAbsAll);
end

%% ===================================================================
%  3. FIGURE (working)
%  ===================================================================
figure('Position', [80 80 1300 500]);

% (a) arch crossing against the cell's mean final depth
ax = subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
hM = gobjects(1,3);
for g = 1:3
    sel = ok & cellModel == MDL(g);
    hM(g) = plot(ax, d_final_mean(sel), z_cross(sel), MRK{g}, 'MarkerSize', 6, 'Color', CO(g,:));
end
% the pooled fitted slope, drawn as a reference line
xs = linspace(min(d_final_mean(ok)), max(d_final_mean(ok)), 50);
Xp = [ones(sum(ok),1) d_final_mean(ok)];  bp = Xp \ z_cross(ok);
hF = plot(ax, xs, bp(1) + bp(2)*xs, 'k--', 'LineWidth', 1.2);
xlabel(ax,'cell mean d_{final} (cm)'); ylabel(ax,'arch zero crossing z_{cross} (cm)');
title(ax, sprintf('(a) arch scaling: slope %+.3f [%+.3f %+.3f]', bp(2), ci_pool));
legend(ax, [hM hF], [cellstr(MDL) {'pooled fit'}], 'Location','best');

% (b) the six-cell window comparison
ax = subplot(1,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
bar(ax, tab(:,2:3));
% one single-line label per group: an embedded newline would be split across ticks
lab = arrayfun(@(q) sprintf('%s  v_0=%.0f', cellModel(pick(q)), tab(q,1)), ...
               (1:6).', 'UniformOutput', false);
set(ax, 'XTick', 1:6, 'XTickLabel', lab, 'XTickLabelRotation', 30);
ylabel(ax,'median stopping time (s)');
title(ax,'(b) tstop window sensitivity');
legend(ax, {'absolute 4-30 cm/s','relative 0.05-0.25 v_0'}, 'Location','best');

exportgraphics(gcf, FIGPATH, 'Resolution', 200);
fprintf('\nfigure written: %s\n', FIGPATH);

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function t0 = local_tstop(d, vlo, vhi)
% The tstop5 rule with a configurable speed window: linear fit of v against t
% over vlo < speed < vhi AND depth > 0.5*d_final, extrapolated to v = 0.
    t0 = NaN;
    w = d.speed > vlo & d.speed < vhi & d.depth > 0.5*d.d_final & isfinite(d.speed);
    if sum(w) < 5, return; end
    p = polyfit(d.t_s(w), d.speed(w), 1);
    if p(1) ~= 0, t0 = -p(2)/p(1); end
end

function zc = local_zero_cross_after_lobe(zg, r)
% The arch's upward zero crossing: the LAST negative-to-positive crossing
% before the arch's positive peak, i.e. the crossing that leads into the peak.
%
% Anchoring on the peak rather than on the global minimum matters. In several
% cells the residual's deepest point is the terminal downturn near the stop,
% not the shallow lobe, so searching forward from the global minimum finds no
% upward crossing at all and the cell is silently lost. Anchoring on the peak
% locates the same feature in every cell that has an arch.
    zc = NaN;
    ok = isfinite(r); zg = zg(ok); r = r(ok);
    if numel(r) < 3, return; end
    [pk, jpk] = max(r);
    if pk <= 0 || jpk < 2, return; end                   % no positive arch at all
    j = find(r(1:jpk-1) < 0 & r(2:jpk) >= 0, 1, 'last'); % last crossing before the peak
    if isempty(j), return; end
    zc = zg(j) + (zg(j+1)-zg(j)) * (-r(j)) / (r(j+1)-r(j));   % linear interpolation
end

function [zz, uu] = local_monotone_descent(depth, speed, speed2)
% Strictly increasing descent segment, so interp1 gets unique sample points.
% Identical to step5_fleet_fits.
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

function p = local_newest(dirPath, pattern)
% Newest file matching a pattern, by modification time.
    f = dir(fullfile(dirPath, pattern));
    if isempty(f), error('step5w:noFile', 'No file matching %s in %s', pattern, dirPath); end
    [~, j] = max([f.datenum]);
    p = fullfile(f(j).folder, f(j).name);
end

function v = local_iqr(x)
% Interquartile range, base MATLAB.
    q = local_prctile(x, [25 75]);
    v = q(2) - q(1);
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

function s = local_tern(cond, a, b)
% Small inline conditional, purely to keep the printf lines readable.
    if cond, s = a; else, s = b; end
end
