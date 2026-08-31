% STEP 5.5 -- does the geometry ladder need different (k, d1)?
%
% Formal test of whether Tight / Default / Wide require separate KD parameters.
% No new force-law terms, no change to the mask, the data or the exclusions.
% The model, the nested solver and the trial-balanced weighting are the ones
% used in step5_fleet_fits, shared from src/ so nothing can drift.
%
%   A  CELL-LEVEL COMPARISON (PRIMARY). Cells are the independent units --
%      18 per geometry -- so the pairwise ratios and their bootstrap CIs are
%      the test. Per-trial fits are not used here.
%   B  NESTED LADDER on the trial-balanced objective: shared, per-geometry k,
%      per-geometry d1, fully separate. Selected by cluster-bootstrap BIC.
%   C  OBSERVABLE CONSEQUENCE. Even if parameters differ, does any predicted
%      depth move by more than the measurement noise floor?
%   D  CONFOUND STATEMENT.
%
% Base MATLAB only.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));  % shared KD helpers

% -- constants --------------------------------------------------------
mass = 65;      % projectile mass, g
grav = 980;     % cm/s^2

MASTER  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.mat';
EXPDIR  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports';
FIGPATH = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step5_geometry_ladder_working.png';

% default mask, identical to checkpoints 5.1-5.4
Z_MIN = 0.1;    % cm
V_MIN = 15;     % cm/s

MDL = ["Tight" "Default" "Wide"];
MRK = {'o' 's' '^'};
CO  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

NB_RATIO = 2000;    % bootstrap samples for the cell-level ratios
NB_MODEL = 500;     % cluster-bootstrap samples for model selection
SD_WITHIN = 0.326;  % cm, pooled within-(model,height) SD of d_final (step4)
V0_TABLE = [100 175 250];

fprintf('\n=== STEP 5.5  GEOMETRY LADDER ===\n');
fprintf('Cells are the independent units: the cell-level ratios in A are the PRIMARY test.\n');

%% ===================================================================
%  LOAD
%  ===================================================================
L = load(MASTER);
T = L.T;  S = L.S;
keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
model = string(K.model);

cellCsv = local_newest(EXPDIR, 'step5_cell_fits_*.csv');
C = readtable(cellCsv);  C.model = string(C.model);
fprintf('\nanalysis set %d trials | cell fits %s (%d cells)\n', nTr, cellCsv, height(C));

%% ===================================================================
%  A. CELL-LEVEL COMPARISON (PRIMARY)
%  ===================================================================
fprintf('\n=== A. CELL-LEVEL COMPARISON (PRIMARY, 18 cells per geometry) ===\n');

ok = isfinite(C.k_fit) & isfinite(C.d1_fit);
fprintf('\n  %-8s %6s %14s %14s %10s %10s\n', 'model', 'cells', 'med k', 'IQR k', 'med d1', 'IQR d1');
for g = 1:3
    sel = ok & C.model == MDL(g);
    fprintf('  %-8s %6d %14.4e %14.4e %10.3f %10.3f\n', MDL(g), sum(sel), ...
        median(C.k_fit(sel)), local_iqr(C.k_fit(sel)), ...
        median(C.d1_fit(sel)), local_iqr(C.d1_fit(sel)));
end

% pairwise ratios, bootstrap resampling cells WITHIN each geometry
pairs = [3 2; 3 1; 1 2];                 % W/D, W/T, T/D  (indices into MDL)
pairName = ["W/D" "W/T" "T/D"];
fprintf('\n  pairwise ratios, cluster bootstrap over cells (B = %d):\n', NB_RATIO);
rng(1);
for p = 1:2
    if p == 1, x = C.k_fit; pname = 'k'; else, x = C.d1_fit; pname = 'd1'; end
    for q = 1:3
        ia = find(ok & C.model == MDL(pairs(q,1)));
        ib = find(ok & C.model == MDL(pairs(q,2)));
        ratio = median(x(ia)) / median(x(ib));
        rb = nan(NB_RATIO,1);
        for b = 1:NB_RATIO
            sa = ia(randi(numel(ia), numel(ia), 1));
            sb = ib(randi(numel(ib), numel(ib), 1));
            rb(b) = median(x(sa)) / median(x(sb));
        end
        ci = local_prctile(rb, [2.5 97.5]);
        % retained only so the distribution can be persisted below; no draw is
        % added or reordered here, so every printed number is unchanged
        ratioBoot{p,q} = rb; ratioCI(p,q,:) = ci; ratioPt(p,q) = ratio; %#ok<SAGROW>
        ratioName(p,q) = pname + " " + pairName(q);                     %#ok<SAGROW>
        if ci(1) <= 1 && ci(2) >= 1
            verdict = 'CI includes 1';
        else
            verdict = sprintf('differs by factor %.3f [%.3f %.3f]', ratio, ci);
        end
        fprintf('    %-3s %-3s ratio = %.3f [%.3f %.3f]   %s\n', pname, pairName(q), ratio, ci, verdict);
    end
end

% Kruskal-Wallis across the three geometries -- DESCRIPTIVE label only
fprintf('\n  Kruskal-Wallis across geometries (DESCRIPTIVE, not the decision criterion):\n');
gidx = double(C.model == MDL(1)) + 2*double(C.model == MDL(2)) + 3*double(C.model == MDL(3));
[Hk, pk] = local_kruskal_wallis(C.k_fit(ok),  gidx(ok));
[Hd, pd] = local_kruskal_wallis(C.d1_fit(ok), gidx(ok));
fprintf('    k  : H = %.3f, p = %.4f\n', Hk, pk);
fprintf('    d1 : H = %.3f, p = %.4f\n', Hd, pd);

%% ===================================================================
%  B. NESTED LADDER on the trial-balanced objective
%  ===================================================================
fprintf('\n=== B. NESTED LADDER (trial-balanced weighting) ===\n');

% pool every masked point, carrying its trial v0, geometry, weight and cell
tagsS = string({S.trialTag});
[cellId, ~, ~] = findgroups(model, K.dropHeight_mm);
nCell = max(cellId);
z_all = []; u_all = []; v0_all = []; w_all = []; g_all = []; c_all = [];
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    depth = S(j).z_cm(:); speed = S(j).v_cm_s(:);
    m = depth > Z_MIN & isfinite(speed) & speed > V_MIN;
    n = sum(m);
    if n < 5, continue; end
    z_all  = [z_all;  depth(m)];                          %#ok<AGROW>
    u_all  = [u_all;  speed(m).^2];                       %#ok<AGROW>
    v0_all = [v0_all; repmat(K.v0_cm_s(i), n, 1)];        %#ok<AGROW>
    w_all  = [w_all;  repmat(1/n, n, 1)];                 %#ok<AGROW>
    g_all  = [g_all;  repmat(find(MDL == model(i)), n, 1)];%#ok<AGROW>
    c_all  = [c_all;  repmat(cellId(i), n, 1)];           %#ok<AGROW>
end
n_eff = sum(w_all);       % = number of contributing trials, by construction
fprintf('  pooled points %d | n_eff = %.1f trials (weights sum to 1 per trial)\n', numel(z_all), n_eff);

% fit the four models on the full data
Mfit = local_fit_ladder(z_all, u_all, v0_all, w_all, g_all, mass, grav);

% convergence check for M2, whose outer search is multidimensional
fprintf('\n  M2 convergence check (restarts from perturbed seeds):\n');
rng(2);
base = Mfit.M2.d1;
for r = 1:3
    seed = base .* (0.6 + 0.8*rand(1,3));
    alt = local_fit_M2(z_all, u_all, v0_all, w_all, g_all, mass, grav, seed);
    fprintf('    seed [%.2f %.2f %.2f] -> d1 [%.3f %.3f %.3f], rss %.6e (base %.6e, rel diff %.2e)\n', ...
        seed, alt.d1, alt.rss, Mfit.M2.rss, abs(alt.rss - Mfit.M2.rss)/Mfit.M2.rss);
end

% BIC with n_eff = 524: the weights sum to one per trial, so a TRIAL is the
% observation unit, not a frame. This is an approximation, stated as such.
names = ["M0 shared" "M1 k per geom" "M2 d1 per geom" "M3 separate"];
npar  = [2 4 4 6];
flds  = ["M0" "M1" "M2" "M3"];
fprintf('\n  BIC uses n_eff = %.0f (weights sum to 1 per trial; trials treated as\n', n_eff);
fprintf('  the observations, frames are not independent within a trial). Approximate.\n');
fprintf('\n  %-16s %5s %14s %14s\n', 'model', 'par', 'wRMSE', 'BIC');
bic = nan(1,4);
for m = 1:4
    F = Mfit.(flds(m));
    bic(m) = n_eff*log(F.rss/n_eff) + npar(m)*log(n_eff);
    fprintf('  %-16s %5d %14.2f %14.2f\n', names(m), npar(m), sqrt(F.rss/n_eff), bic(m));
end
[~, bestFull] = min(bic);
fprintf('  lowest BIC on the full data: %s\n', names(bestFull));

% robust criterion: cluster-bootstrap model selection over the 54 cells
fprintf('\n  cluster-bootstrap model selection (B = %d, resampling cells):\n', NB_MODEL);
cellPts = arrayfun(@(c) find(c_all == c), (1:nCell).', 'UniformOutput', false);
rng(1);
wins = zeros(1,4);
for b = 1:NB_MODEL
    pick = randi(nCell, nCell, 1);
    idx  = vertcat(cellPts{pick});
    zb = z_all(idx); ub = u_all(idx); vb = v0_all(idx); wb = w_all(idx); gb = g_all(idx);
    if numel(unique(gb)) < 3, continue; end        % need all three geometries present
    nb = sum(wb);
    Fb = local_fit_ladder(zb, ub, vb, wb, gb, mass, grav);
    bb = nan(1,4);
    for m = 1:4
        bb(m) = nb*log(Fb.(flds(m)).rss/nb) + npar(m)*log(nb);
    end
    [~, jw] = min(bb);
    wins(jw) = wins(jw) + 1;
    if mod(b, 100) == 0, fprintf('    %d / %d resamples\n', b, NB_MODEL); end
end
frac = wins / sum(wins);
fprintf('\n  %-16s %12s\n', 'model', 'win fraction');
for m = 1:4
    fprintf('  %-16s %12.3f\n', names(m), frac(m));
end
[fbest, jbest] = max(frac);
fprintf('  DECISION: %s wins %.1f%% of cluster resamples.\n', names(jbest), 100*fbest);

%% ===================================================================
%  C. OBSERVABLE-LEVEL CONSEQUENCE
%  ===================================================================
fprintf('\n=== C. OBSERVABLE CONSEQUENCE (does any parameter difference matter?) ===\n');
fprintf('\n  d_pred (cm) at three speeds, M0 shared vs M3 per-geometry:\n');
fprintf('  %-8s %6s %10s %10s %10s\n', 'model', 'v0', 'M0', 'M3', 'diff');
maxdiff = 0;
for g = 1:3
    for q = 1:numel(V0_TABLE)
        v0 = V0_TABLE(q);
        d0 = local_predict_d(v0, Mfit.M0.k(1),  Mfit.M0.d1(1),  mass, grav);
        d3 = local_predict_d(v0, Mfit.M3.k(g),  Mfit.M3.d1(g),  mass, grav);
        fprintf('  %-8s %6d %10.3f %10.3f %+10.3f\n', MDL(g), v0, d0, d3, d3-d0);
        maxdiff = max(maxdiff, abs(d3-d0));
    end
end
fprintf('\n  largest |M3 - M0| = %.3f cm = %.2f x the within-group SD of %.3f cm\n', ...
    maxdiff, maxdiff/SD_WITHIN, SD_WITHIN);
if maxdiff < SD_WITHIN
    fprintf(['  Per-geometry parameters move no predicted depth by as much as one within-group\n' ...
             '  SD, so the parameter difference has no observable consequence at the noise floor.\n']);
else
    fprintf(['  Per-geometry parameters move at least one predicted depth by more than the\n' ...
             '  within-group SD, so the difference is observable above the noise floor.\n']);
end

%% ===================================================================
%  D. CONFOUND STATEMENT
%  ===================================================================
fprintf('\n=== D. CONFOUND ===\n');
fprintf(['  Geometry is aliased with capture campaign and bed handling in this dataset\n' ...
         '  (separate campaigns per model, documented release-height offsets). Any\n' ...
         '  between-geometry parameter difference is an upper bound on a splay effect;\n' ...
         '  it cannot be attributed to splay alone.\n']);

%% ===================================================================
%  E. FIGURE (working)
%  ===================================================================
figure('Position', [60 60 1450 900]);

% (a) cell k by geometry
ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for g = 1:3
    sel = ok & C.model == MDL(g);
    xj = g + 0.14*(rand(sum(sel),1)-0.5);        % small jitter so points separate
    plot(ax, xj, C.k_fit(sel), MRK{g}, 'MarkerSize', 5, 'Color', CO(g,:));
    plot(ax, g+[-0.25 0.25], median(C.k_fit(sel))*[1 1], 'k-', 'LineWidth', 2);
end
set(ax,'XTick',1:3,'XTickLabel',cellstr(MDL)); xlim(ax,[0.5 3.5]);
ylabel(ax,'cell k_{fit} (g/s^2)'); title(ax,'(a) cell k by geometry');

% (b) cell d1 by geometry
ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for g = 1:3
    sel = ok & C.model == MDL(g);
    xj = g + 0.14*(rand(sum(sel),1)-0.5);
    plot(ax, xj, C.d1_fit(sel), MRK{g}, 'MarkerSize', 5, 'Color', CO(g,:));
    plot(ax, g+[-0.25 0.25], median(C.d1_fit(sel))*[1 1], 'k-', 'LineWidth', 2);
end
set(ax,'XTick',1:3,'XTickLabel',cellstr(MDL)); xlim(ax,[0.5 3.5]);
ylabel(ax,'cell d_1 (cm)'); title(ax,'(b) cell d_1 by geometry');

% (c) joint (k, d1): shows which direction the two parameters trade off
ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
hG = gobjects(1,3);
for g = 1:3
    sel = ok & C.model == MDL(g);
    hG(g) = plot(ax, C.d1_fit(sel), C.k_fit(sel), MRK{g}, 'MarkerSize', 5, 'Color', CO(g,:));
end
xlabel(ax,'cell d_1 (cm)'); ylabel(ax,'cell k_{fit} (g/s^2)');
title(ax,'(c) joint (k, d_1): trade-off orientation');
legend(ax, hG, cellstr(MDL), 'Location','best');

% (d) predicted depth under the shared and the separate model
ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
vg = 60:5:290;
d0c = arrayfun(@(v) local_predict_d(v, Mfit.M0.k(1), Mfit.M0.d1(1), mass, grav), vg);
h0 = plot(ax, vg, d0c, 'k-', 'LineWidth', 2);
hG = gobjects(1,3);
for g = 1:3
    d3c = arrayfun(@(v) local_predict_d(v, Mfit.M3.k(g), Mfit.M3.d1(g), mass, grav), vg);
    hG(g) = plot(ax, vg, d3c, '-', 'LineWidth', 1.3, 'Color', CO(g,:));
end
xlabel(ax,'v_0 (cm/s)'); ylabel(ax,'d_{pred} (cm)');
title(ax,'(d) d_{pred}: M0 shared (black) vs M3 separate');
legend(ax, [h0 hG], [{'M0 shared'} cellstr(MDL)], 'Location','northwest');

exportgraphics(gcf, FIGPATH, 'Resolution', 200);
fprintf('\nfigure written: %s\n', FIGPATH);

%% ===================================================================
%  PERSIST THE RESAMPLE DISTRIBUTION (added for fig_dynamics rev2)
%  ===================================================================
% Downstream figures need the per-geometry MEDIAN distribution, not just the
% pairwise ratios, so it is generated here with this script's own scheme --
% resample the geometry's cells with replacement, take the median -- and its
% own seed. This block runs AFTER every print above and re-seeds explicitly,
% so it cannot perturb any number this script reports.
rng(1);
medBoot = struct('k', nan(NB_RATIO,3), 'c', nan(NB_RATIO,3), 'd1', nan(NB_RATIO,3));
for g = 1:3
    ig = find(ok & C.model == MDL(g));
    for b = 1:NB_RATIO
        s = ig(randi(numel(ig), numel(ig), 1));
        medBoot.k(b,g)  = median(C.k_fit(s));
        medBoot.d1(b,g) = median(C.d1_fit(s));
        % c is formed per cell as mass/d1 and medianed after, as the figure requires
        medBoot.c(b,g)  = median(mass ./ C.d1_fit(s));
    end
end
LADDER_BOOT = struct('scheme', 'resample cells within geometry, with replacement', ...
    'B', NB_RATIO, 'seed', 1, 'models', MDL, 'pairName', pairName, ...
    'ratioBoot', {ratioBoot}, 'ratioCI', ratioCI, 'ratioPt', ratioPt, ...
    'ratioName', ratioName, 'medBoot', medBoot, ...
    'cellCsv', cellCsv, 'mass', mass);
bootPath = fullfile(EXPDIR, 'step5_ladder_bootstrap.mat');
save(bootPath, 'LADDER_BOOT');
fprintf('resample distribution written: %s\n', bootPath);

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function F = local_fit_ladder(z, u, v0, w, g, mass, grav)
% All four models of the ladder, on the trial-balanced objective.
%   M0 shared (k, d1)        M1 per-geometry k, shared d1
%   M2 shared k, per-geom d1 M3 fully separate
    % M0: one d1, one k
    o0 = @(d1) local_rss_given(z, u, v0, w, g, d1, true, mass, grav);
    d0 = fminbnd(o0, 0.1, 30);
    [r0, k0] = local_rss_given(z, u, v0, w, g, d0, true, mass, grav);
    F.M0 = struct('d1', d0, 'k', k0, 'rss', r0);

    % M1: one d1, three k
    o1 = @(d1) local_rss_given(z, u, v0, w, g, d1, false, mass, grav);
    d1s = fminbnd(o1, 0.1, 30);
    [r1, k1] = local_rss_given(z, u, v0, w, g, d1s, false, mass, grav);
    F.M1 = struct('d1', d1s, 'k', k1, 'rss', r1);

    % M2: three d1, one k -- multidimensional outer search
    F.M2 = local_fit_M2(z, u, v0, w, g, mass, grav, d0*[1 1 1]);

    % M3: fully separate, which is just three independent weighted fits
    k3 = nan(1,3); d3 = nan(1,3); r3 = 0;
    for q = 1:3
        s = g == q;
        [kq, dq, rq] = kd_fit_nested(z(s), u(s), v0(s), mass, grav, w(s));
        k3(q) = kq; d3(q) = dq; r3 = r3 + rq;
    end
    F.M3 = struct('d1', d3, 'k', k3, 'rss', r3);
end

function F = local_fit_M2(z, u, v0, w, g, mass, grav, seed)
% Shared k with a per-geometry d1. The outer search is over three lengths, so
% fminsearch is used; log-parameterising keeps every d1 positive without bounds.
    obj = @(L) local_rss_given(z, u, v0, w, g, exp(L), true, mass, grav);
    Lo  = fminsearch(obj, log(seed), optimset('TolX',1e-6,'TolFun',1e-8,'MaxFunEvals',2000));
    d2  = exp(Lo);
    [r2, k2] = local_rss_given(z, u, v0, w, g, d2, true, mass, grav);
    F = struct('d1', d2, 'k', k2, 'rss', r2);
end

function [rss, ks] = local_rss_given(z, u, v0, w, g, d1, sharedK, mass, grav)
% Weighted RSS with k solved in closed form at the given d1 (or d1 triple).
% d1 scalar -> one length for every point; d1 of length 3 -> per geometry.
% sharedK true -> a single k pooled over all points; false -> one k per group.
    if numel(d1) == 1, d1p = repmat(d1, size(z)); else, d1p = d1(g).'; d1p = d1p(:); end
    Ez = exp(-2*z ./ d1p);
    A  = v0.^2 .* Ez + grav*d1p .* (1-Ez);                       % k-independent part
    B  = -( (d1p/mass) .* z - (d1p.^2/(2*mass)) .* (1-Ez) );     % coefficient of k
    ks = nan(1,3);
    if sharedK
        k = sum(w .* B .* (u - A)) / sum(w .* B.^2);
        ks(:) = k;
        res = u - A - k*B;
    else
        res = nan(size(u));
        for q = 1:3
            s = g == q;
            if ~any(s), continue; end
            kq = sum(w(s) .* B(s) .* (u(s) - A(s))) / sum(w(s) .* B(s).^2);
            ks(q) = kq;
            res(s) = u(s) - A(s) - kq*B(s);
        end
    end
    rss = sum(w .* res.^2, 'omitnan');
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

function [H, p] = local_kruskal_wallis(x, g)
% Kruskal-Wallis H with the standard tie correction. For 3 groups the null
% distribution is chi-square with 2 df, whose survival function is exactly
% exp(-H/2) -- so no Statistics Toolbox is needed.
    x = x(:); g = g(:);
    ok = isfinite(x); x = x(ok); g = g(ok);
    N = numel(x);
    r = local_tiedrank(x);
    grp = unique(g);
    Hs = 0;
    for q = 1:numel(grp)
        s = g == grp(q);
        Hs = Hs + numel(x(s)) * (mean(r(s)) - (N+1)/2)^2;
    end
    H = 12/(N*(N+1)) * Hs;
    % tie correction
    [xs, ~] = sort(x);
    tie = 0; i = 1;
    while i <= N
        j = i;
        while j < N && xs(j+1) == xs(i), j = j + 1; end
        t = j - i + 1;
        if t > 1, tie = tie + t^3 - t; end
        i = j + 1;
    end
    if tie > 0, H = H / (1 - tie/(N^3 - N)); end
    p = exp(-H/2);                       % chi-square, 2 df
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

function p = local_newest(dirPath, pattern)
% Newest file matching a pattern, by modification time.
    f = dir(fullfile(dirPath, pattern));
    if isempty(f), error('step5g:noFile', 'No file matching %s in %s', pattern, dirPath); end
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
