% STEP 5.2 -- the published Katsuragi-Durian law fitted across the full fleet.
%
% Question: are k and d1 CONSTANTS of the substrate, as the published law
% asserts, or do they move with impact speed and foot geometry? The force law
% itself is not modified anywhere in this script. Nothing writes to the data,
% the exclusions or the pipeline outputs; the only files produced are two CSV
% exports and one working figure.
%
% INTERPRETATION HIERARCHY -- this governs how every number below is read:
%
%   PRIMARY     the 54 (geometry x height) CELL fits. Averaging within a cell
%               before fitting suppresses per-trial frame noise, so these are
%               the test of whether k and d1 are constants.
%   DIAGNOSTIC  the 524 PER-TRIAL fits. Distributions and failure flags only.
%               A per-trial median is not evidence about constancy; it is
%               evidence about the spread of single-trial fits.
%   SUMMARY     the GLOBAL fits. The trial-balanced estimator is primary and
%               the concatenated-point estimator is secondary, always labelled.
%
% The numerics -- kd_speed2_model, the nested solver (closed-form k inside,
% fminbnd over d1 outside) and the default mask z_min = 0.1 cm, v_min = 15 cm/s
% -- are the ones verified in scripts/step5_exemplar_fits.m and now shared from
% src/ so the two checkpoints cannot drift apart.
%
% Base MATLAB only. Runtime is several minutes.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));  % shared KD helpers

% -- constants --------------------------------------------------------
mass = 65;      % projectile mass, g   (as track_tracers_2 / export_master_dataset)
grav = 980;     % cm/s^2               (as get_calibration)
k_seed = 1.70e5;  % energy-argument seed for k, g/s^2, carried for comparison only

MASTER  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.mat';
OUTDIR  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports';
FIGPATH = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step5_fleet_fits_working.png';

% default mask, identical to checkpoint 5.1
Z_MIN = 0.1;    % cm
V_MIN = 15;     % cm/s

% geometry identity, fixed across every panel and table (project standard)
MDL = ["Tight" "Default" "Wide"];
MRK = {'o' 's' '^'};

GRID_STEP = 0.025;   % cm, the common depth grid inside a cell
NBOOT = 1000;        % cluster-bootstrap samples

% geometry landmarks drawn on the residual panels (reference lines ONLY,
% they are not cutoffs and nothing is masked at them)
LANDMARKS = [0.80 1.12];

stamp = datestr(now, 'yyyymmdd_HHMMSS');

fprintf('\n=== STEP 5.2  KD FLEET FITS ===\n');
fprintf('HIERARCHY: 54 cell fits are PRIMARY; 524 per-trial fits are DIAGNOSTIC only;\n');
fprintf('           global fits are SUMMARY (trial-balanced primary, concatenated secondary).\n');

%% ===================================================================
%  LOAD -- analysis set and per-trial series
%  ===================================================================
L = load(MASTER);
T = L.T;  S = L.S;

% the reviewed, non-quarantined set: this selection is not re-litigated here
keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
fprintf('\nanalysis set: %d trials (keep_reviewed & ~isZeroDrop)\n', nTr);

tagsS = string({S.trialTag});
model = string(K.model);
h_mm  = K.dropHeight_mm;

% one struct per trial: the series, the default mask, and the scalars
D = struct('tag', cell(nTr,1));
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    if isempty(j)
        error('step5b:missingSeries', 'No series for %s.', K.trialTag(i));
    end
    D(i).tag     = string(K.trialTag(i));
    D(i).model   = model(i);
    D(i).h_mm    = h_mm(i);
    D(i).depth   = S(j).z_cm(:);
    D(i).speed   = S(j).v_cm_s(:);
    D(i).speed2  = D(i).speed.^2;
    D(i).t_s     = S(j).t_s(:);
    D(i).v0      = K.v0_cm_s(i);
    D(i).d_final = K.d_final_cm(i);
    % the default mask, identical to 5.1
    D(i).mask    = D(i).depth > Z_MIN & isfinite(D(i).speed) & D(i).speed > V_MIN;
    D(i).nused   = sum(D(i).mask);
end

% cell identity: one cell per (geometry x height)
[cellId, cellModel, cellH] = findgroups(model, h_mm);
nCell = max(cellId);
fprintf('cells: %d (geometry x height)\n', nCell);

%% ===================================================================
%  A. PER-TRIAL FITS (524) -- DIAGNOSTIC ONLY
%  ===================================================================
fprintf('\n=== A. PER-TRIAL FITS (DIAGNOSTIC ONLY) ===\n');

k_fit_tr  = nan(nTr,1);  d1_fit_tr = nan(nTr,1);
rmse_tr   = nan(nTr,1);  dpred_tr  = nan(nTr,1);
kneg_tr   = false(nTr,1); bound_tr = false(nTr,1); ill_tr = false(nTr,1);

for i = 1:nTr
    z = D(i).depth(D(i).mask);  u = D(i).speed2(D(i).mask);
    if numel(z) < 5, continue; end            % too few frames to fit at all
    [kf, df, ~, rs] = kd_fit_nested(z, u, D(i).v0, mass, grav);
    k_fit_tr(i) = kf;  d1_fit_tr(i) = df;  rmse_tr(i) = rs;

    % flags: each is a fact about the fit, none of them change the fit
    kneg_tr(i)  = kf < 0;                                      % closed-form k < 0 at the optimum
    bound_tr(i) = df <= 0.1*1.01 || df >= 30*0.99;             % pinned within 1% of a bound
    ill_tr(i)   = local_ill_cond(df, z, u, D(i).v0, mass, grav);  % flat outer objective

    dpred_tr(i) = local_predict_d(D(i).v0, kf, df, mass, grav, D(i).d_final);

    if mod(i,100) == 0, fprintf('  fitted %d / %d trials\n', i, nTr); end
end
fprintf('  fitted %d / %d trials\n', nTr, nTr);

c_fit_tr = mass ./ d1_fit_tr;
flagged_tr = kneg_tr | bound_tr | ill_tr;

fprintf('\n  flag counts (of %d trials):\n', nTr);
fprintf('    k_neg       : %d\n', sum(kneg_tr));
fprintf('    d1_at_bound : %d\n', sum(bound_tr));
fprintf('    ill_cond    : %d\n', sum(ill_tr));
fprintf('    any flag    : %d   (clean %d)\n', sum(flagged_tr), sum(~flagged_tr));

% per-geometry medians and IQRs, flagged trials excluded and counted separately
fprintf('\n  per-geometry per-trial distributions (flagged excluded from medians):\n');
fprintf('  %-8s %6s %7s %12s %12s %10s %10s\n', ...
    'model', 'n', 'flagged', 'med k_fit', 'IQR k_fit', 'med d1', 'IQR d1');
for g = 1:3
    sel = model == MDL(g) & ~flagged_tr;
    nfl = sum(model == MDL(g) & flagged_tr);
    fprintf('  %-8s %6d %7d %12.4e %12.4e %10.3f %10.3f\n', MDL(g), sum(sel), nfl, ...
        median(k_fit_tr(sel)), local_iqr(k_fit_tr(sel)), ...
        median(d1_fit_tr(sel)), local_iqr(d1_fit_tr(sel)));
end

% export -- with a header comment that states what these rows are for
csvA = fullfile(OUTDIR, sprintf('step5_pertrial_fits_%s.csv', stamp));
A = table(string({D.tag}).', model, h_mm, [D.v0].', k_fit_tr, d1_fit_tr, c_fit_tr, ...
          rmse_tr, [D.nused].', dpred_tr, kneg_tr, bound_tr, ill_tr, ...
    'VariableNames', {'trialTag','model','dropHeight_mm','v0_cm_s','k_fit','d1_fit', ...
                      'c_fit','rmse_speed2','n_frames_used','d_pred','k_neg', ...
                      'd1_at_bound','ill_cond'});
fid = fopen(csvA, 'w');
fprintf(fid, '# step5_pertrial_fits: diagnostic distributions only; the 54 cell fits are the primary constancy test\n');
fclose(fid);
writetable(A, csvA, 'WriteMode', 'append', 'WriteVariableNames', true);
fprintf('\n  wrote %s\n', csvA);

%% ===================================================================
%  B. ENSEMBLE-MEDIAN CELL FITS (PRIMARY)
%  ===================================================================
fprintf('\n=== B. CELL FITS (PRIMARY, %d cells) ===\n', nCell);

k_fit_c = nan(nCell,1);  d1_fit_c = nan(nCell,1);  rmse_c = nan(nCell,1);
v0_c    = nan(nCell,1);  n_c      = nan(nCell,1);
kneg_c  = false(nCell,1); bound_c = false(nCell,1); ill_c = false(nCell,1);
resid_c = cell(nCell,1);  grid_c  = cell(nCell,1);

for c = 1:nCell
    idx = find(cellId == c);
    n_c(c)  = numel(idx);
    v0_c(c) = mean([D(idx).v0]);

    % common depth grid for the cell, out to its deepest trial
    zg = (0 : GRID_STEP : max([D(idx).d_final])).';
    M  = nan(numel(zg), numel(idx));
    for q = 1:numel(idx)
        % strictly increasing descent segment, so interp1 has unique samples
        [zz, uu] = local_monotone_descent(D(idx(q)).depth, D(idx(q)).speed, D(idx(q)).speed2);
        if numel(zz) < 2, continue; end
        M(:,q) = interp1(zz, uu, zg, 'linear', NaN);     % linear, no extrapolation
    end

    % median across the cell's trials, kept only where at least half support it
    support = sum(isfinite(M), 2);
    med = median(M, 2, 'omitnan');
    med(support < ceil(n_c(c)/2)) = NaN;

    % the default mask, applied to the median curve
    mk = zg > Z_MIN & isfinite(med) & med > V_MIN^2;
    if sum(mk) < 5, continue; end

    % the cell's mean v0 is the boundary condition for the median curve
    [kf, df, ~, rs] = kd_fit_nested(zg(mk), med(mk), v0_c(c), mass, grav);
    k_fit_c(c) = kf;  d1_fit_c(c) = df;  rmse_c(c) = rs;
    kneg_c(c)  = kf < 0;
    bound_c(c) = df <= 0.1*1.01 || df >= 30*0.99;
    ill_c(c)   = local_ill_cond(df, zg(mk), med(mk), v0_c(c), mass, grav);

    % residual curve for the arch analysis (section E)
    grid_c{c}  = zg(mk);
    resid_c{c} = med(mk) - kd_speed2_model(zg(mk), v0_c(c), kf, df, mass, grav);
end

c_fit_c = mass ./ d1_fit_c;
flagged_c = kneg_c | bound_c | ill_c;
fprintf('  cell flags: k_neg %d | d1_at_bound %d | ill_cond %d | any %d\n', ...
    sum(kneg_c), sum(bound_c), sum(ill_c), sum(flagged_c));

fprintf('\n  per-geometry CELL distributions (PRIMARY):\n');
fprintf('  %-8s %6s %12s %12s %10s %10s\n', 'model', 'cells', 'med k_fit', 'IQR k_fit', 'med d1', 'IQR d1');
for g = 1:3
    sel = cellModel == MDL(g) & isfinite(k_fit_c);
    fprintf('  %-8s %6d %12.4e %12.4e %10.3f %10.3f\n', MDL(g), sum(sel), ...
        median(k_fit_c(sel)), local_iqr(k_fit_c(sel)), ...
        median(d1_fit_c(sel)), local_iqr(d1_fit_c(sel)));
end

csvB = fullfile(OUTDIR, sprintf('step5_cell_fits_%s.csv', stamp));
Bt = table(cellModel, cellH, n_c, v0_c, k_fit_c, d1_fit_c, c_fit_c, rmse_c, ...
           kneg_c, bound_c, ill_c, ...
    'VariableNames', {'model','h_mm','n_trials','v0_mean','k_fit','d1_fit','c_fit', ...
                      'rmse_speed2','k_neg','d1_at_bound','ill_cond'});
fid = fopen(csvB, 'w');
fprintf(fid, '# step5_cell_fits: PRIMARY constancy test -- ensemble-median curve per (geometry x height) cell\n');
fclose(fid);
writetable(Bt, csvB, 'WriteMode', 'append', 'WriteVariableNames', true);
fprintf('  wrote %s\n', csvB);

%% ===================================================================
%  C. GLOBAL FITS (SUMMARY) -- trial-balanced primary, concatenated secondary
%  ===================================================================
fprintf('\n=== C. GLOBAL FITS (SUMMARY) ===\n');

% pool every masked point, carrying its own trial v0, weight and cell
z_all = []; u_all = []; v0_all = []; w_all = []; cell_all = []; geom_all = [];
for i = 1:nTr
    if D(i).nused < 5, continue; end
    m = D(i).mask;
    z_all   = [z_all;   D(i).depth(m)];                       %#ok<AGROW>
    u_all   = [u_all;   D(i).speed2(m)];                      %#ok<AGROW>
    v0_all  = [v0_all;  repmat(D(i).v0, D(i).nused, 1)];      %#ok<AGROW>
    % trial-balanced weight: every trial contributes total weight 1, so a
    % long record cannot outvote a short one. Frame count correlates with v0
    % and geometry, which is exactly the axis being tested.
    w_all   = [w_all;   repmat(1/D(i).nused, D(i).nused, 1)]; %#ok<AGROW>
    cell_all = [cell_all; repmat(cellId(i), D(i).nused, 1)];  %#ok<AGROW>
    geom_all = [geom_all; repmat(find(MDL == D(i).model), D(i).nused, 1)]; %#ok<AGROW>
end
fprintf('  pooled masked points: %d\n', numel(z_all));

% point indices per cell, so the bootstrap resamples cells cheaply
cellPts = arrayfun(@(c) find(cell_all == c), (1:nCell).', 'UniformOutput', false);

setName = [MDL "Pooled"];
kb = nan(4,1); db = nan(4,1); cb = nan(4,1);       % trial-balanced (PRIMARY)
kc = nan(4,1); dc = nan(4,1); cc = nan(4,1);       % concatenated  (secondary)
for s = 1:4
    if s < 4, sel = geom_all == s; else, sel = true(size(geom_all)); end
    % PRIMARY: weighted so every trial counts once
    [kb(s), db(s)] = kd_fit_nested(z_all(sel), u_all(sel), v0_all(sel), mass, grav, w_all(sel));
    cb(s) = mass / db(s);
    % SECONDARY: unweighted concatenation of points
    [kc(s), dc(s)] = kd_fit_nested(z_all(sel), u_all(sel), v0_all(sel), mass, grav);
    cc(s) = mass / dc(s);
end

% cluster bootstrap on the trial-balanced estimator only
fprintf('  cluster bootstrap (B = %d, resampling cells) ...\n', NBOOT);
rng(1);
kCI = nan(4,2); dCI = nan(4,2); cCI = nan(4,2);
for s = 1:4
    if s < 4, pool = find(cellModel == MDL(s)); else, pool = (1:nCell).'; end
    kbs = nan(NBOOT,1); dbs = nan(NBOOT,1);
    for b = 1:NBOOT
        pick = pool(randi(numel(pool), numel(pool), 1));
        idx  = vertcat(cellPts{pick});
        [kbs(b), dbs(b)] = kd_fit_nested(z_all(idx), u_all(idx), v0_all(idx), ...
                                         mass, grav, w_all(idx));
    end
    kCI(s,:) = local_prctile(kbs, [2.5 97.5]);
    dCI(s,:) = local_prctile(dbs, [2.5 97.5]);
    cCI(s,:) = local_prctile(mass./dbs, [2.5 97.5]);
    fprintf('    %s done\n', setName(s));
end

fprintf('\n  PRIMARY trial-balanced (cluster CIs from %d cell resamples):\n', NBOOT);
for s = 1:4
    fprintf('    %-8s k = %.4e [%.4e %.4e] | d1 = %7.3f [%7.3f %7.3f] cm | c = %8.3f [%8.3f %8.3f] g/cm\n', ...
        setName(s), kb(s), kCI(s,1), kCI(s,2), db(s), dCI(s,1), dCI(s,2), cb(s), cCI(s,1), cCI(s,2));
end
fprintf('\n  SECONDARY concatenated points (no CIs; secondary by construction):\n');
for s = 1:4
    fprintf('    %-8s k = %.4e | d1 = %7.3f cm | c = %8.3f g/cm\n', setName(s), kc(s), dc(s), cc(s));
end
% a disagreement here is frame-count / v0 confounding, not a physical result
fprintf('\n  balanced vs concatenated (trial-balanced stands where they differ):\n');
for s = 1:4
    fprintf('    %-8s k differs by %+6.1f%% | d1 differs by %+6.1f%%%s\n', setName(s), ...
        100*(kc(s)-kb(s))/kb(s), 100*(dc(s)-db(s))/db(s), ...
        local_tern(abs((dc(s)-db(s))/db(s)) > 0.10 || abs((kc(s)-kb(s))/kb(s)) > 0.10, ...
                   '   <- frame-count/v0 confounding', ''));
end

%% ===================================================================
%  D. CONSTANCY DIAGNOSTICS (HEADLINE) -- computed on the CELL fits
%  ===================================================================
fprintf('\n=== D. CONSTANCY DIAGNOSTICS (on the PRIMARY cell fits) ===\n');

% (i) Spearman rank correlation -- DESCRIPTIVE ONLY, never a decision criterion
fprintf('\n  (i) Spearman rank correlation vs v0_mean -- DESCRIPTIVE ONLY:\n');
for s = 1:4
    if s < 4, sel = cellModel == MDL(s) & isfinite(k_fit_c); else, sel = isfinite(k_fit_c); end
    fprintf('    %-8s rho(d1, v0) = %+.3f   rho(k, v0) = %+.3f   (n = %d cells)\n', ...
        setName(s), local_spearman(d1_fit_c(sel), v0_c(sel)), ...
        local_spearman(k_fit_c(sel), v0_c(sel)), sum(sel));
end

% (ii) effect size with uncertainty: high/low v0 terciles of the cell fits
fprintf('\n  (ii) tercile effect size, cluster-bootstrapped (B = %d):\n', NBOOT);
rng(1);
% terciles are formed WITHIN geometry; the pooled row combines those terciles
loIdx = cell(3,1); hiIdx = cell(3,1);
for g = 1:3
    sel = find(cellModel == MDL(g) & isfinite(k_fit_c));
    [~, ord] = sort(v0_c(sel));
    nt = floor(numel(sel)/3);
    loIdx{g} = sel(ord(1:nt));
    hiIdx{g} = sel(ord(end-nt+1:end));
end

fprintf('  %-8s %-6s %22s %12s %s\n', 'model', 'param', 'ratio high/low [95% CI]', 'noise floor', 'decision');
for s = 1:4
    if s < 4, gset = s; else, gset = 1:3; end
    lo = vertcat(loIdx{gset});  hi = vertcat(hiIdx{gset});
    for p = 1:2
        if p == 1, x = d1_fit_c; pname = 'd1'; else, x = k_fit_c; pname = 'k'; end
        ratio = median(x(hi)) / median(x(lo));
        % bootstrap resamples cells WITHIN each tercile, preserving the split
        rb = nan(NBOOT,1);
        for b = 1:NBOOT
            lb = lo(randi(numel(lo), numel(lo), 1));
            hb = hi(randi(numel(hi), numel(hi), 1));
            rb(b) = median(x(hb)) / median(x(lb));
        end
        ci = local_prctile(rb, [2.5 97.5]);
        % noise floor: the relative IQR of the cell fits within the terciles
        nf = mean([local_iqr(x(lo))/abs(median(x(lo))), local_iqr(x(hi))/abs(median(x(hi)))]);
        % a failure claim needs the CI to exclude 1 AND the effect to clear the floor
        if ci(1) <= 1 && ci(2) >= 1
            decision = 'consistent with constant';
        elseif abs(ratio-1) <= nf
            decision = sprintf('varies by factor %.2f [%.2f %.2f] but within the noise floor', ratio, ci);
        else
            decision = sprintf('varies by factor %.2f [%.2f %.2f]', ratio, ci);
        end
        fprintf('  %-8s %-6s %8.3f [%.3f %.3f] %11.3f   %s\n', ...
            setName(s), pname, ratio, ci(1), ci(2), nf, decision);
    end
end

%% ===================================================================
%  E. RESIDUAL ARCH -- does the 5.1 exemplar arch recur at fleet level?
%  ===================================================================
fprintf('\n=== E. RESIDUAL ARCH (cell residual curves) ===\n');

medResid = cell(3,1); medGrid = cell(3,1);
for g = 1:3
    sel = find(cellModel == MDL(g) & isfinite(k_fit_c));
    zmax = max(cellfun(@(v) max(v), grid_c(sel)));
    zg = (0 : GRID_STEP : zmax).';
    R = nan(numel(zg), numel(sel));
    for q = 1:numel(sel)
        c = sel(q);
        % normalised by v0^2 so cells of different impact speed are comparable
        R(:,q) = interp1(grid_c{c}, resid_c{c}/v0_c(c)^2, zg, 'linear', NaN);
    end
    medGrid{g}  = zg;
    medResid{g} = median(R, 2, 'omitnan');

    % where the median arch crosses zero and where it peaks
    [zc, zp, pk] = local_arch_features(zg, medResid{g});
    fprintf('  %-8s zero crossing at %.3f cm | positive peak at %.3f cm (%.4f of v0^2)\n', ...
        MDL(g), zc, zp, pk);
end

%% ===================================================================
%  F. STEP5-DERIVED a_stop AND t_stop (pipeline outputs untouched)
%  ===================================================================
fprintf('\n=== F. STEP5-DERIVED a_stop / t_stop ===\n');
fprintf('  Pipeline a_stop is NaN for 590/600 trials, so these are recomputed HERE and\n');
fprintf('  are labelled step5-derived. Nothing in the pipeline is modified.\n');

astop5 = nan(nTr,1); tstop5 = nan(nTr,1);
for i = 1:nTr
    % window: past half depth, speed between 4 and 30 cm/s, at least 5 frames
    w = D(i).speed > 4 & D(i).speed < 30 & D(i).depth > 0.5*D(i).d_final & isfinite(D(i).speed);
    if sum(w) < 5, continue; end
    p = polyfit(D(i).t_s(w), D(i).speed(w), 1);
    astop5(i) = p(1);                       % signed slope, expected negative
    if p(1) ~= 0, tstop5(i) = -p(2)/p(1); end  % linear extrapolation to v = 0
end
k_astop5 = mass*(grav - astop5) ./ [D.d_final].';

fprintf('\n  %-8s %9s %7s %14s %12s %14s %14s\n', ...
    'model', 'recovered', 'NaN', 'med astop5', 'med tstop5', 'med k_astop5', 'IQR k_astop5');
for g = 1:3
    sel = model == MDL(g);
    ok  = sel & isfinite(astop5);
    fprintf('  %-8s %9d %7d %14.1f %12.5f %14.4e %14.4e\n', MDL(g), sum(ok), sum(sel & ~isfinite(astop5)), ...
        median(astop5(ok)), median(tstop5(ok), 'omitnan'), ...
        median(k_astop5(ok)), local_iqr(k_astop5(ok)));
end
fprintf('\n  three-k comparison per geometry (cell-fit k is the primary one):\n');
fprintf('  %-8s %14s %14s %14s\n', 'model', 'cell-fit k', 'k_astop5', 'k_seed');
for g = 1:3
    selc = cellModel == MDL(g) & isfinite(k_fit_c);
    ok   = model == MDL(g) & isfinite(astop5);
    fprintf('  %-8s %14.4e %14.4e %14.4e\n', MDL(g), median(k_fit_c(selc)), median(k_astop5(ok)), k_seed);
end

%% ===================================================================
%  G. FIGURE (working, minimal styling)
%  ===================================================================
figure('Position', [60 60 1500 900]);

% (a) cell-fit k against cell mean v0
ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for g = 1:3
    sel = cellModel == MDL(g) & isfinite(k_fit_c);
    plot(ax, v0_c(sel), k_fit_c(sel), MRK{g}, 'MarkerSize', 6);
end
xlabel(ax,'cell mean v_0 (cm/s)'); ylabel(ax,'k_{fit} (g/s^2)');
title(ax,'(a) cell fits: k vs v_0'); legend(ax, cellstr(MDL), 'Location','best');

% (b) cell-fit d1 against cell mean v0, log y because d1 spans decades
ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for g = 1:3
    sel = cellModel == MDL(g) & isfinite(d1_fit_c);
    plot(ax, v0_c(sel), d1_fit_c(sel), MRK{g}, 'MarkerSize', 6);
end
set(ax,'YScale','log');
% widen the limits so no cell is clipped at the axis edge
ylim(ax, [0.8*min(d1_fit_c(isfinite(d1_fit_c))), 1.25*max(d1_fit_c(isfinite(d1_fit_c)))]);
xlabel(ax,'cell mean v_0 (cm/s)'); ylabel(ax,'d_1 (cm)');
title(ax,'(b) cell fits: d_1 vs v_0'); legend(ax, cellstr(MDL), 'Location','best');

% (c) per-geometry median residual curves, with the geometry landmarks
ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
yline(ax, 0, '-');
for zL = LANDMARKS, xline(ax, zL, ':', sprintf('%.2f cm', zL)); end
hCurve = gobjects(1,3);
for g = 1:3
    hCurve(g) = plot(ax, medGrid{g}, medResid{g}, '-', 'LineWidth', 1.5);
end
xlabel(ax,'depth (cm)'); ylabel(ax,'median residual speed^2 / v_0^2');
% legend the curves explicitly, or it would label the landmark lines instead
title(ax,'(c) median residual arch'); legend(ax, hCurve, cellstr(MDL), 'Location','best');

% (d) per-trial d1 distribution by geometry, log x
ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
edges = logspace(log10(0.1), log10(30), 30);
for g = 1:3
    sel = model == MDL(g) & isfinite(d1_fit_tr);
    histogram(ax, d1_fit_tr(sel), edges, 'DisplayStyle','stairs', 'LineWidth', 1.5);
end
set(ax,'XScale','log');
xlabel(ax,'per-trial d_1 (cm)'); ylabel(ax,'trials');
title(ax,'(d) per-trial d_1 (DIAGNOSTIC)'); legend(ax, cellstr(MDL), 'Location','best');

exportgraphics(gcf, FIGPATH, 'Resolution', 200);
fprintf('\nfigure written: %s\n', FIGPATH);
fprintf('\nREMINDER: cell fits are PRIMARY; per-trial fits are diagnostic; globals are summaries.\n');

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function [zz, uu] = local_monotone_descent(depth, speed, speed2)
% The descent segment with strictly increasing depth, so interp1 gets unique
% sample points. Runs from peak speed (impact) to deepest point.
    [~, i1] = max(speed);
    [~, i2] = max(depth);
    zz = []; uu = [];
    if i2 <= i1, return; end
    z = depth(i1:i2); u = speed2(i1:i2);
    ok = isfinite(z) & isfinite(u);
    z = z(ok); u = u(ok);
    % keep each point only if it is deeper than every point kept before it
    keep = false(numel(z),1); last = -inf;
    for q = 1:numel(z)
        if z(q) > last, keep(q) = true; last = z(q); end
    end
    zz = z(keep); uu = u(keep);
end

function tf = local_ill_cond(d1, z, u, v0, mass, grav)
% Flat outer objective: RSS changes by less than 1% over a 2x range of d1
% centred on the optimum, evaluated on a coarse grid.
    lo = max(0.1, d1/sqrt(2));  hi = min(30, d1*sqrt(2));
    dg = linspace(lo, hi, 9);
    rg = arrayfun(@(dd) local_rss_at(dd, z, u, v0, mass, grav), dg);
    tf = (max(rg) - min(rg)) / min(rg) < 0.01;
end

function r = local_rss_at(d1, z, u, v0, mass, grav)
% RSS after the inner closed-form k solve at a fixed d1.
    [~, r] = kd_inner_k(d1, z, u, v0, mass, grav);
end

function d_pred = local_predict_d(v0, k, d1, mass, grav, d_final)
% Predicted final depth = the root of kd_speed2_model(z) = 0.
    f = @(z) kd_speed2_model(z, v0, k, d1, mass, grav);
    zs = linspace(1e-6, max(6*d_final, 1), 400);
    fs = arrayfun(f, zs);
    j  = find(fs(1:end-1) > 0 & fs(2:end) <= 0, 1);
    try
        if ~isempty(j), d_pred = fzero(f, [zs(j) zs(j+1)]);
        else,           d_pred = fzero(f, d_final);
        end
    catch
        d_pred = NaN;                       % model never reaches zero speed
    end
end

function [zc, zp, pk] = local_arch_features(zg, r)
% Depth of the median curve's first negative-to-positive zero crossing and of
% its positive peak.
    zc = NaN; zp = NaN; pk = NaN;
    ok = isfinite(r);
    zg = zg(ok); r = r(ok);
    if numel(r) < 3, return; end
    j = find(r(1:end-1) < 0 & r(2:end) >= 0, 1);
    if ~isempty(j)
        zc = zg(j) + (zg(j+1)-zg(j)) * (-r(j)) / (r(j+1)-r(j));   % linear interpolation
    end
    [pk, jp] = max(r);
    zp = zg(jp);
end

function v = local_iqr(x)
% Interquartile range, base MATLAB (no Statistics Toolbox).
    q = local_prctile(x, [25 75]);
    v = q(2) - q(1);
end

function y = local_prctile(x, p)
% Percentiles by linear interpolation of the empirical CDF at midpoints,
% matching MATLAB's prctile convention. Base MATLAB.
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0, y = nan(size(p)); return; end
    if n == 1, y = repmat(x, size(p)); return; end
    q = (0.5:1:(n-0.5)) / n * 100;
    y = interp1(q, x, p, 'linear', 'extrap');
    y = min(max(y, x(1)), x(end));           % clamp to the observed range
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

function s = local_tern(cond, a, b)
% Small inline conditional, purely to keep the printf lines readable.
    if cond, s = a; else, s = b; end
end
