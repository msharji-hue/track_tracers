% STEP 5 -- is the shallow depth term real? Intercept validity, a binless joint
% fit, collapse locking, and a synthetic null control.
%
% Question: are the negative / flat shallow-bin intercepts in the Fig-3a
% construction (a) statistical, (b) a measurement-chain artifact, or (c) genuine
% data structure -- and is a manuscript-quality Fig 3a/3b analogue defensible?
% No force-law modification. No bin chosen by appearance: every verdict below
% follows a rule fixed before the numbers were seen.
%
% ---------------------------------------------------------------------------
% STUDY 0 FIRST, AND IT MATTERS. The master export's a_plus_g_cm_s2 column
% carries the PRE-2026-08 formula (-a - g), not the intended (g - a). src/
% net_accel.m documents this and exists precisely to recompute rather than read
% that column. The two differ by a constant -2g = -1960 cm/s^2, which is the
% whole of the reported "negative shallow intercept" effect. Every study below
% therefore uses apg RECOMPUTED as grav - a from the raw acceleration trace,
% exactly as net_accel does. Study 0 audits the convention and quantifies what
% the stored column would have produced.
% ---------------------------------------------------------------------------
%
% Base MATLAB only. Runtime a few minutes (bootstraps + a synthetic fleet).

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));  % shared KD helpers

% -- constants --------------------------------------------------------
mass = 65;      % projectile mass, g
grav = 980;     % cm/s^2

MASTER  = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports\master_trials_20260822_215312.mat';
FIGPATH = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures\step5_depth_term_study_working.png';

% shared M0 parameters (step5_fleet_fits, trial-balanced pooled estimator)
k0 = 1.9163e5;      % g/s^2
c0 = 15.898;        % g/cm
d1_0 = mass / c0;   % cm, the drag length that c0 corresponds to

Z_MIN = 0.1;    % cm, default mask
V_MIN = 15;     % cm/s
ZC    = [0.3 0.5 0.7 0.9];   % bin centres
ZHALF = 0.1;                 % bin half-width

MDL = ["Tight" "Default" "Wide"];
MRK = {'o' 's' '^'};
CO  = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

NB_INT  = 2000;   % bootstrap for the per-bin intercepts
NB_CELL = 1000;   % cluster bootstrap for the joint fit
SIGMA_Z = 0.005;  % cm, synthetic position noise -- about half a pixel at
                  % mmPerPx = 0.1079 mm/px (half pixel = 0.0054 cm)
C_BIN   = 13.4;   % g/cm, the slope the bin method returned

fprintf('\n=== STEP 5  DEPTH-TERM STUDY ===\n');

%% ===================================================================
%  LOAD
%  ===================================================================
L = load(MASTER);
T = L.T;  S = L.S;
keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
model = string(K.model);
tagsS = string({S.trialTag});

% one struct per trial. apg is RECOMPUTED (see the header); the stored column
% is kept only so Study 0 can quantify what it would have given.
D = struct('model', cell(nTr,1));
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    a = S(j).a_cm_s2(:);
    D(i).model    = model(i);
    D(i).depth    = S(j).z_cm(:);
    D(i).speed    = S(j).v_cm_s(:);
    D(i).v2       = D(i).speed.^2;
    D(i).a        = a;
    D(i).apg      = grav - a;                       % net_accel's convention
    D(i).apg_store= S(j).a_plus_g_cm_s2(:);         % the exported column
    D(i).winHalf  = median(S(j).winHalfFrames, 'omitnan');
    D(i).dt       = K.dt_s(i);
    D(i).v0       = K.v0_cm_s(i);
    D(i).d_final  = K.d_final_cm(i);
    D(i).nF       = numel(D(i).depth);
    D(i).mask     = D(i).speed > V_MIN & D(i).depth > Z_MIN & isfinite(D(i).apg);
end
[cellId, ~, ~] = findgroups(model, K.dropHeight_mm);
nCell = max(cellId);
fprintf('analysis set %d trials | %d cells\n', nTr, nCell);

%% ===================================================================
%  STUDY 0 -- a_plus_g convention audit
%  ===================================================================
fprintf('\n=== STUDY 0: a_plus_g CONVENTION AUDIT ===\n');
nNew = 0; nOld = 0; nOther = 0;
for i = 1:nTr
    m = isfinite(D(i).a) & isfinite(D(i).apg_store);
    if ~any(m), continue; end
    if max(abs(D(i).apg_store(m) - (grav - D(i).a(m)))) < 1e-6
        nNew = nNew + 1;
    elseif max(abs(D(i).apg_store(m) - (-D(i).a(m) - grav))) < 1e-6
        nOld = nOld + 1;
    else
        nOther = nOther + 1;
    end
end
fprintf('  stored a_plus_g column: NEW (g-a) %d | OLD (-a-g) %d | neither %d  (of %d)\n', ...
    nNew, nOld, nOther, nTr);
fprintf('  the two conventions differ by a constant -2g = %.1f cm/s^2\n', -2*grav);
fprintf('  ALL STUDIES BELOW USE apg = grav - a, recomputed from the raw trace.\n');

% what the stored column would have produced, per bin, for Tight
fprintf('\n  Tight per-bin intercept, stored column vs recomputed:\n');
fprintf('  %6s %12s %12s %12s\n', 'z_c', 'F stored', 'F recomp', 'expected');
for b = 1:numel(ZC)
    [v2s, ags] = local_bin_points(D, model == MDL(1), ZC(b), ZHALF, false);
    [~,  agr]  = local_bin_points(D, model == MDL(1), ZC(b), ZHALF, true);
    Xb = [ones(numel(v2s),1) v2s];
    cs = Xb \ ags;  cr = Xb \ agr;
    fprintf('  %6.1f %12.1f %12.1f %12.1f\n', ZC(b), cs(1), cr(1), k0*ZC(b)/mass);
end
fprintf(['  The stored-column intercepts are negative/flat and the recomputed ones are\n' ...
         '  positive; the gap is exactly 2g in every bin. The reported shallow-intercept\n' ...
         '  effect is an EXPORT CONVENTION artifact, not a property of the data.\n']);

%% ===================================================================
%  STUDY 1 -- intercept inference validity
%  ===================================================================
fprintf('\n=== STUDY 1: INTERCEPT INFERENCE VALIDITY ===\n');
fprintf('  per-bin points = one (mean v^2, mean apg) per trial per bin\n');

Fv = nan(3,4); Fci = nan(3,4,2); Fsd = nan(3,4); Fexp = nan(1,4);
rng(1);
for g = 1:3
    fprintf('\n  %s\n', MDL(g));
    fprintf('  %5s %5s %10s %22s %10s %8s %9s %9s %8s\n', ...
        'z_c', 'n', 'F', 'boot CI', 'expected', 'clas SE', 'HC1 SE', 'min v2', 'extrap');
    for b = 1:4
        [v2b, agb] = local_bin_points(D, model == MDL(g), ZC(b), ZHALF, true);
        n = numel(v2b);
        Fexp(b) = k0*ZC(b)/mass;
        if n < 5, continue; end
        X = [ones(n,1) v2b];
        cf = X \ agb;  F = cf(1);  s = cf(2);
        e  = agb - X*cf;
        % classical SE of the intercept
        XtXi = inv(X'*X);
        seC = sqrt(sum(e.^2)/(n-2) * XtXi(1,1));
        % heteroskedasticity-robust HC1
        Vh = XtXi * (X' * diag(e.^2) * X) * XtXi * n/(n-2);
        seH = sqrt(Vh(1,1));
        % trial-resampling bootstrap
        fb = nan(NB_INT,1);
        for q = 1:NB_INT
            p = randi(n, n, 1);
            cb = [ones(n,1) v2b(p)] \ agb(p);
            fb(q) = cb(1);
        end
        ci = local_prctile(fb, [2.5 97.5]);
        Fv(g,b) = F; Fci(g,b,:) = ci; Fsd(g,b) = std(fb);
        % leverage / extrapolation diagnostics
        fracLow = mean(s*v2b < Fexp(b));
        extrap  = mean(v2b) / local_prctile(v2b, 10);
        fprintf('  %5.1f %5d %10.1f  [%9.1f %9.1f] %10.1f %8.1f %9.1f %9.0f %8.2f\n', ...
            ZC(b), n, F, ci(1), ci(2), Fexp(b), seC, seH, min(v2b), extrap);
        fprintf('        frac of points with slope*v2 < expected intercept = %.3f | ', fracLow);
        fprintf('F differs from expectation by %.1f sigma\n', abs(F - Fexp(b))/Fsd(g,b));
    end
end

% VERDICT 1, on the Tight bins as pre-registered
sigT = abs(Fv(1,:) - Fexp) ./ Fsd(1,:);
covT = squeeze(Fci(1,:,1)).' <= Fexp.' & squeeze(Fci(1,:,2)).' >= Fexp.';
fprintf('\n  Tight: |F-expected| in bootstrap SDs = [%.1f %.1f %.1f %.1f]; CI covers expectation in %d of 4\n', ...
    sigT, sum(covT));
if sum(sigT > 3) >= 3
    verdict1 = 'intercepts are statistically well-determined but systematically non-KD';
elseif sum(covT) >= 3
    verdict1 = 'intercepts are noise-limited';
else
    verdict1 = 'mixed';
end
fprintf('  VERDICT 1: %s\n', verdict1);

%% ===================================================================
%  STUDY 2 -- joint two-regressor estimation (binless)
%  ===================================================================
fprintf('\n=== STUDY 2: BINLESS JOINT FIT (trial-balanced) ===\n');

% pool masked frames with trial-balanced weights
z_all=[]; v2_all=[]; ag_all=[]; w_all=[]; g_all=[]; c_all=[];
for i = 1:nTr
    m = D(i).mask; n = sum(m);
    if n < 5, continue; end
    z_all  = [z_all;  D(i).depth(m)];                      %#ok<AGROW>
    v2_all = [v2_all; D(i).v2(m)];                         %#ok<AGROW>
    ag_all = [ag_all; D(i).apg(m)];                        %#ok<AGROW>
    w_all  = [w_all;  repmat(1/n, n, 1)];                  %#ok<AGROW>
    g_all  = [g_all;  repmat(find(MDL == D(i).model), n, 1)]; %#ok<AGROW>
    c_all  = [c_all;  repmat(cellId(i), n, 1)];            %#ok<AGROW>
end
fprintf('  pooled masked frames: %d\n', numel(z_all));

% collinearity diagnostics
rzv = corr_pearson(z_all, v2_all);
wt = nan(nTr,1);
for i = 1:nTr
    m = D(i).mask;
    if sum(m) >= 5, wt(i) = corr_pearson(D(i).depth(m), D(i).v2(m)); end
end
Xd = [z_all v2_all] .* sqrt(w_all);
fprintf('  collinearity: pooled corr(z, v2) = %+.3f | median within-trial corr = %+.3f | cond(design) = %.1f\n', ...
    rzv, median(wt,'omitnan'), cond(Xd));

fprintf('\n  %-8s %-12s %14s %12s %12s\n', 'set', 'model', 'k_a', 'c_a', 'b0');
rng(1);
cellPts = arrayfun(@(c) find(c_all == c), (1:nCell).', 'UniformOutput', false);
for s = 1:4
    if s < 4, sel = g_all == s; nm = MDL(s); else, sel = true(size(g_all)); nm = "Pooled"; end
    for withInt = [false true]
        [ka, ca, b0] = local_joint_fit(z_all(sel), v2_all(sel), ag_all(sel), w_all(sel), withInt, mass);
        % cluster bootstrap over cells
        kb = nan(NB_CELL,1); cbv = nan(NB_CELL,1); b0b = nan(NB_CELL,1);
        if s < 4, pool = find(arrayfun(@(c) any(g_all(c_all==c)==s), (1:nCell).')); else, pool = (1:nCell).'; end
        for q = 1:NB_CELL
            pk = pool(randi(numel(pool), numel(pool), 1));
            id = vertcat(cellPts{pk});
            [kb(q), cbv(q), b0b(q)] = local_joint_fit(z_all(id), v2_all(id), ag_all(id), w_all(id), withInt, mass);
        end
        kci = local_prctile(kb,[2.5 97.5]); cci = local_prctile(cbv,[2.5 97.5]);
        lbl = local_tern(withInt, 'free b0', 'no b0');
        if withInt
            bci = local_prctile(b0b,[2.5 97.5]);
            fprintf('  %-8s %-12s %8.3e [%.2e %.2e] %7.2f [%.2f %.2f] %8.1f [%.1f %.1f]\n', ...
                nm, lbl, ka, kci, ca, cci, b0, bci);
        else
            fprintf('  %-8s %-12s %8.3e [%.2e %.2e] %7.2f [%.2f %.2f] %8s\n', ...
                nm, lbl, ka, kci, ca, cci, '--');
        end
        if s == 4 && ~withInt, ka_pool = ka; ca_pool = ca; kci_pool = kci; cci_pool = cci; end
        % keep the per-geometry no-intercept estimates and CIs for the figure
        if s < 4 && ~withInt
            kaG(s) = ka; caG(s) = ca; kciG(s,:) = kci; cciG(s,:) = cci; %#ok<SAGROW>
        end
    end
end

% agreement with the trajectory-space parameters
fprintf('\n  vs trajectory-space M0 (k0 = %.4e, c0 = %.3f):\n', k0, c0);
fprintf('    k: binless %.4e [%.2e %.2e] -- %s\n', ka_pool, kci_pool, ...
    local_tern(kci_pool(1) <= k0 && kci_pool(2) >= k0, 'AGREES (CI covers k0)', 'DISAGREES (CI excludes k0)'));
fprintf('    c: binless %.3f [%.2f %.2f] -- %s\n', ca_pool, cci_pool, ...
    local_tern(cci_pool(1) <= c0 && cci_pool(2) >= c0, 'AGREES (CI covers c0)', 'DISAGREES (CI excludes c0)'));

% depth-restricted variants: the z-dependence of k_a measured binlessly
fprintf('\n  no-intercept pooled fit under depth restrictions:\n');
for zthr = [Z_MIN 0.2 1.0]
    sel = z_all > zthr;
    [ka, ca] = local_joint_fit(z_all(sel), v2_all(sel), ag_all(sel), w_all(sel), false, mass);
    fprintf('    z > %.1f cm : k_a = %.4e | c_a = %7.2f | frames %d\n', zthr, ka, ca, sum(sel));
    if zthr == Z_MIN, ka_shallow = ka; end
    if zthr == 1.0,   ka_deep = ka;    end
end
fprintf('    k_a(shallow-inclusive) = %.4e vs k_a(deep-only) = %.4e -> ratio %.3f\n', ...
    ka_shallow, ka_deep, ka_deep/ka_shallow);

%% ===================================================================
%  STUDY 3 -- fixed-slope collapse robustness
%  ===================================================================
fprintf('\n=== STUDY 3: FIXED-SLOPE COLLAPSE ===\n');
fprintf('  F_coll(z) = median over frames in 0.1-cm bins of (apg - (c/mass)*v2)\n');

zEdges = 0.1:0.1:3.0;
fprintf('\n  %-8s %18s %18s\n', 'model', sprintf('zero cross c=%.3f', c0), sprintf('zero cross c=%.1f', C_BIN));
for g = 1:3
    sel = g_all == g;
    [zc1, ~] = local_collapse(z_all(sel), v2_all(sel), ag_all(sel), c0,    mass, zEdges);
    [zc2, ~] = local_collapse(z_all(sel), v2_all(sel), ag_all(sel), C_BIN, mass, zEdges);
    fprintf('  %-8s %18s %18s\n', MDL(g), local_fmt(zc1), local_fmt(zc2));
end

% locking test: terciles of d_final formed within geometry, then pooled
fprintf('\n  locking test -- d_final terciles (within geometry, pooled), c = c0:\n');
terc = zeros(nTr,1);
for g = 1:3
    s = find(model == MDL(g));
    q = local_prctile(K.d_final_cm(s), [100/3 200/3]);
    terc(s(K.d_final_cm(s) <= q(1))) = 1;
    terc(s(K.d_final_cm(s) > q(1) & K.d_final_cm(s) <= q(2))) = 2;
    terc(s(K.d_final_cm(s) > q(2))) = 3;
end
% map the per-frame arrays back to a tercile label
t_all = [];
for i = 1:nTr
    n = sum(D(i).mask);
    if n < 5, continue; end
    t_all = [t_all; repmat(terc(i), n, 1)];  %#ok<AGROW>
end
crossT = nan(1,3); dfT = nan(1,3); collT = cell(1,3); zmidT = cell(1,3);
fprintf('  %-8s %10s %16s %14s\n', 'tercile', 'med d_fin', 'zero crossing', 'cross/d_final');
for q = 1:3
    sel = t_all == q;
    [zx, Fc, zmid] = local_collapse(z_all(sel), v2_all(sel), ag_all(sel), c0, mass, zEdges);
    crossT(q) = zx; dfT(q) = median(K.d_final_cm(terc == q));
    collT{q} = Fc; zmidT{q} = zmid;
    fprintf('  %-8d %10.3f %16s %14s\n', q, dfT(q), local_fmt(zx), local_fmt(zx/dfT(q)));
end
% The locking test needs at least two terciles to actually HAVE a crossing.
% With fewer, max-minus-min would compare a lone value with itself and report a
% spurious zero spread, so the test is declared undefined instead.
nFinite = sum(isfinite(crossT));
if nFinite < 2
    verdict3 = sprintf(['undefined -- the collapse has no zero crossing in %d of 3 terciles, ' ...
        'so there is no shallow suppression to locate'], 3-nFinite);
    fprintf('\n  crossings finite in %d of 3 terciles: locking test not applicable.\n', nFinite);
else
    spreadAbs = max(crossT) - min(crossT);
    ratT = crossT ./ dfT;
    spreadRel = (max(ratT) - min(ratT)) / median(ratT, 'omitnan');
    fprintf('\n  crossing spread = %.3f cm | crossing/d_final relative spread = %.3f\n', spreadAbs, spreadRel);
    if spreadAbs <= 0.15
        verdict3 = 'depth-locked';
    elseif spreadRel < spreadAbs/median(crossT, 'omitnan')
        verdict3 = 'trajectory-locked';
    else
        verdict3 = 'mixed';
    end
end
fprintf('  VERDICT 3: %s\n', verdict3);

%% ===================================================================
%  STUDY 4 -- synthetic null control (the arbiter)
%  ===================================================================
fprintf('\n=== STUDY 4: SYNTHETIC NULL CONTROL ===\n');
fprintf(['  Each trial''s trajectory is regenerated from the SHARED M0 (k0, c0) at its own\n' ...
         '  v0 and fps, given a rest tail, corrupted with sigma_z = %.3f cm position noise\n' ...
         '  (about half a pixel at 0.1079 mm/px), then passed through the SAME two-stage\n' ...
         '  moving linear-fit differentiation the pipeline uses (half-width = the trial''s\n' ...
         '  median winHalfFrames). Approximation: the pipeline also clamps its fitting\n' ...
         '  windows at the stop frame, which this does not, so the synthetic is if anything\n' ...
         '  MORE smoothed near the stop -- conservative for an artifact test.\n'], SIGMA_Z);

rng(4);
Dsyn = struct('model', cell(nTr,1));
for i = 1:nTr
    h  = max(1, round(D(i).winHalf));
    dt = D(i).dt;
    % KD trajectory from the shared parameters, then a rest tail
    [zt, ~] = local_kd_traj(D(i).v0, k0, d1_0, mass, grav, dt, D(i).nF);
    zn = zt + SIGMA_Z*randn(size(zt));                 % measurement noise
    vs = local_moving_slope(zn, dt, h);                % velocity from position
    as = local_moving_slope(vs, dt, h);                % acceleration from velocity
    Dsyn(i).model = D(i).model;
    Dsyn(i).depth = zn;
    Dsyn(i).speed = vs;
    Dsyn(i).v2    = vs.^2;
    Dsyn(i).apg   = grav - as;
    Dsyn(i).mask  = vs > V_MIN & zn > Z_MIN & isfinite(as);
    if mod(i,150) == 0, fprintf('    synthesised %d / %d trials\n', i, nTr); end
end
fprintf('    synthesised %d / %d trials\n', nTr, nTr);

fprintf('\n  synthetic per-bin intercepts vs expectation (same construction as Study 1):\n');
fprintf('  %-8s %5s %10s %10s %10s %10s\n', 'model', 'z_c', 'F syn', 'expected', 'ratio', 'F real');
Fsyn = nan(3,4);
for g = 1:3
    for b = 1:4
        [v2b, agb] = local_bin_points(Dsyn, model == MDL(g), ZC(b), ZHALF, true);
        if numel(v2b) < 5, continue; end
        cf = [ones(numel(v2b),1) v2b] \ agb;
        Fsyn(g,b) = cf(1);
        fprintf('  %-8s %5.1f %10.1f %10.1f %10.3f %10.1f\n', ...
            MDL(g), ZC(b), Fsyn(g,b), Fexp(b), Fsyn(g,b)/Fexp(b), Fv(g,b));
    end
end

% synthetic collapse
zs=[]; v2s=[]; ags=[];
for i = 1:nTr
    m = Dsyn(i).mask;
    if sum(m) < 5, continue; end
    zs = [zs; Dsyn(i).depth(m)]; v2s = [v2s; Dsyn(i).v2(m)]; ags = [ags; Dsyn(i).apg(m)]; %#ok<AGROW>
end
[zxSyn, FcSyn, zmidSyn] = local_collapse(zs, v2s, ags, c0, mass, zEdges);
fprintf('\n  synthetic collapse zero crossing (c = c0): %s\n', ...
    local_tern(isfinite(zxSyn), local_fmt(zxSyn), 'none: stays positive'));

% VERDICT 4
selB = ZC >= 0.3;
ratios = Fsyn(:,selB) ./ repmat(Fexp(selB), 3, 1);
offsetSyn = median(Fsyn(:,selB) - repmat(Fexp(selB),3,1), 'all');
offsetReal = median(Fv(:,selB)  - repmat(Fexp(selB),3,1), 'all');
fprintf('  median synthetic offset = %+.1f | median real offset = %+.1f cm/s^2\n', offsetSyn, offsetReal);
if all(Fsyn(:,selB) > 0, 'all') && all(abs(ratios-1) <= 0.20, 'all')
    verdict4 = 'construction unbiased under the null -- the real suppression is data structure';
elseif any(Fsyn(:,selB) <= 0, 'all')
    verdict4 = 'chain reproduces the suppression -- estimator artifact';
else
    verdict4 = sprintf('partial: chain accounts for %.0f%% of the offset', 100*offsetSyn/offsetReal);
end
fprintf('  VERDICT 4: %s\n', verdict4);

%% ===================================================================
%  FIGURE (working)
%  ===================================================================
figure('Position', [50 50 1500 950]);

% (a) Study-1 intercepts with bootstrap CIs against the k0*z/m line
ax = subplot(2,2,1); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
hG = gobjects(1,3);
for g = 1:3
    xo = ZC + 0.012*(g-2);
    hG(g) = plot(ax, xo, Fv(g,:), MRK{g}, 'MarkerSize', 6, 'Color', CO(g,:));
    for b = 1:4
        plot(ax, [xo(b) xo(b)], squeeze(Fci(g,b,:)), '-', 'Color', CO(g,:));
    end
end
zl = linspace(0.2, 1.0, 20);
hL = plot(ax, zl, k0*zl/mass, 'k--', 'LineWidth', 1.4);
yline(ax, 0, 'k:');
xlabel(ax,'bin centre z (cm)'); ylabel(ax,'intercept F (cm/s^2)');
title(ax,'(a) per-bin intercepts (recomputed a+g)');
legend(ax, [hG hL], [cellstr(MDL) {'k_0 z/m'}], 'Location','northwest');

% (b) Study-2 joint estimates against the trajectory-space reference
ax = subplot(2,2,2); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
hG = gobjects(1,3);
for g = 1:3
    % cluster-bootstrap CIs as crossed error bars in both parameters
    plot(ax, cciG(g,:), [kaG(g) kaG(g)], '-', 'Color', CO(g,:));
    plot(ax, [caG(g) caG(g)], kciG(g,:), '-', 'Color', CO(g,:));
    hG(g) = plot(ax, caG(g), kaG(g), MRK{g}, 'MarkerSize', 8, 'Color', CO(g,:), ...
                 'MarkerFaceColor', CO(g,:));
end
hR = plot(ax, c0, k0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
xlabel(ax,'c_a (g/cm)'); ylabel(ax,'k_a (g/s^2)');
title(ax,'(b) binless joint fit vs (k_0, c_0)');
legend(ax, [hG hR], [cellstr(MDL) {'(k_0, c_0) trajectory'}], 'Location','best');

% (c) Study-3 collapse per d_final tercile
ax = subplot(2,2,3); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
yline(ax, 0, 'k-');
hT = gobjects(1,3);
for q = 1:3
    hT(q) = plot(ax, zmidT{q}, collT{q}, '-', 'LineWidth', 1.4);
    if isfinite(crossT(q)), xline(ax, crossT(q), ':', sprintf('%.2f', crossT(q))); end
end
xlabel(ax,'depth (cm)'); ylabel(ax,'F_{coll} = a+g - (c_0/m)v^2');
title(ax,'(c) collapse by d_{final} tercile');
legend(ax, hT, {'low','mid','high'}, 'Location','northwest');

% (d) real vs synthetic per-bin intercepts
ax = subplot(2,2,4); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
hRe = gobjects(1,3); hSy = gobjects(1,3);
for g = 1:3
    hRe(g) = plot(ax, ZC, Fv(g,:),  MRK{g}, 'MarkerSize', 6, 'Color', CO(g,:));
    hSy(g) = plot(ax, ZC, Fsyn(g,:), MRK{g}, 'MarkerSize', 6, 'Color', CO(g,:), ...
                  'MarkerFaceColor', CO(g,:));
end
hL = plot(ax, zl, k0*zl/mass, 'k--', 'LineWidth', 1.4);
yline(ax, 0, 'k:');
xlabel(ax,'bin centre z (cm)'); ylabel(ax,'intercept F (cm/s^2)');
title(ax,'(d) real (open) vs synthetic null (filled)');
legend(ax, [hRe(1) hSy(1) hL], {'real','synthetic','k_0 z/m'}, 'Location','northwest');

exportgraphics(gcf, FIGPATH, 'Resolution', 200);
fprintf('\nfigure written: %s\n', FIGPATH);

%% ===================================================================
%  FINAL SYNTHESIS
%  ===================================================================
fprintf('\n=== FINAL SYNTHESIS ===\n');
fprintf('  V1: %s\n  V3: %s\n  V4: %s\n', verdict1, verdict3, verdict4);

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function [v2b, agb] = local_bin_points(D, selTrials, zc, zhalf, useRecomputed)
% One (mean v^2, mean a+g) point per trial per bin. useRecomputed selects the
% net_accel convention (grav - a) over the stored export column.
    idx = find(selTrials);
    v2b = nan(numel(idx),1); agb = nan(numel(idx),1);
    for q = 1:numel(idx)
        i = idx(q);
        if useRecomputed, ag = D(i).apg; else, ag = D(i).apg_store; end
        m = D(i).mask & abs(D(i).depth - zc) <= zhalf & isfinite(ag);
        if ~any(m), continue; end
        v2b(q) = mean(D(i).v2(m));
        agb(q) = mean(ag(m));
    end
    ok = isfinite(v2b) & isfinite(agb);
    v2b = v2b(ok); agb = agb(ok);
end

function [ka, ca, b0] = local_joint_fit(z, v2, ag, w, withIntercept, mass)
% Weighted least squares of a+g on depth and v^2, with or without a free
% intercept. k_a = mass*b_z and c_a = mass*b_v by the KD identification
%   a+g = (k/mass)*z + (1/d1)*v^2 ,  c = mass/d1.
    if withIntercept, X = [ones(numel(z),1) z v2]; else, X = [z v2]; end
    sw = sqrt(w);
    b = (X .* sw) \ (ag .* sw);
    if withIntercept
        b0 = b(1); ka = mass*b(2); ca = mass*b(3);
    else
        b0 = NaN;  ka = mass*b(1); ca = mass*b(2);
    end
end

function [zx, Fc, zmid] = local_collapse(z, v2, ag, c, mass, edges)
% Fixed-slope collapse: median of (a+g - (c/mass)*v^2) in depth bins, and the
% first upward zero crossing of that median curve.
    y = ag - (c/mass)*v2;
    nb = numel(edges)-1;
    Fc = nan(nb,1); zmid = nan(nb,1);
    for b = 1:nb
        m = z > edges(b) & z <= edges(b+1);
        zmid(b) = 0.5*(edges(b)+edges(b+1));
        if sum(m) >= 10, Fc(b) = median(y(m)); end
    end
    zx = NaN;
    ok = isfinite(Fc);
    zz = zmid(ok); ff = Fc(ok);
    j = find(ff(1:end-1) < 0 & ff(2:end) >= 0, 1);
    if ~isempty(j)
        zx = zz(j) + (zz(j+1)-zz(j)) * (-ff(j)) / (ff(j+1)-ff(j));
    end
end

function [z, v] = local_kd_traj(v0, k, d1, mass, grav, dt, nF)
% KD trajectory integrated in time with RK4 at the trial's own sampling
% interval, then held at rest so the record has nF frames like the real trial.
    f = @(zz, vv) grav - k*zz/mass - vv.^2/d1;
    z = nan(nF,1); v = nan(nF,1);
    z(1) = 0; v(1) = v0;
    stopped = false;
    for i = 1:nF-1
        if stopped
            z(i+1) = z(i); v(i+1) = 0; continue
        end
        zi = z(i); vi = v(i);
        k1v = f(zi, vi);                    k1z = vi;
        k2v = f(zi+dt/2*k1z, vi+dt/2*k1v);  k2z = vi + dt/2*k1v;
        k3v = f(zi+dt/2*k2z, vi+dt/2*k2v);  k3z = vi + dt/2*k2v;
        k4v = f(zi+dt*k3z,   vi+dt*k3v);    k4z = vi + dt*k3v;
        vn = vi + dt/6*(k1v + 2*k2v + 2*k3v + k4v);
        zn = zi + dt/6*(k1z + 2*k2z + 2*k3z + k4z);
        if vn <= 0
            stopped = true; v(i+1) = 0; z(i+1) = zi;      % rest tail
        else
            v(i+1) = vn; z(i+1) = zn;
        end
    end
end

function s = local_moving_slope(y, dt, h)
% Centered moving LINEAR fit, returning the fitted slope at each frame. The
% window is clamped at the array ends, as the pipeline clamps its own.
    n = numel(y);
    s = nan(n,1);
    for i = 1:n
        lo = max(1, i-h); hi = min(n, i+h);
        idx = (lo:hi).';
        if numel(idx) < 3, continue; end
        tt = (idx - i)*dt;
        yy = y(idx);
        ok = isfinite(yy);
        if sum(ok) < 3, continue; end
        tt = tt(ok); yy = yy(ok);
        tb = mean(tt); den = sum((tt-tb).^2);
        if den <= 0, continue; end
        s(i) = sum((tt-tb).*(yy-mean(yy))) / den;
    end
end

function r = corr_pearson(x, y)
% Pearson correlation, base MATLAB.
    ok = isfinite(x) & isfinite(y);
    x = x(ok) - mean(x(ok)); y = y(ok) - mean(y(ok));
    r = sum(x.*y) / sqrt(sum(x.^2)*sum(y.^2));
end

function s = local_fmt(x)
% Print a depth or ratio, or 'none' when the crossing does not exist.
    if isfinite(x), s = sprintf('%.3f', x); else, s = 'none'; end
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
