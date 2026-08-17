function R = depth_scaling(root, varargin)
% DEPTH_SCALING  Penetration-depth scaling analysis.
%
%   'form'  'ambroso'    d = (d0^2 * H)^(1/3)         (default)
%           'powerlaw'   d ~ v0^alpha, alpha fitted
%           'literature' comparison against published forms
%           'velocity'   d vs measured v0
%   'model','condition','exclude'  passed through to load_kinematics_set
%
%   Replaces ambroso_collapse, depth_scaling_forms, depth_scaling_literature
%   and depth_velocity_scaling, which loaded the same data four times and fitted
%   variants of the same law. The fitting logic is preserved; the data source,
%   the height column and the entry point are what changed.
%
%   THE DRIVING VARIABLE IS THE MEASURED v0
%   The published laws are written in the total drop distance H = h + d, which
%   presumes a frictionless fall so that the release height alone sets the
%   energy delivered. The carriage here runs on a rail whose friction has not
%   been measured, so h is a control setting, not a measured input. Substituting
%   the free-fall relation H = v0^2/2g moves every law into v0 with no change of
%   form; a 1/3 power in distance is a 2/3 power in speed. A sqrt(2gh) column is
%   carried as v0_ff and reported as a cross-check ratio, but it is never the
%   fitted predictor.
%
%   THE HEIGHT AXIS IS dropHeight_mm
%   The drop height recorded on disk is the physical height. trialTag is an
%   identifier only and must not reach a physical axis.
%
%   FITS ARE TAKEN ON PER-HEIGHT MEANS, which is what the literature plots
%   (Seguin average about ten experiments per point). Replicates are grouped by
%   their drop setting, which is a LABEL for which trials are repeats and never
%   enters a fit. Grouping by rounded v0 instead splits repeats across bins and
%   starves conditions with few settings. Individual-trial fits are reported as
%   a robustness column; their R^2 is limited by real trial-to-trial granular
%   scatter, not by the model.
%
%   THE REPLICATE-GROUP KEY IS (model, condition, dropHeight_mm). Model is
%   part of the key because the three feet are different geometries: averaging a
%   Tight and a Default trial released from the same height would blend two
%   projectiles into one point. Fits and figures follow the same grouping, so
%   models are never pooled back together downstream. FIT and MEANS both carry a
%   model column, and the group mean v0 is the representative x with std(v0) as
%   its x-uncertainty (v0_mean / v0_sd).
%
%   The summary prints the measured v0 span per (model, condition), excluding
%   zero-drop trials. Read the fitted exponents against it: a narrow span leaves
%   a power law poorly determined however good its R^2 looks.
%
%   USAGE
%       R = depth_scaling('D:\ME_GRANULAB\JerboaImpact');
%       R = depth_scaling(root,'form','powerlaw','condition','GB/full');
%
%   OPTIONS
%       'form'       'ambroso' | 'powerlaw' | 'literature' | 'velocity'
%       'model'      '' | 'Default' | 'Tight' | 'Wide'   (default 'Default')
%       'condition'  '' (all present) | 'GB/full' | ...
%       'exclude'    trialTags to drop; honoured rather than re-derived
%       'MinRep'     min trials per drop setting to form a mean (default 3)
%       'RefV0'      reference speed for the t_stop comparison (default 250)
%       'RefCond'    reference condition for the density checks
%                    (default 'GB/full')
%       'DepthCut'   drop h > 0 trials deeper than this, cm ([] = off)
%       'OutDir'     default <root>/03_RESULTS/_batch_logs
%       'Figures'    'show' (default) | 'none'
%
%   Base MATLAB only; no toolbox required.

G = 980;                                   % cm/s^2

opt.form      = 'ambroso';
opt.model     = 'Default';
opt.condition = '';
opt.exclude   = strings(0,1);
opt.MinRep    = 3;
opt.RefV0     = 250;
opt.RefCond   = "GB/full";
opt.DepthCut  = [];
opt.OutDir    = '';
opt.Figures   = 'show';
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

form = lower(char(opt.form));
if ~ismember(form, {'ambroso','powerlaw','literature','velocity'})
    error('depth_scaling:badForm', ...
        'form must be ambroso | powerlaw | literature | velocity, got "%s".', form);
end
if isempty(opt.OutDir), opt.OutDir = fullfile(root,'03_RESULTS','_batch_logs'); end
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end

fprintf('\n=== depth_scaling (form: %s) ===\n', form);

% ── 1) data ──────────────────────────────────────────────────────────────
K = load_kinematics_set(root, 'model', opt.model, 'condition', opt.condition, ...
                              'exclude', opt.exclude, 'depthCut', opt.DepthCut);

% Fits need a real fall and a real penetration; zero-drop trials are kept by
% the loader for the measured-d0 report but cannot enter a v0 fit. Count them
% before the filter, since nothing downstream can see them afterwards.
nZeroDrop = sum(K.isZeroDrop);
K = K(~K.isZeroDrop & isfinite(K.v0_cm_s) & isfinite(K.d_final_cm) & ...
      K.v0_cm_s > 0 & K.d_final_cm > 0, :);
if height(K) < 3
    error('depth_scaling:tooFewTrials','Only %d usable trials.', height(K));
end

% ── 2) derived columns ───────────────────────────────────────────────────
% rho_g is the BULK density of the medium, phi * rho_particle. Newhall use
% 1.51 g/cm^3 for beads at phi ~ 0.63.
K.rho_g = nan(height(K),1);
cs = unique(K.condition);
for c = 1:numel(cs)
    parts = split(cs(c), "/");
    sub   = get_substrate_properties(parts(1), parts(end));
    if ~sub.ok
        error('depth_scaling:unknownCond', ...
              'No finalised phi for condition "%s": %s', cs(c), sub.reason);
    end
    K.rho_g(K.condition==cs(c)) = sub.rho_bulk_g_cm3;
end

K.h_cm      = K.dropHeight_mm / 10;
K.v0_ff     = sqrt(2*G*K.h_cm);                  % cross-check only
K.v0_ratio  = K.v0_cm_s ./ K.v0_ff;
K.H_cm      = K.v0_cm_s.^2/(2*G) + K.d_final_cm; % total travel, from v0 not h
K.d0_cm     = sqrt(2*G) * K.d_final_cm.^1.5 ./ K.v0_cm_s;   % = sqrt(d^3/H_ff)
K.dNorm     = K.d_final_cm .* sqrt(K.rho_g);     % Uehara Eq.(2), fixed foot
K.comp      = K.d_final_cm ./ K.v0_cm_s.^(2/3);  % flat iff the 2/3 law holds

fprintf('usable trials : %d   (driving variable: measured v_0)\n', height(K));
fprintf('v0/sqrt(2gh)  : median %.2f  [%.2f, %.2f]\n', ...
        median(K.v0_ratio), min(K.v0_ratio), max(K.v0_ratio));

% ── 3) per-height means ──────────────────────────────────────────────────
condOrder = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
present   = unique(K.condition,'stable');
conds     = [condOrder(ismember(condOrder,present)), ...
             reshape(present(~ismember(present,condOrder)),1,[])];

% The replicate-group key is (model, condition, dropHeight_mm). Model is in
% the key because the three feet are different geometries: averaging a Tight and
% a Default trial released from the same height would blend two different
% projectiles into one point. Before campaign 2 only Default was analysable, so
% the key could omit it without effect; from campaign 2 every trial carries a
% model suffix and it cannot.
models = unique(K.model,'stable');
mrows  = {};
for mi = 1:numel(models)
for c = 1:numel(conds)
    S  = K(K.model==models(mi) & K.condition==conds(c),:);
    if isempty(S), continue; end
    hr = round(S.h_cm,2);
    hs = unique(hr);
    for k = 1:numel(hs)
        g = S(hr==hs(k),:);
        if height(g) < opt.MinRep, continue; end
        mrows{end+1} = table(models(mi), conds(c), hs(k), height(g), ...
            mean(g.d_final_cm), std(g.d_final_cm), ...
            mean(g.H_cm), mean(g.v0_cm_s), std(g.v0_cm_s), ...
            mean(g.t_stop_s), std(g.t_stop_s), ...
            mean(g.dNorm), std(g.dNorm), ...
            mean(g.d0_cm), std(g.d0_cm), ...
            mean(g.d_final_cm)/mean(g.v0_cm_s)^(2/3), ...
            std(g.d_final_cm)/mean(g.v0_cm_s)^(2/3), ...
            g.rho_g(1), ...
            'VariableNames',{'model','condition','h_cm','nRep','d_mean','d_sd', ...
                'H_mean','v0_mean','v0_sd','tstop_mean','tstop_sd', ...
                'dNorm_mean','dNorm_sd','d0_mean','d0_sd', ...
                'comp','comp_sd','rho_g'}); %#ok<AGROW>
    end
end
end
if isempty(mrows)
    error('depth_scaling:noMeans', ...
          'No (model, condition, height) group had at least MinRep = %d trials.', ...
          opt.MinRep);
end
MEANS = vertcat(mrows{:});

fprintf('\nreplicate groups per (model, condition) (>= %d trials each):\n', opt.MinRep);
for mi = 1:numel(models)
    for c = 1:numel(conds)
        ng = sum(MEANS.model==models(mi) & MEANS.condition==conds(c));
        if ng == 0, continue; end
        fprintf('  %-8s %-16s %2d groups%s\n', models(mi), conds(c), ng, ...
            local_tern(ng < 3,'   <-- too few to fit',''));
    end
end

% Measured impact-speed span per (model, condition). Zero-drop trials are
% excluded: their v0 is 0 by protocol, so including them would report a
% meaningless floor. This is the range the fitted exponents are constrained
% over -- a narrow span makes a power-law exponent poorly determined however
% good its R^2 looks.
fprintf('\nmeasured v0 range, cm/s (zero-drop excluded):\n');
for mi = 1:numel(models)
    for c = 1:numel(conds)
        m = K.model==models(mi) & K.condition==conds(c);
        if ~any(m), continue; end
        vv = K.v0_cm_s(m);
        fprintf('  %-8s %-16s %6.1f - %6.1f   (%4.1fx span, n = %3d)\n', ...
            models(mi), conds(c), min(vv), max(vv), max(vv)/max(min(vv),eps), sum(m));
    end
end
if nZeroDrop > 0
    fprintf(['  %d zero-drop trial(s) excluded from the fits (v0 = 0 by protocol;\n' ...
             '    measured values are in v0_meas_cm_s, and load_kinematics_set\n' ...
             '    reports the measured d0 over them)\n'], nZeroDrop);
end

% ── 4) form-specific fit ─────────────────────────────────────────────────
% Fits follow the same (model, condition) grouping as the means. Pooling models
% back together at this stage would undo the point of keying the bins by model.
rows = {};
for mi = 1:numel(models)
for c = 1:numel(conds)
    M = MEANS(MEANS.model==models(mi) & MEANS.condition==conds(c),:);
    S = K(K.model==models(mi) & K.condition==conds(c),:);
    if isempty(S), continue; end
    if height(M) < 3
        fprintf('  SKIPPED %s / %s: only %d group(s) with >= %d repeats\n', ...
                models(mi), conds(c), height(M), opt.MinRep);
        continue
    end
    switch form
        case 'ambroso',    r = local_fit_ambroso(conds(c), S, M, G);
        case 'powerlaw',   r = local_fit_powerlaw(conds(c), S, M);
        case 'literature', r = local_fit_literature(conds(c), S, M, opt, G);
        case 'velocity',   r = local_fit_velocity(conds(c), S, M);
    end
    r.model = models(mi);
    rows{end+1} = movevars(r, 'model', 'Before', 1); %#ok<AGROW>
end
end
if isempty(rows)
    error('depth_scaling:noFits', ...
          'No (model, condition) group had enough replicate groups to fit.');
end
FIT = vertcat(rows{:});

% Density cross-checks that need the whole table, not one condition.
if strcmp(form,'literature')
    % Both predictions are "fixed projectile, vary the medium", so the reference
    % must be the SAME model as the row being predicted. Anchoring every model
    % to one model's reference condition would fold a geometry difference into
    % what is meant to be a pure density comparison.
    FIT.F2_alpha_pred = nan(height(FIT),1);
    FIT.t0_pred_rho14 = nan(height(FIT),1);
    for mi = 1:numel(models)
        rowsM = find(FIT.model==models(mi));
        iRef  = rowsM(FIT.condition(rowsM)==string(opt.RefCond));
        if isempty(iRef), continue; end
        iRef = iRef(1);
        % Goldman & Umbanhowar: fixed projectile -> alpha ~ rho_g^(-1/2)
        FIT.F2_alpha_pred(rowsM) = FIT.F2_alpha_s(iRef) * ...
            sqrt(FIT.rho_g(iRef)./FIT.rho_g(rowsM));
        % Their Eq.(8): t0 ~ (rho_s/rho_g)^(1/4) sqrt(R/g) -> t0 ~ rho_g^(-1/4).
        % Evaluated at the COMMON reference speed so the comparison does not
        % depend on each condition's v0 range.
        FIT.t0_pred_rho14(rowsM) = FIT.t0_at_ref_s(iRef) * ...
            (FIT.rho_g(iRef)./FIT.rho_g(rowsM)).^(1/4);
    end
    FIT.alpha_obs_over_pred = FIT.F2_alpha_s  ./ FIT.F2_alpha_pred;
    FIT.t0_obs_over_pred    = FIT.t0_at_ref_s ./ FIT.t0_pred_rho14;
end

% ── 5) console ───────────────────────────────────────────────────────────
local_report(form, FIT, opt);

% ── 6) figure ────────────────────────────────────────────────────────────
if ~strcmpi(opt.Figures,'none')
    local_figure(form, K, MEANS, FIT, conds, G);
end

% ── 7) outputs ───────────────────────────────────────────────────────────
stamp = sprintf('depth_scaling_%s', form);
writetable(FIT,   fullfile(opt.OutDir, [stamp '_fit.csv']));
writetable(MEANS, fullfile(opt.OutDir, [stamp '_means.csv']));
writetable(removevars(K, intersect({'kinPath','reason'}, K.Properties.VariableNames)), ...
           fullfile(opt.OutDir, [stamp '_trials.csv']));
fprintf('wrote %s_{fit,means,trials}.csv to %s\n\n', stamp, opt.OutDir);

R = struct('form', form, 'FIT', FIT, 'MEANS', MEANS, 'DATA', K, 'opt', opt);
end

% ═════════════════════════════════════════════════════════════════════════
% FORM 1 -- Ambroso et al. (2005) PRE 71 051305 Fig 4; Katsuragi & Durian
% (2007) Nat. Phys. 3 420 Fig 2a.  d = (d0^2 H)^(1/3), recast in v0:
%     d = (d0^2 v0^2 / 2g)^(1/3)      d0 = sqrt(2g) d^(3/2) / v0
%     d/d0 = (v0/v*)^(2/3),  v* = sqrt(2 g d0)
%     d / v0^(2/3) is CONSTANT iff the law holds
% One consequence, stated plainly: in the original form d0 is exactly the
% penetration at zero drop height, because the +d inside H survives as v0 -> 0.
% Dropping to pure v0 removes that term, so here d0 is a characteristic length
% fitted from d and v0 rather than the h = 0 depth itself. The two agree closely
% when the penetration is small next to the fall. Compare against the MEASURED
% h = 0 depth that load_kinematics_set reports.
% ═════════════════════════════════════════════════════════════════════════
function R = local_fit_ambroso(name, S, M, G)
d0  = mean(M.d0_mean);
vS  = sqrt(2*G*d0);
x   = M.v0_mean/vS;   y = M.d_mean/d0;
dev = 100*(y - x.^(2/3)) ./ x.^(2/3);
C = local_ols(log10(M.v0_mean), M.comp);            % compensated flatness
P = local_ols(log10(M.v0_mean), log10(M.d_mean));   % free exponent
R = table(string(name), height(S), height(M), d0, std(M.d0_mean), vS, ...
    min(x), max(x), sqrt(mean(dev.^2)), mean(dev), ...
    C.b(2), C.se(2), abs(C.b(2)) < 2*C.se(2), ...
    P.b(2), P.ci(2,1), P.ci(2,2), (P.ci(2,1)<=2/3)&&(P.ci(2,2)>=2/3), ...
    'VariableNames',{'condition','n_trials','n_speeds','d0_cm','d0_sd','vStar', ...
        'x_min','x_max','rms_dev_pct','mean_dev_pct', ...
        'comp_slope','comp_slope_se','comp_flat', ...
        'exponent','exp_lo','exp_hi','exp_incl_target'});
end

% ═════════════════════════════════════════════════════════════════════════
% FORM 2 -- free power law, d ~ v0^alpha, alpha fitted on per-height means.
% Newhall & Durian (2003) Eq.(4) reduced for fixed geometry predicts 2/3.
% Reported against the individual-trial fit as a robustness column.
% ═════════════════════════════════════════════════════════════════════════
function R = local_fit_powerlaw(name, S, M)
P  = local_ols(log10(M.v0_mean),  log10(M.d_mean));     % on means
Q  = local_ols(log10(S.v0_cm_s),  log10(S.d_final_cm)); % on trials
predP = 10.^(P.b(1) + P.b(2)*log10(M.v0_mean));
% within-height relative SD: the scatter the trial-level R^2 is limited by,
% which is granular, not a failure of the model
relSD = mean(M.d_sd ./ M.d_mean, 'omitnan');
R = table(string(name), height(S), height(M), ...
    P.b(2), P.ci(2,1), P.ci(2,2), (P.ci(2,1)<=2/3)&&(P.ci(2,2)>=2/3), ...
    P.b(1), 10^P.b(1), P.R2, sqrt(mean((M.d_mean-predP).^2)), ...
    Q.b(2), Q.ci(2,1), Q.ci(2,2), Q.R2, ...
    100*relSD, median(S.v0_ratio), ...
    'VariableNames',{'condition','n_trials','n_heights', ...
        'exponent','exp_lo95','exp_hi95','exp_includes_2_3', ...
        'log10K','K','R2_means','rmse_cm', ...
        'exponent_trials','exp_tr_lo95','exp_tr_hi95','R2_trials', ...
        'within_height_relSD_pct','v0ratio_med'});
end

% ═════════════════════════════════════════════════════════════════════════
% FORM 3 -- competing published forms, head to head on identical means.
%   F1  Ambroso / Katsuragi & Durian   d = (d0^2 H)^(1/3)      [1 parameter]
%       d0 = sqrt(d^3/H) is a measurable length, not a fitted constant.
%   F2  de Bruyn & Walsh; Goldman & Umbanhowar (2008)
%       d = d0 + alpha*v0                                       [2 parameters]
%       de Bruyn & Walsh model it as a Bingham fluid with a yield stress,
%       which permits a NEGATIVE intercept.
%   F3  Newhall & Durian Eq.(4) reduced, free exponent
%       d = K v0^n                                              [2 parameters]
% F1 is handicapped, one parameter against two. If it still wins, that is a
% strong result. All three are mu-free.
% Also: an F-test of the velocity forms against a CONSTANT. If velocity
% dependence is not established the fitted parameters mean nothing, and the
% flag says so rather than leaving the reader to notice. And t_stop vs v0 --
% KD report the counterintuitive signature that stopping time DECREASES with
% impact speed, so deeper penetration takes less time.
% Panel A of the old literature script (Uehara Eq.(2), d*sqrt(rho_g) ~ v0^(2/3),
% C ~ 1/mu) is carried here as the ueh_* columns.
% ═════════════════════════════════════════════════════════════════════════
function R = local_fit_literature(name, S, M, opt, G)
dm = M.d_mean; vm = M.v0_mean; Hm = M.H_mean;
nH = height(M);
SStot = sum((dm-mean(dm)).^2);

% F1 -- Ambroso, ONE parameter
d0A   = mean(sqrt(dm.^3 ./ Hm));
predA = (d0A^2 * Hm).^(1/3);
rmseA = sqrt(mean((dm-predA).^2));
R2A   = 1 - sum((dm-predA).^2)/SStot;

% F2 -- linear in v0, TWO parameters
L = local_ols(vm, dm);
SSE2 = sum((dm-(L.b(1)+L.b(2)*vm)).^2);
F2stat = ((SStot - SSE2)/1) / (SSE2/(nH-2));

% F3 -- free power law in v0, TWO parameters
P     = local_ols(log10(vm), log10(dm));
predP = 10.^(P.b(1) + P.b(2)*log10(vm));
rmseP = sqrt(mean((dm-predP).^2));

% Uehara Eq.(2) normalised depth, fixed foot: d*sqrt(rho_g) ~ v0^(2/3), C ~ 1/mu
U = local_ols(log10(vm), log10(M.dNorm_mean));

% stopping time trend, in PHYSICAL units. No normalisation and no plateau fit:
% KD themselves plot t_stop [s] vs v0 [cm/s] and merely DRAW the two projectile
% size scales as reference lines.
Tst     = local_ols(vm, M.tstop_mean);
t0ref   = Tst.b(1) + Tst.b(2)*opt.RefV0;
t0refSE = sqrt(Tst.s2*(1/Tst.n + (opt.RefV0-Tst.xbar)^2/Tst.Sxx));
extrap  = opt.RefV0 < min(vm) || opt.RefV0 > max(vm);

vStar = sqrt(d0A*G);   tStar = sqrt(d0A/G);

R = table(string(name), height(S), nH, M.rho_g(1), median(S.v0_ratio), ...
    d0A, rmseA, R2A, ...
    L.b(1), L.ci(1,1), L.ci(1,2), L.b(2), L.ci(2,1), L.ci(2,2), L.R2, L.rmse, ...
    F2stat, local_fcrit(nH-2), F2stat > local_fcrit(nH-2), ...
    P.b(2), P.ci(2,1), P.ci(2,2), rmseP, ...
    U.b(2), U.ci(2,1), U.ci(2,2), 10^U.b(1), ...
    Tst.b(2), Tst.ci(2,1), Tst.ci(2,2), ...
    t0ref, t0ref-Tst.tc*t0refSE, t0ref+Tst.tc*t0refSE, extrap, ...
    vStar, tStar, min(vm)/vStar, ...
    local_best(rmseA, L.rmse, rmseP), ...
    'VariableNames',{'condition','n_trials','n_heights','rho_g','v0_ratio_median', ...
        'F1_d0_cm','F1_rmse','F1_R2', ...
        'F2_d0_cm','F2_d0_lo','F2_d0_hi','F2_alpha_s','F2_a_lo','F2_a_hi', ...
        'F2_R2','F2_rmse','F_vs_const','F_crit05','vel_dep_established', ...
        'F3_n','F3_n_lo','F3_n_hi','F3_rmse', ...
        'ueh_n','ueh_n_lo','ueh_n_hi','ueh_C', ...
        'dtstop_dv0','dt_lo','dt_hi', ...
        't0_at_ref_s','t0_ref_lo','t0_ref_hi','t0_extrapolated', ...
        'vStar_cm_s','tStar_s','v0min_over_vStar','best_form'});
end

% ═════════════════════════════════════════════════════════════════════════
% FORM 4 -- L0 baseline: depth vs MEASURED impact velocity, per trial.
% Tests only the fixed-geometry reduction of Newhall & Durian Eq.(4),
%     d = K v0^(2/3)
% No force-law fitting, no k or C, no A(z), no area normalisation, no pooling
% across substrates. Fit A frees the exponent; fit B fixes it at 2/3 and reads
% the intercept as a d0.
% ═════════════════════════════════════════════════════════════════════════
function R = local_fit_velocity(name, S, M)
v = S.v0_cm_s;  d = S.d_final_cm;  r = S.v0_ratio;
A = local_ols(log10(v), log10(d));    % free exponent
B = local_ols(v.^(2/3), d);           % exponent fixed at 2/3
R = table(string(name), numel(v), height(M), ...
    min(v), max(v), min(d), max(d), ...
    A.b(2), A.ci(2,1), A.ci(2,2), (A.ci(2,1)<=2/3)&&(A.ci(2,2)>=2/3), ...
    A.b(1), A.R2, ...
    B.b(2), B.b(1), B.ci(1,1), B.ci(1,2), B.R2, B.rmse, ...
    median(r), min(r), max(r), sum(r < 0.7 | r > 1.3), ...
    'VariableNames',{'condition','n','n_heights', ...
        'v0_min','v0_max','d_min','d_max', ...
        'exponent','exp_lo95','exp_hi95','exp_includes_2_3', ...
        'log10K','R2_loglog', ...
        'slope_a','intercept_d0','d0_lo95','d0_hi95','R2_linear','rmse_cm', ...
        'v0ratio_med','v0ratio_min','v0ratio_max','n_v0ratio_off'});
end

% ═════════════════════════════════════════════════════════════════════════
function local_report(form, FIT, opt)
nc  = height(FIT);
% Rows are keyed by (model, condition); label them so the two are never
% conflated in the printed summary.
lbl = FIT.model + " / " + FIT.condition;
switch form
    case 'ambroso'
        fprintf('\n--- d0, the single length (from d and v_0 only) ---\n');
        for c = 1:nc
            fprintf('  %-26s d0 = %5.3f +/- %5.3f cm   v* = %5.1f cm/s   (%2d speeds, %3d trials)\n', ...
                lbl(c), FIT.d0_cm(c), FIT.d0_sd(c), FIT.vStar(c), ...
                FIT.n_speeds(c), FIT.n_trials(c));
        end
        fprintf('\n--- Collapse quality: deviation from d/d0 = (v_0/v*)^(2/3) ---\n');
        for c = 1:nc
            fprintf('  %-26s v_0/v* spans %5.2f - %5.2f (%4.1fx)   rms deviation %5.1f%%\n', ...
                lbl(c), FIT.x_min(c), FIT.x_max(c), ...
                FIT.x_max(c)/FIT.x_min(c), FIT.rms_dev_pct(c));
        end
        fprintf('\n--- Compensated test: d/v_0^{2/3} is constant iff the law holds ---\n');
        for c = 1:nc
            fprintf('  %-26s slope on log10(v_0) = %+7.4f +/- %6.4f   %s\n', ...
                lbl(c), FIT.comp_slope(c), FIT.comp_slope_se(c), ...
                local_tern(FIT.comp_flat(c),'FLAT','SLOPED -- law does not hold here'));
        end
        fprintf('\n--- Free exponent, cross-check (target 2/3 = 0.667) ---\n');
        for c = 1:nc
            fprintf('  %-26s n = %5.3f [%5.3f, %5.3f]%s\n', lbl(c), ...
                FIT.exponent(c), FIT.exp_lo(c), FIT.exp_hi(c), ...
                local_tern(FIT.exp_incl_target(c),'   2/3 in CI',''));
        end

    case 'powerlaw'
        fprintf('\n--- d ~ v0^n on per-height means (target n = 2/3 = 0.667) ---\n');
        for c = 1:nc
            fprintf('  %-26s n = %6.3f [%6.3f, %6.3f]  K = %6.4f  R2 = %.3f  rmse = %.3f cm%s\n', ...
                lbl(c), FIT.exponent(c), FIT.exp_lo95(c), FIT.exp_hi95(c), ...
                FIT.K(c), FIT.R2_means(c), FIT.rmse_cm(c), ...
                local_tern(FIT.exp_includes_2_3(c),'  <- 2/3 in CI',''));
        end
        fprintf('\n--- robustness: same fit on individual trials ---\n');
        for c = 1:nc
            fprintf('  %-26s n = %6.3f [%6.3f, %6.3f]  R2 = %.3f   within-height relative SD %.1f%%\n', ...
                lbl(c), FIT.exponent_trials(c), FIT.exp_tr_lo95(c), ...
                FIT.exp_tr_hi95(c), FIT.R2_trials(c), FIT.within_height_relSD_pct(c));
        end
        fprintf(['  The trial-level R2 is limited by real granular scatter (the column\n' ...
                 '  above), not by the model.\n']);

    case 'literature'
        fprintf('\n--- Head to head on identical per-height means (lower rmse wins) ---\n');
        fprintf('  %-26s %-22s %-22s %-22s %s\n','condition', ...
                'F1 Ambroso (1 par)','F2 linear v0 (2 par)','F3 power v0 (2 par)','best');
        for c = 1:nc
            fprintf('  %-26s rmse %6.4f R2 %5.3f   rmse %6.4f R2 %5.3f   rmse %6.4f n %5.3f   %s\n', ...
                lbl(c), FIT.F1_rmse(c), FIT.F1_R2(c), ...
                FIT.F2_rmse(c), FIT.F2_R2(c), FIT.F3_rmse(c), FIT.F3_n(c), ...
                FIT.best_form(c));
        end
        fprintf('\n--- Is velocity dependence even established? F(1,n-2) vs a constant ---\n');
        for c = 1:nc
            fprintf('  %-26s F = %8.2f  (crit %.2f)  %s\n', lbl(c), ...
                FIT.F_vs_const(c), FIT.F_crit05(c), ...
                local_tern(FIT.vel_dep_established(c),'established', ...
                           'NOT established -- fitted parameters are meaningless here'));
        end
        fprintf('\n--- Uehara Eq.(2): d*sqrt(rho_g) ~ v0^n, C ~ 1/mu (target n = 2/3) ---\n');
        for c = 1:nc
            fprintf('  %-26s n = %6.3f [%6.3f, %6.3f]   C = %7.4f   rho_g = %.3f g/cm^3\n', ...
                lbl(c), FIT.ueh_n(c), FIT.ueh_n_lo(c), FIT.ueh_n_hi(c), ...
                FIT.ueh_C(c), FIT.rho_g(c));
        end
        fprintf('\n--- alpha ~ rho_g^(-1/2) (Goldman & Umbanhowar, fixed projectile) ---\n');
        for c = 1:nc
            fprintf('  %-26s alpha = %7.5f s   predicted %7.5f   obs/pred = %5.2f\n', ...
                lbl(c), FIT.F2_alpha_s(c), FIT.F2_alpha_pred(c), ...
                FIT.alpha_obs_over_pred(c));
        end
        fprintf('\n--- t_stop vs v0 at the common reference speed %g cm/s ---\n', opt.RefV0);
        for c = 1:nc
            fprintf('  %-26s dt/dv0 = %+9.3e s/(cm/s) [%+9.3e, %+9.3e]   t0 = %.4f s%s\n', ...
                lbl(c), FIT.dtstop_dv0(c), FIT.dt_lo(c), FIT.dt_hi(c), ...
                FIT.t0_at_ref_s(c), local_tern(FIT.t0_extrapolated(c),'  (EXTRAPOLATED)',''));
        end
        fprintf(['  KD report that stopping time DECREASES with impact speed, so a\n' ...
                 '  negative dt/dv0 reproduces their signature.\n']);

    case 'velocity'
        fprintf('\n--- counts ---\n');
        for c = 1:nc
            fprintf('  %-26s n = %3d   over %d heights\n', ...
                lbl(c), FIT.n(c), FIT.n_heights(c));
        end
        fprintf('\n--- Fit A: free exponent, d ~ v0^n   (target n = 2/3 = 0.667) ---\n');
        for c = 1:nc
            fprintf('  %-26s n = %6.3f  [%6.3f, %6.3f]  R2 = %.3f %s\n', ...
                lbl(c), FIT.exponent(c), FIT.exp_lo95(c), FIT.exp_hi95(c), ...
                FIT.R2_loglog(c), local_tern(FIT.exp_includes_2_3(c),' <- 2/3 in CI',''));
        end
        fprintf('\n--- Fit B: d = a*v0^(2/3) + d0 ---\n');
        for c = 1:nc
            fprintf('  %-26s a = %7.4f   d0 = %6.3f cm [%6.3f, %6.3f]   R2 = %.3f   rmse = %.3f cm\n', ...
                lbl(c), FIT.slope_a(c), FIT.intercept_d0(c), ...
                FIT.d0_lo95(c), FIT.d0_hi95(c), FIT.R2_linear(c), FIT.rmse_cm(c));
        end
        fprintf('\n--- v0 reliability, measured / sqrt(2 g h_true) ---\n');
        for c = 1:nc
            fprintf('  %-26s median %.2f  range %.2f - %.2f   %d trial(s) outside 0.7-1.3\n', ...
                lbl(c), FIT.v0ratio_med(c), FIT.v0ratio_min(c), ...
                FIT.v0ratio_max(c), FIT.n_v0ratio_off(c));
        end
end
fprintf('\n');
end

% ═════════════════════════════════════════════════════════════════════════
function local_figure(form, K, MEANS, FIT, conds, G) %#ok<INUSD>
% Two panels: physical units on the left, the form's collapse on the right.
% Colour encodes CONDITION, marker shape encodes MODEL, so the three feet stay
% visually distinct now that the bins are keyed by (model, condition, height).
% CHIN is drawn with open markers because at phi ~ 0.28 and 63-110 um it sits
% outside the dry, noncohesive regime these laws were established on; its
% exponent is under-constrained rather than demonstrably different.
col = containers.Map({'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {[0 .45 .74],[.30 .75 .93],[.85 .33 .10],[.64 .08 .18]});
mkModel = containers.Map({'Default','Tight','Wide'}, {'o','^','s'});
isGB = @(s) startsWith(s,"GB");

% Every (model, condition) pair actually present, in MEANS order.
[grp, gModel, gCond] = findgroups(MEANS.model, MEANS.condition);
nG = max(grp);

fig = figure('Color','w','Position',[60 60 1180 500]); %#ok<NASGU>
tl  = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
title(tl, sprintf('depth scaling -- form: %s', form), 'FontWeight','bold');

% ── Panel A: d vs v0, physical units, log-log ────────────────────────────
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log','YScale','log');
hA = gobjects(nG,1); lA = strings(nG,1);
for g = 1:nG
    M = MEANS(grp==g,:);
    if isempty(M), continue; end
    [cc, mkr, fc] = local_style(gCond(g), gModel(g), col, mkModel, isGB);
    errorbar(ax1, M.v0_mean, M.d_mean, M.d_sd, 'LineStyle','none','Color',cc, ...
        'LineWidth',1.1,'CapSize',3,'HandleVisibility','off');
    hA(g) = plot(ax1, M.v0_mean, M.d_mean, mkr, 'LineStyle','none', ...
        'MarkerSize',7, 'MarkerEdgeColor',cc, 'MarkerFaceColor',fc, 'LineWidth',1.2);
    lA(g) = gModel(g) + " / " + gCond(g);
end
% reference slope 2/3, anchored on the data
vv = MEANS.v0_mean;  dd = MEANS.d_mean;
xr = [min(vv) max(vv)];
anchor = median(dd ./ vv.^(2/3));
plot(ax1, xr, anchor*xr.^(2/3), 'k--', 'LineWidth',1.2);
xlabel(ax1,'measured v_0  (cm/s)'); ylabel(ax1,'penetration depth d  (cm)');
title(ax1,'d vs measured v_0  (dashed: slope 2/3)');
ok = isgraphics(hA);
legend(ax1, [hA(ok); plot(ax1,NaN,NaN,'k--')], [lA(ok); "v_0^{2/3}"], ...
       'Location','northwest','Box','off');

% ── Panel B: the form's collapse ─────────────────────────────────────────
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
switch form
    case {'ambroso','literature'}
        % d/d0 vs v0/v*, the collapse the law predicts. d0 is read from the fit
        % row for the SAME (model, condition), not the first matching condition.
        set(ax2,'XScale','log','YScale','log');
        for g = 1:nG
            i = find(FIT.model==gModel(g) & FIT.condition==gCond(g), 1);
            if isempty(i), continue; end
            if strcmp(form,'ambroso'), d0 = FIT.d0_cm(i); else, d0 = FIT.F1_d0_cm(i); end
            vS = sqrt(2*G*d0);
            M  = MEANS(grp==g,:);
            [cc, mkr, fc] = local_style(gCond(g), gModel(g), col, mkModel, isGB);
            plot(ax2, M.v0_mean/vS, M.d_mean/d0, mkr, 'LineStyle','none', ...
                'MarkerSize',7,'MarkerEdgeColor',cc,'MarkerFaceColor',fc, ...
                'LineWidth',1.2, 'DisplayName',char(gModel(g)+" / "+gCond(g)));
        end
        xl = xlim(ax2);
        plot(ax2, xl, xl.^(2/3), 'k--','LineWidth',1.2,'DisplayName','(v_0/v^*)^{2/3}');
        xlabel(ax2,'v_0 / v^*'); ylabel(ax2,'d / d_0');
        title(ax2,'Ambroso collapse');
    otherwise
        % compensated plot: d/v0^(2/3) is FLAT iff the law holds
        set(ax2,'XScale','log');
        for g = 1:nG
            M = MEANS(grp==g,:);
            if isempty(M), continue; end
            [cc, mkr, fc] = local_style(gCond(g), gModel(g), col, mkModel, isGB);
            errorbar(ax2, M.v0_mean, M.comp, M.comp_sd, 'LineStyle','none', ...
                'Color',cc,'LineWidth',1.1,'CapSize',3,'HandleVisibility','off');
            plot(ax2, M.v0_mean, M.comp, mkr, 'LineStyle','none','MarkerSize',7, ...
                'MarkerEdgeColor',cc,'MarkerFaceColor',fc,'LineWidth',1.2, ...
                'DisplayName',char(gModel(g)+" / "+gCond(g)));
        end
        xlabel(ax2,'measured v_0  (cm/s)');
        ylabel(ax2,'d / v_0^{2/3}   (cm / (cm/s)^{2/3})');
        title(ax2,'Compensated: flat iff d \propto v_0^{2/3}');
end
legend(ax2,'Location','best','Box','off');
end

function [cc, mkr, fc] = local_style(condName, modelName, col, mkModel, isGB)
k = char(condName);
if isKey(col,k), cc = col(k); else, cc = [0.3 0.3 0.3]; end
km = char(modelName);
if isKey(mkModel,km), mkr = mkModel(km); else, mkr = 'd'; end
if isGB(condName), fc = cc; else, fc = [1 1 1]; end
end

% ═════════════════════════════════════════════════════════════════════════
function S = local_ols(x,y)
x=x(:); y=y(:); ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
n=numel(x); X=[ones(n,1),x];
if n<3
    S.b=[NaN;NaN]; S.se=[NaN;NaN]; S.ci=nan(2,2); S.R2=NaN; S.rmse=NaN;
    S.n=n; S.xbar=NaN; S.Sxx=NaN; S.s2=NaN; S.tc=NaN; return
end
b=X\y; r=y-X*b; s2=(r'*r)/(n-2); se=sqrt(diag(s2*((X'*X)\eye(2))));
tc=local_t95(n-2);
S.b=b; S.se=se; S.ci=[b-tc*se, b+tc*se];
S.R2=1-(r'*r)/sum((y-mean(y)).^2); S.rmse=sqrt(s2);
S.n=n; S.xbar=mean(x); S.Sxx=sum((x-mean(x)).^2); S.s2=s2; S.tc=tc;
end

function t = local_t95(dof), try, t=tinv(0.975,dof); catch, t=1.96; end, end

function f = local_fcrit(dof2)
% F(1, dof2) upper 5% critical value
tbl=[1 161;2 18.5;3 10.13;4 7.71;5 6.61;6 5.99;7 5.59;8 5.32;9 5.12;10 4.96; ...
     12 4.75;14 4.60;16 4.49;18 4.41;20 4.35;25 4.24;30 4.17];
if dof2>=30, f=4.0; else, f=interp1(tbl(:,1),tbl(:,2),max(dof2,1),'linear','extrap'); end
end

function s = local_best(a,b,c)
[~,i]=min([a b c]); n=["F1 Ambroso","F2 linear","F3 power"]; s=n(i);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
