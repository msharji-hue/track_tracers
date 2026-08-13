function [FIT, DATA, MEANS] = depth_scaling_literature(root, varargin)
%DEPTH_SCALING_LITERATURE  Penetration-depth scaling against published laws.
%
%   Two literature-facing tests, one figure, fixed foot (+/-7.92 deg splay).
%
%   PANEL A -- Uehara et al. (2003), reproduced as Eq.(2) of Newhall & Durian
%   (2003):
%
%       d = 0.14 * [ (rho_p / (rho_g mu^2))^(3/2) * D_p^2 * H ]^(1/3)
%
%   For a single foot rho_p and D_p are constant, so
%
%       d * sqrt(rho_g)  =  C * v0^(2/3),     C ~ 1/mu
%
%   driven by the MEASURED impact speed v0 (the published H = h + d is
%   rho_g the BULK density of the medium. Newhall use rho_g = 1.51 g/cm^3
%   for beads at a volume fraction of about 63%, i.e. rho_g = phi*rho_particle.
%   Reference slopes drawn at 2/3 (Uehara) and 0.80 (Seguin et al. 2008,
%   who report alpha = 0.40 +/- 0.04 and a field range 0.3-0.5).
%
%   PANEL B -- Newhall & Durian Eq.(4) footprint form, unnormalised:
%
%       d = K * v0^(2/3)
%
%   which is Eq.(2) after substituting H -> v0^2/2g. Kept in the original
%   variables so the PI-facing question (does d scale with v0^(2/3)?) is
%   answered directly. Reference slope drawn at 2/3.
%
%   FITTING
%   Fits are taken on PER-HEIGHT MEANS, which is what the literature plots
%   (Seguin average about ten experiments per point). Individual-trial fits
%   are also reported in the table as a robustness column; their R^2 is
%   limited by real trial-to-trial granular scatter, not by the model, so
%   the within-height relative SD is reported alongside.
%
%   CHIN is drawn on the SAME axes with open markers and dashed fits, and
%   its exponent CI is printed in the legend, because at phi ~ 0.28 and
%   63-110 um it sits outside the dry, noncohesive regime these laws were
%   established on, and its exponent is under-constrained rather than
%   demonstrably different.
%
%   MU IS A MATERIAL PROPERTY. In Eq.(2) the packing state enters only through
%   rho_g, which is already on the y axis, so the two packing states of one
%   material MUST invert to the same mu -- a powder does not change friction
%   angle because it was poured differently. mu is therefore inverted once per
%   MATERIAL (GB, CHIN), never per condition, and the within-material spread of
%   C is reported as the diagnostic: if it is not ~1, rho_g alone has failed to
%   absorb the packing change and Eq.(2) is incomplete for that material.
%
%   NOT DONE HERE: no force-law fitting (no k, no C_drag, no alpha_KD), no
%   A(z), no projected area, no mu in the plotted variable.
%
%   USAGE
%       root = 'D:\ME_GRANULAB\JerboaImpact';
%       [FIT, DATA, MEANS] = depth_scaling_literature(root);
%
%   OPTIONS (name/value)
%       'OutDir'    output folder (default <root>/03_RESULTS/_batch_logs)
%       'MinRep'    minimum trials per height to form a mean (default 3)
%       'MuRef'     mu for the reference substrate (default 0.47, GB)
%       'RefCond'   reference condition for mu inversion (default "GB/full")
%       'Select'    list of trialTags to restrict the fits to. Use this to
%                   drive them from a curated set, e.g. the cleaned Default
%                   GB/full dataset returned by load_default_gb.
%       'DriveVar'  'v0' (default) uses H_v = v0^2/2g + d, reconstructed from
%                   the MEASURED impact speed, because the carriage rail has
%                   unmeasured friction so the control height h is not the
%                   delivered energy, so v0 is used throughout. Equations
%                   are identical either way; only the value of H changes.
%
%   Base MATLAB only; no toolbox required.

% ------------------------------------------------------------------ options
if nargin < 1 || isempty(root)
    error('depth_scaling_literature:noRoot', 'Supply the JerboaImpact data root.');
end

opt.OutDir  = fullfile(root, '03_RESULTS', '_batch_logs');
opt.MinRep  = 3;
opt.MuRef   = 0.47;
opt.DriveVar= 'v0';   % 'v0' (default) or 'height'              % ring-shear internal friction, glass beads
opt.RefCond = "GB/full";
opt.Select  = [];     % list of trialTags to restrict to ([] = all found)

for i = 1:2:numel(varargin)
    f = varargin{i};
    if ~isfield(opt, f)
        error('depth_scaling_literature:badOption', 'Unknown option: %s', f);
    end
    opt.(f) = varargin{i+1};
end

FOOT_LABEL = '\pm7.92\circ foot';
M_PROJ_G   = 65;                 % total projectile mass, g (Newhall Eq.1 only)
G_CM_S2    = 980;

% ---------------------------------------------------- finalised substrate set
% phi: mean of five independent preparations, +/- one sample SD.
% Reported +/- are RANDOM (preparation/measurement) only; assumed particle
% densities are treated as exact, so absolute phi carries an additional
% systematic that cancels between conditions but not against literature phi.
% These supersede the phi values written into _kin_scalars.csv.
SUB = struct( ...
  'name',   {"GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"}, ...
  'phi',    {   0.624,      0.643,        0.280,          0.402   }, ...
  'phiSD',  {   0.004,      0.009,        0.004,          0.004   }, ...
  'rho_p',  {   2.50,       2.50,         2.35,           2.35    }, ...
  'family', {   "GB",       "GB",         "CHIN",         "CHIN"  });

fprintf('\n=== depth_scaling_literature | %s ===\n', FOOT_LABEL);
fprintf('phi source: finalised (5 preparations); CSV phi ignored.\n');

% ------------------------------------------------------- 1. gather scalars
D = dir(fullfile(root, '03_RESULTS', '**', '*_kin_scalars.csv'));
D = D(~[D.isdir]);
if isempty(D)
    error('depth_scaling_literature:noFiles', ...
          'No *_kin_scalars.csv found under %s', fullfile(root,'03_RESULTS'));
end

n = numel(D);
trialTag = strings(n,1);  condition = strings(n,1);
hTrue_mm = nan(n,1);      v0_cm_s   = nan(n,1);   d_final = nan(n,1);

for i = 1:n
    try
        T = readtable(fullfile(D(i).folder, D(i).name));
    catch ME
        warning('depth_scaling_literature:readFail','Skipping %s (%s)', ...
                D(i).name, ME.message);
        continue
    end
    if height(T) < 1, continue; end
    trialTag(i)  = string(erase(D(i).name, '_kin_scalars.csv'));
    condition(i) = local_getstr(T, 'condition');
    hTrue_mm(i)  = local_getnum(T, 'dropHeight_true_mm');
    v0_cm_s(i)   = local_getnum(T, 'v0_cm_s');
    d_final(i)   = local_getnum(T, 'd_final_cm');
end

nFound = sum(trialTag ~= "");
fprintf('kin_scalars files found : %d\n', nFound);

% Rows without a true height cannot enter an H-based fit. These are the
% stragglers whose kinematics survived the rewrite; report them by name.
noH = (trialTag ~= "") & ~isfinite(hTrue_mm);
if any(noH)
    fprintf('dropped, no dropHeight_true_mm (%d):\n', sum(noH));
    fprintf('   %s\n', trialTag(noH));
end

keep = (trialTag ~= "") & isfinite(hTrue_mm) & isfinite(v0_cm_s) & ...
       isfinite(d_final) & v0_cm_s > 0 & d_final > 0;
if ~isempty(opt.Select)
    fprintf('Select: %d of %d requested tags found\n', ...
            sum(ismember(string(opt.Select), trialTag)), numel(opt.Select));
    keep = keep & ismember(trialTag, string(opt.Select));
end
fprintf('usable trials           : %d\n', sum(keep));

% ------------------------------------------------------- 2. derived columns
DATA = table(trialTag(keep), condition(keep), hTrue_mm(keep)/10, ...
             abs(v0_cm_s(keep)), d_final(keep), ...
    'VariableNames', {'trialTag','condition','dropSetting_cm','v0_cm_s','d_final_cm'});

% finalised phi / rho_g, matched by condition
DATA.phi    = nan(height(DATA),1);
DATA.rho_g  = nan(height(DATA),1);
DATA.family = strings(height(DATA),1);
for s = 1:numel(SUB)
    m = DATA.condition == SUB(s).name;
    DATA.phi(m)    = SUB(s).phi;
    DATA.rho_g(m)  = SUB(s).phi * SUB(s).rho_p;      % bulk density of medium
    DATA.family(m) = SUB(s).family;
end
if any(~isfinite(DATA.rho_g))
    u = unique(DATA.condition(~isfinite(DATA.rho_g)));
    error('depth_scaling_literature:unknownCond', ...
          'No finalised phi for condition(s): %s', strjoin(cellstr(u), ', '));
end

% Driving variable: the MEASURED impact speed.
% Driving variable is the MEASURED impact speed. The published laws are
% written in drop distance H, but the rail friction is unmeasured so the
% release height is not the delivered energy. Substituting H = v0^2/2g moves
% every law into v0 with no change of form; a 1/3 power in distance is a 2/3
% power in speed.
fprintf('driving variable: measured v_0\n');
% Uehara Eq.(2) normalised depth, fixed foot
DATA.dNorm   = DATA.d_final_cm .* sqrt(DATA.rho_g);
% Newhall Eq.(1) energy balance: <F> d = m g H. Algebra, not a fit.
DATA.Fmean_N = (M_PROJ_G/1000) * (G_CM_S2/100) .* (DATA.v0_cm_s.^2/(2*G_CM_S2)) ./ DATA.d_final_cm;
% Ambroso et al. (2005), recast as d/d0 = (H/d0)^(1/3), i.e. d = (d0^2 H)^(1/3);
% same form used by Katsuragi & Durian (2007). Inverting, d0 = sqrt(d^3/H).
% d0 is the h=0 penetration depth -- a MEASURABLE length, not a fit parameter.
% Uses only d and H: no mu, no rho_g, no area, no mass.
DATA.d0_cm = sqrt(2*G_CM_S2) * DATA.d_final_cm.^1.5 ./ DATA.v0_cm_s;

% ------------------------------------------------------- 3. per-height means
condOrder = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
present   = unique(DATA.condition, 'stable');
conds     = [condOrder(ismember(condOrder, present)), ...
             reshape(present(~ismember(present, condOrder)), 1, [])];
nc = numel(conds);

mrows = {};
for c = 1:nc
    S  = DATA(DATA.condition == conds(c), :);
    % Group replicates by drop setting -- a LABEL for which trials are repeats.
    % It never enters a fit; everything fitted uses measured v0. Grouping by
    % rounded v0 instead splits repeats across bins and starves conditions
    % that have few settings (it dropped GB/shallow entirely).
    hr = round(S.dropSetting_cm, 2);
    hs = unique(hr);
    for k = 1:numel(hs)
        g = S(hr == hs(k), :);
        if height(g) < opt.MinRep, continue; end
        mrows{end+1} = table(conds(c), hs(k), height(g), ...
            mean(g.d_final_cm), std(g.d_final_cm), ...
            mean(g.v0_cm_s), ...
            mean(g.dNorm),  std(g.dNorm), g.rho_g(1), ...
            mean(g.d0_cm), std(g.d0_cm), ...
            mean(g.d_final_cm)/mean(g.v0_cm_s)^(2/3), ...
            std(g.d_final_cm)/mean(g.v0_cm_s)^(2/3), ...
            'VariableNames', {'condition','v0_grp','nRep', ...
                'd_mean','d_sd','v0_mean', ...
                'dNorm_mean','dNorm_sd','rho_g', ...
                'd0_mean','d0_sd','comp','comp_sd'}); %#ok<AGROW>
    end
end
if isempty(mrows)
    error('depth_scaling_literature:noMeans', ...
          'No drop height had at least MinRep = %d trials.', opt.MinRep);
end
MEANS = vertcat(mrows{:});
fprintf('\nreplicate groups per condition (>= %d trials each):\n', opt.MinRep);
for c = 1:numel(conds)
    ng = sum(MEANS.condition == conds(c));
    fprintf('  %-16s %2d groups%s\n', conds(c), ng, ...
        local_tern(ng < 3, '   <-- too few to fit', ''));
end

% ------------------------------------------------------------- 4. fits
rows = cell(nc,1);
for c = 1:nc
    S = DATA(DATA.condition == conds(c), :);
    Mn = MEANS(MEANS.condition == conds(c), :);

    if height(Mn) >= 3
        A = local_ols(log10(Mn.v0_mean), log10(Mn.dNorm_mean));   % Uehara, in v0
        B = local_ols(log10(Mn.v0_mean), log10(Mn.d_mean));       % Eq.(4)
        Cpre = mean(Mn.dNorm_mean ./ Mn.v0_mean.^(2/3));          % C at n=2/3
    else
        A = local_nanfit(); B = local_nanfit(); Cpre = NaN;
    end
    At = local_ols(log10(S.v0_cm_s), log10(S.dNorm));             % trials
    Bt = local_ols(log10(S.v0_cm_s), log10(S.d_final_cm));

    relSD = NaN;
    if ~isempty(Mn), relSD = mean(Mn.d_sd ./ Mn.d_mean, 'omitnan'); end

    rows{c} = table(conds(c), height(S), height(Mn), ...
        S.phi(1), S.rho_g(1), ...
        A.b(2), A.ci(2,1), A.ci(2,2), A.R2, ...
        (A.ci(2,1)<=2/3)&&(A.ci(2,2)>=2/3), ...
        (A.ci(2,1)<=0.80)&&(A.ci(2,2)>=0.80), ...
        At.b(2), At.ci(2,1), At.ci(2,2), At.R2, ...
        B.b(2), B.ci(2,1), B.ci(2,2), B.R2, ...
        (B.ci(2,1)<=2/3)&&(B.ci(2,2)>=2/3), ...
        Bt.b(2), Bt.ci(2,1), Bt.ci(2,2), Bt.R2, ...
        Cpre, 100*relSD, mean(S.Fmean_N), ...
        'VariableNames', {'condition','n_trials','n_heights','phi','rho_g', ...
            'a_means','a_lo95','a_hi95','R2_means','a_incl_1_3','a_incl_0_40', ...
            'a_trials','a_tr_lo95','a_tr_hi95','R2_trials', ...
            'n_v_means','n_v_lo95','n_v_hi95','R2_v_means','n_incl_2_3', ...
            'n_v_trials','n_v_tr_lo95','n_v_tr_hi95','R2_v_trials', ...
            'C_prefactor','withinH_relSD_pct','Fmean_N'});
end
FIT = vertcat(rows{:});

% ---- Ambroso length scale d0 and the COMPENSATED flatness test.
% If d = (d0^2 v0^2/2g)^(1/3) holds, then d/v0^(2/3) is a CONSTANT.
% Regressing that ratio on log10(H) should give slope 0. This is a sharper
% discriminator than the log-log exponent CI, which is wide enough to admit
% both 1/3 and 0.40.
FIT.d0_mean        = nan(height(FIT),1);
FIT.d0_sd          = nan(height(FIT),1);
FIT.d0_relspread   = nan(height(FIT),1);
FIT.comp_slope     = nan(height(FIT),1);
FIT.comp_slope_se  = nan(height(FIT),1);
FIT.comp_flat      = false(height(FIT),1);
for c = 1:height(FIT)
    Mn = MEANS(MEANS.condition == FIT.condition(c), :);
    if isempty(Mn), continue; end
    FIT.d0_mean(c) = mean(Mn.d0_mean);
    if height(Mn) > 1
        FIT.d0_sd(c)      = std(Mn.d0_mean);
        FIT.d0_relspread(c) = 100*FIT.d0_sd(c)/FIT.d0_mean(c);
    end
    if height(Mn) >= 3
        Cf = local_ols(log10(Mn.v0_mean), Mn.comp);
        FIT.comp_slope(c)    = Cf.b(2);
        FIT.comp_slope_se(c) = Cf.se(2);
        FIT.comp_flat(c)     = abs(Cf.b(2)) < 2*Cf.se(2);
    end
end

% ---- mu is a MATERIAL property, not a per-packing-state one.
% In Eq.(2) the packing state enters ONLY through rho_g, which is already on
% the y axis. So the two packing states of one material must invert to the
% SAME mu. One mu per material; the within-material spread of C is the test.
FIT.family    = strings(height(FIT),1);
for c = 1:height(FIT)
    s = find([SUB.name] == FIT.condition(c), 1);
    if ~isempty(s), FIT.family(c) = SUB(s).family; end
end

fams = unique(FIT.family, 'stable');
FIT.C_material     = nan(height(FIT),1);
FIT.C_ratio_within = nan(height(FIT),1);
for k = 1:numel(fams)
    m  = FIT.family == fams(k);
    Mm = MEANS(ismember(MEANS.condition, FIT.condition(m)), :);
    if ~isempty(Mm)
        FIT.C_material(m) = mean(Mm.dNorm_mean ./ Mm.v0_mean.^(2/3));
    end
    cc = FIT.C_prefactor(m & isfinite(FIT.C_prefactor));
    if numel(cc) > 1, FIT.C_ratio_within(m) = max(cc)/min(cc); end
end

% reference material = family of RefCond
refFam = "";
iRef = find(FIT.condition == opt.RefCond, 1);
if ~isempty(iRef), refFam = FIT.family(iRef); end
Cref = NaN;
if refFam ~= ""
    Cref = FIT.C_material(find(FIT.family == refFam, 1));
end
FIT.mu_material        = opt.MuRef * Cref ./ FIT.C_material;
FIT.theta_material_deg = atand(FIT.mu_material);

% ------------------------------------------------------------- 5. console
fprintf('\n--- Uehara Eq.(2) in v_0:  d*sqrt(rho_g) = C * v0^n   (targets 2/3, 0.80) ---\n');
for c = 1:nc
    fprintf('  %-16s means  a = %6.3f [%6.3f,%6.3f]  R2 = %.3f   %s%s\n', ...
        FIT.condition(c), FIT.a_means(c), FIT.a_lo95(c), FIT.a_hi95(c), ...
        FIT.R2_means(c), local_tern(FIT.a_incl_1_3(c),'2/3 in CI  ',''), ...
        local_tern(FIT.a_incl_0_40(c),'0.80 in CI',''));
    fprintf('  %-16s trials a = %6.3f [%6.3f,%6.3f]  R2 = %.3f   (within-height SD %.0f%% of mean)\n', ...
        '', FIT.a_trials(c), FIT.a_tr_lo95(c), FIT.a_tr_hi95(c), ...
        FIT.R2_trials(c), FIT.withinH_relSD_pct(c));
end

fprintf('\n--- Newhall Eq.(4):  d = K * v0^n   (target n = 2/3) ---\n');
for c = 1:nc
    fprintf('  %-16s means  n = %6.3f [%6.3f,%6.3f]  R2 = %.3f %s\n', ...
        FIT.condition(c), FIT.n_v_means(c), FIT.n_v_lo95(c), FIT.n_v_hi95(c), ...
        FIT.R2_v_means(c), local_tern(FIT.n_incl_2_3(c),' <- 2/3 in CI',''));
end

fprintf('\n--- Eq.(2) prefactor: per condition, then per MATERIAL ---\n');
for c = 1:nc
    fprintf('  %-16s [%-4s]  C = %6.3f\n', FIT.condition(c), FIT.family(c), FIT.C_prefactor(c));
end
fprintf('\n  Within-material consistency of C. Eq.(2) puts the packing state ONLY in\n');
fprintf('  rho_g, which is already on the y axis, so two packing states of the same\n');
fprintf('  material must give the SAME C. Ratio C_max/C_min:\n');
for k = 1:numel(fams)
    i1 = find(FIT.family == fams(k), 1);
    r  = FIT.C_ratio_within(i1);
    fprintf('    %-6s ratio = %.2f   %s\n', fams(k), r, ...
        local_tern(isfinite(r) && r < 1.15, 'consistent', ...
                   '<-- INCONSISTENT: rho_g alone does not absorb the packing change'));
end
fprintf('\n  One mu per MATERIAL (reference material %s, mu = %.2f):\n', refFam, opt.MuRef);
for k = 1:numel(fams)
    i1 = find(FIT.family == fams(k), 1);
    fprintf('    %-6s C = %6.3f   mu = %5.2f   (draining repose %4.1f deg)%s\n', ...
        fams(k), FIT.C_material(i1), FIT.mu_material(i1), FIT.theta_material_deg(i1), ...
        local_tern(FIT.theta_material_deg(i1) > 50, '   <-- above any physical repose angle', ''));
end
fprintf(['\n  NOTE: mu is INVERTED FROM THE FIT, not measured, and is reported per\n' ...
         '        MATERIAL only -- the same powder cannot have two friction angles just\n' ...
         '        because it was packed differently. Uehara define mu as tan(draining\n' ...
         '        repose angle); the %.2f reference is ring-shear internal friction, a\n' ...
         '        different quantity. Treat these as predictions to check with a funnel.\n'], opt.MuRef);

fprintf('\n--- Ambroso Eq.: d = (d0^2 H)^(1/3),  d0 = sqrt(d^3/H)  [uses d and H only] ---\n');
for c = 1:nc
    fprintf('  %-16s d0 = %5.3f +/- %5.3f cm  (spread %3.0f%% over %2d heights)\n', ...
        FIT.condition(c), FIT.d0_mean(c), FIT.d0_sd(c), ...
        FIT.d0_relspread(c), FIT.n_heights(c));
end
fprintf('\n  COMPENSATED test: if the law holds, d/v0^(2/3) is constant, so its\n');
fprintf('  regression on log10(v0) has slope 0. Sharper than the exponent CI.\n');
for c = 1:nc
    fprintf('    %-16s slope = %+6.3f +/- %5.3f   %s\n', FIT.condition(c), ...
        FIT.comp_slope(c), FIT.comp_slope_se(c), ...
        local_tern(FIT.comp_flat(c), 'FLAT  (consistent with 1/3)', 'SLOPED (law fails)'));
end
fprintf(['\n  d0 is the h = 0 penetration depth and is DIRECTLY MEASURABLE: rest the\n' ...
         '  foot on the prepared bed with the lowest toe just touching, release, and\n' ...
         '  record the depth (Ambroso do exactly this). Comparing measured d0 against\n' ...
         '  sqrt(d^3/H) tests the law with a parameter that was never fitted.\n']);

fprintf('\n--- Newhall Eq.(1) energy balance, <F> = m g H / d  (m = %d g) ---\n', M_PROJ_G);
for c = 1:nc
    fprintf('  %-16s <F> = %.2f N\n', FIT.condition(c), FIT.Fmean_N(c));
end
fprintf('\n');

% -------------------------------------------------------------- 6. figure
col = containers.Map( ...
    {'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {[0.00 0.45 0.74],[0.30 0.75 0.93],[0.85 0.33 0.10],[0.64 0.08 0.18]});
mk  = containers.Map( ...
    {'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {'o','s','^','d'});
isGB = @(s) startsWith(s, "GB");     % GB filled, CHIN open

fig = figure('Color','w','Position',[40 60 1680 520]);
tl  = tiledlayout(fig, 1, 3, 'Padding','compact', 'TileSpacing','compact');

% ---- Panel A: Uehara Eq.(2) normalised axis
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log','YScale','log');
hA = gobjects(nc,1); lblA = strings(nc,1);
for c = 1:nc
    key = char(conds(c)); cc = local_mapget(col,key,[.4 .4 .4]);
    S  = DATA(DATA.condition == conds(c), :);
    Mn = MEANS(MEANS.condition == conds(c), :);
    if isempty(Mn), continue; end
    errorbar(ax1, Mn.v0_mean, Mn.dNorm_mean, Mn.dNorm_sd, 'LineStyle','none', ...
        'Color',cc, 'LineWidth',1.1, 'CapSize',3, 'HandleVisibility','off');
    fc = local_tern(isGB(conds(c)), cc, [1 1 1]);
    hA(c) = plot(ax1, Mn.v0_mean, Mn.dNorm_mean, local_mapget(mk,key,'o'), ...
        'LineStyle','none', 'MarkerSize',7, 'MarkerEdgeColor',cc, ...
        'MarkerFaceColor',fc, 'LineWidth',1.2);
    if isfinite(FIT.a_means(c))
        xx = linspace(min(Mn.v0_mean), max(Mn.v0_mean), 50);
        b0 = mean(log10(Mn.dNorm_mean)) - FIT.a_means(c)*mean(log10(Mn.v0_mean));
        plot(ax1, xx, 10.^(FIT.a_means(c)*log10(xx) + b0), ...
            local_tern(isGB(conds(c)),'-','--'), 'Color',cc, 'LineWidth',1.7, ...
            'HandleVisibility','off');
    end
    lblA(c) = sprintf('%s   a = %.3f \\pm %.3f', conds(c), FIT.a_means(c), ...
                      (FIT.a_hi95(c)-FIT.a_lo95(c))/2);
end
% reference slopes anchored at the GB/full centroid
iG = find(conds == "GB/full", 1); if isempty(iG), iG = 1; end
MG = MEANS(MEANS.condition == conds(iG), :);
if ~isempty(MG)
    xg = exp(mean(log(MG.v0_mean)));  yg = exp(mean(log(MG.dNorm_mean)));
    xr = [min(DATA.v0_cm_s) max(DATA.v0_cm_s)];
    r1 = plot(ax1, xr, yg*(xr/xg).^(2/3),  'k--',  'LineWidth',1.3);
    r2 = plot(ax1, xr, yg*(xr/xg).^(0.80), 'k:',   'LineWidth',1.3);
else
    r1 = gobjects(1); r2 = gobjects(1);
end
xlabel(ax1,'impact velocity  v_0   (cm s^{-1})');
ylabel(ax1,'d\cdot\rho_g^{1/2}   (cm\cdot(g cm^{-3})^{1/2})');
title(ax1,'A.  Uehara Eq.(2):  d\cdot\rho_g^{1/2} = C\cdotv_0^{n}');
legend(ax1, [hA(isgraphics(hA)); r1; r2], ...
    [cellstr(lblA(isgraphics(hA))); {'slope 2/3 (Uehara)'}; {'slope 0.80 (Seguin)'}], ...
    'Location','southeast','Box','off','FontSize',8);

% ---- Panel B: Newhall Eq.(4), LINEARISED.
% Plotting d against v0^(2/3) on linear axes turns the power law into a
% straight line, so agreement is judged by straightness rather than by slope
% on a log-log plot -- where a 0.6 and a 0.75 exponent look nearly identical.
% The fitted intercept is the depth the line extrapolates to at v0 = 0.
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
hB = gobjects(nc,1); lblB = strings(nc,1);
for c = 1:nc
    key = char(conds(c)); cc = local_mapget(col,key,[.4 .4 .4]);
    S  = DATA(DATA.condition == conds(c), :);
    Mn = MEANS(MEANS.condition == conds(c), :);
    if isempty(Mn), continue; end
    xm = Mn.v0_mean.^(2/3);
    errorbar(ax2, xm, Mn.d_mean, Mn.d_sd, 'LineStyle','none', ...
        'Color',cc, 'LineWidth',1.1, 'CapSize',3, 'HandleVisibility','off');
    fc = local_tern(isGB(conds(c)), cc, [1 1 1]);
    hB(c) = plot(ax2, xm, Mn.d_mean, local_mapget(mk,key,'o'), ...
        'LineStyle','none', 'MarkerSize',7, 'MarkerEdgeColor',cc, ...
        'MarkerFaceColor',fc, 'LineWidth',1.2);
    L = local_ols(xm, Mn.d_mean);
    if isfinite(L.b(2))
        xx = linspace(min(xm), max(xm), 50);
        plot(ax2, xx, L.b(1) + L.b(2)*xx, local_tern(isGB(conds(c)),'-','--'), ...
            'Color',cc, 'LineWidth',1.7, 'HandleVisibility','off');
    end
    lblB(c) = sprintf('%s   R^2 = %.2f, intercept %.2f cm', conds(c), L.R2, L.b(1));
end
xlabel(ax2,'v_0^{2/3}   ((cm s^{-1})^{2/3})');
ylabel(ax2,'penetration depth  d\_final   (cm)');
title(ax2,'B.  Newhall Eq.(4) linearised:  d = K\cdotv_0^{2/3}');
legend(ax2, hB(isgraphics(hB)), cellstr(lblB(isgraphics(hB))), ...
    'Location','northwest','Box','off','FontSize',8);

% ---- Panel C: COMPENSATED test. d/v0^(2/3) is constant iff the law holds.
% Linear y axis on purpose: a compensated plot is read by flatness, and a log
% axis would visually flatten everything and defeat the point.
ax3 = nexttile(tl); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
set(ax3,'XScale','log');
hC = gobjects(nc,1); lblC = strings(nc,1);
for c = 1:nc
    key = char(conds(c)); cc = local_mapget(col,key,[.4 .4 .4]);
    S  = DATA(DATA.condition == conds(c), :);
    Mn = MEANS(MEANS.condition == conds(c), :);
    if isempty(Mn), continue; end
    errorbar(ax3, Mn.v0_mean, Mn.comp, Mn.comp_sd, 'LineStyle','none', ...
        'Color',cc, 'LineWidth',1.1, 'CapSize',3, 'HandleVisibility','off');
    fc = local_tern(isGB(conds(c)), cc, [1 1 1]);
    hC(c) = plot(ax3, Mn.v0_mean, Mn.comp, local_mapget(mk,key,'o'), ...
        'LineStyle','none', 'MarkerSize',7, 'MarkerEdgeColor',cc, ...
        'MarkerFaceColor',fc, 'LineWidth',1.2);
    % horizontal reference at this condition's mean -- flat if the law holds
    plot(ax3, [min(Mn.v0_mean) max(Mn.v0_mean)], mean(Mn.comp)*[1 1], ...
        local_tern(isGB(conds(c)),'-','--'), 'Color',cc, 'LineWidth',1.5, ...
        'HandleVisibility','off');
    lblC(c) = sprintf('%s   d_0 = %.2f cm, slope %+.3f\\pm%.3f', conds(c), ...
        FIT.d0_mean(c), FIT.comp_slope(c), FIT.comp_slope_se(c));
end
xlabel(ax3,'impact velocity  v_0   (cm s^{-1})');
ylabel(ax3,'d / v_0^{2/3}   ( constant if the law holds )');
title(ax3,'C.  Ambroso:  d = (d_0^2H)^{1/3}, compensated');
legend(ax3, hC(isgraphics(hC)), cellstr(lblC(isgraphics(hC))), ...
    'Location','northeast','Box','off','FontSize',8);

title(tl, sprintf(['Penetration-depth scaling vs published laws  |  %s  |  %d trials, ' ...
    '%d height means  |  small = trials, large = per-height mean \\pm 1 SD  |  ' ...
    'open markers = CHIN (outside validated regime)'], ...
    FOOT_LABEL, height(DATA), height(MEANS)), 'FontWeight','bold', 'FontSize',9);

% ---------------------------------------------------------------- 7. save
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
p = @(s) fullfile(opt.OutDir, sprintf(s, stamp));

exportgraphics(fig, p('depth_scaling_literature_%s.png'), 'Resolution', 200);
savefig(fig,        p('depth_scaling_literature_%s.fig'));
writetable(FIT,     p('depth_scaling_fits_%s.csv'));
writetable(MEANS,   p('depth_scaling_means_%s.csv'));
writetable(DATA,    p('depth_scaling_data_%s.csv'));

fprintf('wrote:\n  %s\n  %s\n  %s\n  %s\n  %s\n\n', ...
    p('depth_scaling_literature_%s.png'), p('depth_scaling_literature_%s.fig'), ...
    p('depth_scaling_fits_%s.csv'), p('depth_scaling_means_%s.csv'), ...
    p('depth_scaling_data_%s.csv'));

end % ============================================================= main end


% ------------------------------------------------------------------ helpers

function S = local_ols(x, y)
%LOCAL_OLS  y = b(1) + b(2)*x, ordinary least squares, 95% CIs.
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
n = numel(x); X = [ones(n,1), x];
if n < 3, S = local_nanfit(); S.n = n; return; end
b = X \ y; res = y - X*b; dof = n - 2;
s2 = (res'*res)/dof;
se = sqrt(diag(s2 * ((X'*X) \ eye(2))));
tc = local_tcrit95(dof);
S.b    = b;
S.se   = se;
S.ci   = [b - tc*se, b + tc*se];
S.R2   = 1 - (res'*res)/sum((y - mean(y)).^2);
S.rmse = sqrt(s2);
S.n    = n;
end

function S = local_nanfit()
S.b = [NaN;NaN]; S.se = [NaN;NaN]; S.ci = nan(2,2); S.R2 = NaN; S.rmse = NaN; S.n = 0;
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

function v = local_mapget(M, key, default)
if isKey(M, key), v = M(key); else, v = default; end
end

function s = local_tern(cond, a, b)
if cond, s = a; else, s = b; end
end
