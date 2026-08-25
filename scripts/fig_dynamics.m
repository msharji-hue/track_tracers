% FIG_DYNAMICS (rev2) -- the four dynamic-study figures, rebuilt to follow the
% constructions in Katsuragi & Durian, Nature Physics 3, 420 (2007).
%
% NO REFITS. Every fitted parameter used here is a locked result read from a
% stored file or hard-coded with its origin. The computations this script does
% perform are display support required to reproduce KD's own constructions:
%   - per-depth intercepts at a FIXED slope (KD Fig. 3a's F(z_i)/m),
%   - the pooled common slope across depth groups (KD's own d1 procedure),
%   - ODE integration of the locked law for the t_stop curve,
%   - resample percentiles for error bars,
%   - an R^2 prediction score for the locked law.
% None of these is a new fitted result and none is presented as one.
%
% REUSE. Named in the report: the fixed-depth bin extraction of
% step5_depth_term_study, the KD ODE integrator of step5_closure_diagnostics,
% the persisted cluster resamples of step5_geometry_ladder, the ensemble-median
% + support gating of step5_fleet_fits (fig_kinematics gates on an absolute
% MinReplicates = 3 instead), the geometry identity of fig_scaling, net_accel.
% MATLAB local functions are not importable across files, so the integrator and
% bin extraction are reproduced verbatim with their source named at each site.
%
% NOTATION. c = m/d1 in public text. d1 appears only in fig_force_decomposition,
% which reproduces KD Fig. 3 and defines it in that figure's caption text. Two
% distinct c's exist and are never conflated:
%   c   = c0    trajectory-fit law parameter (fig_tstop, fig_scaling_models)
%   c_a = m/d1  acceleration-space common slope from Fig. 4A (all of Fig. 4)
%
% Base MATLAB only.

clear; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));

% ---- display switches -------------------------------------------------
SHOW_M0_REF       = false;      % ladder: dashed k0/c0 reference + CI band
D1_SOURCE         = 'panel';    % 'panel' (c_a, primary) | 'trajectory' (c0)
POOL_4B_BY_HEIGHT = false;      % 4B: pool geometries per height rank if unreadable

% ---- locked constants -------------------------------------------------
mass = 65;        % g
grav = 980;       % cm/s^2  (project constant, as get_calibration / net_accel)
k0   = 1.9163e5;  % g/s^2 = dyn/cm   trial-balanced pooled M0 (step5_fleet_fits C)
c0   = 15.898;    % g/cm             pooled M0
K0CI = [1.7769e5 2.0428e5];
C0CI = [15.025 17.022];
L0   = 0.862;     % cm   fitted d0 of d = (d0^2 H)^(1/3)   (step 3/4)
DS   = 0.751;     % cm   measured surface-release penetration, +/-0.170, n = 23
A_BARE = 2.122;   % cm^2 bare cross-section, identical for all three feet
ZC     = [0.3 0.5 0.7 0.9];   % cm, bin centres (bin-choice study)
ZHALF  = 0.1;                 % cm
V_MIN  = 15;                  % cm/s, the Step-5 frame mask
GRID_STEP = 0.025;            % cm, common depth grid
NB_CELL   = 500;              % cluster resamples (project scheme)
SEED      = 1;                % project seed

STAMP  = '20260824_221720';
EXPDIR = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_exports';
OUTDIR = 'D:\ME_GRANULAB\JerboaImpact\03_RESULTS\_figures';
MASTER = fullfile(EXPDIR, ['master_trials_' STAMP '.mat']);

% geometry identity -- fig_scaling.m lines 101-106, verified identical
MDL  = ["Tight" "Default" "Wide"];
MARK = ["o" "s" "^"];
COL  = [0.8500 0.3250 0.0980;
        0      0.4470 0.7410;
        0.9290 0.6940 0.1250];
MSZ  = 4;
KDCOL = [0.4940 0.1840 0.5560];

fprintf('\n=== fig_dynamics rev2 ===\n');

%% ===================================================================
%  SECTION 0 -- GUARDS
%  ===================================================================
fprintf('\n--- guards ---\n');

% (g1) the stamp must be the corrected export
if ~isfile(MASTER)
    error('figdyn:noMaster', 'Corrected export %s not found.', MASTER);
end
fprintf('  [g1] stamp %s asserted\n', STAMP);

% (g2) the BRANCH's exporter source must recompute via net_accel. The stamp
% assertion alone is not sufficient: a pre-#22 branch can sit beside a
% correct file and would regenerate the -2g bug on the next export.
expSrc = fullfile(fileparts(mfilename('fullpath')), 'export_master_dataset.m');
srcTxt = fileread(expSrc);
if ~contains(srcTxt, 'net_accel(kin')
    error('figdyn:staleExporter', ...
        ['%s does not call net_accel: this branch predates the a_plus_g fix ' ...
         '(PR #22) and would regenerate the -2g bug. Rebase before running.'], expSrc);
end
if contains(srcTxt, 'apg = kin.a_plus_g(:)')
    error('figdyn:staleExporter', '%s still reads kin.a_plus_g directly.', expSrc);
end
fprintf('  [g2] branch exporter calls net_accel (source checked, not just the stamp)\n');

L = load(MASTER);
T = L.T;  S = L.S;
keep = T.keep_reviewed & ~T.isZeroDrop;
K = T(keep, :);
nTr = height(K);
model = string(K.model);
tagsS = string({S.trialTag});
flagDisp = K.d_final_cm < 1 | K.v0_cm_s < 55;

% (g3) ONE source of net acceleration: the corrected export column.
fprintf('  [g3] net acceleration source = the corrected a_plus_g_cm_s2 column\n');
fprintf('       (net_accel is not re-called; the exporter already applied it)\n');

% (g4) sanity on three random trials
rng(SEED);
pick3 = randsample_local(nTr, 3);
fprintf('  [g4] sanity on 3 random trials:\n');
for q = 1:3
    i = pick3(q);
    j = find(tagsS == string(K.trialTag(i)), 1);
    ag = S(j).a_plus_g_cm_s2(:); a = S(j).a_cm_s2(:); z = S(j).z_cm(:);
    o = isfinite(ag) & isfinite(a);
    dId = max(abs(ag(o) - (grav - a(o))));
    nPre = sum(isfinite(a) & z < -0.05);
    fprintf('      %-24s |a+g-(g-a)| = %.2e | finite a at z<-0.05: %d | finite a: %d = finite a+g: %d\n', ...
        K.trialTag(i), dId, nPre, sum(isfinite(a)), sum(isfinite(ag)));
    if dId > 1e-9
        error('figdyn:apgOffset', ...
            ['%s: a+g does not equal g - a. The -1960 cm/s^2 offset is present; ' ...
             'this export is stale.'], K.trialTag(i));
    end
end
% The prescribed pre-contact (~0) and post-stop (~+g) checks CANNOT run on this
% export: kd_kinematics masks BOTH a and a+g outside [impact_index, stopFrame],
% so the retained record contains neither free fall nor rest -- the counts
% printed above show finite(a) == finite(a+g) and zero pre-contact samples.
% The check that actually excludes the -2g offset is the identity a+g = g - a,
% which the old convention (-a - g) fails by exactly 2g. It is enforced above
% per trial and fleet-wide below.
fprintf('       NOTE: pre-contact ~0 / post-stop ~+g are NOT checkable here (both a and\n');
fprintf('             a+g are NaN outside [impact, stop]); the g - a identity is the guard\n');
nExact = 0;
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    ag = S(j).a_plus_g_cm_s2(:); a = S(j).a_cm_s2(:);
    o = isfinite(ag) & isfinite(a);
    if any(o) && max(abs(ag(o) - (grav - a(o)))) < 1e-9, nExact = nExact + 1; end
end
fprintf('       fleet-wide: a+g == g - a on %d / %d trials\n', nExact, nTr);
if nExact ~= nTr, error('figdyn:apgOffset', 'a+g identity fails on %d trials.', nTr-nExact); end

% (g5) locked constants, asserted where a stored file can reproduce them
bp   = (K.v0_cm_s.^(2/3)) \ K.d_final_cm;      % pure 2/3 through the origin
l0rep = bp^(3/2) * sqrt(2*grav);                % H-form: H = v0^2/(2g)
fprintf('  [g5] constants:\n');
fprintf('       m       %.0f g            (master mass_g: %.0f)\n', mass, median(K.mass_g));
fprintf('       A_bare  %.3f cm^2       (master: %.4f)\n', A_BARE, mean(K.A_bare_cm2));
fprintf('       l0      %.3f cm         (reproduced %.4f from the pure-2/3 fit)\n', L0, l0rep);
local_assert3sf('mass',   median(K.mass_g),   mass);
local_assert3sf('A_bare', mean(K.A_bare_cm2), A_BARE);
local_assert3sf('l0',     l0rep,              L0);
fprintf('       k0, c0, d_s hard-coded (no stored file carries them; origins in the header)\n');

% (g6) the ladder resample distribution, and the six recorded ratio CIs
bootPath = fullfile(EXPDIR, 'step5_ladder_bootstrap.mat');
if ~isfile(bootPath)
    error('figdyn:noBoot', ['%s not found. Run scripts/step5_geometry_ladder.m ' ...
        '(it now persists the distribution) before this figure.'], bootPath);
end
BB = load(bootPath); BB = BB.LADDER_BOOT;
RECORDED = [1.196 0.993 1.456; 1.072 0.918 1.157; 1.115 0.940 1.399; ...
            1.108 0.998 1.281; 1.016 0.909 1.158; 1.091 0.938 1.285];
fprintf('  [g6] ladder resamples: scheme "%s", B = %d, seed %d\n', BB.scheme, BB.B, BB.seed);
r = 0; okAll = true;
for p = 1:2
    for q = 1:3
        r = r + 1;
        got = [BB.ratioPt(p,q) squeeze(BB.ratioCI(p,q,:)).'];
        if max(abs(round(got,3) - RECORDED(r,:))) > 1e-9
            fprintf('       MISMATCH %s: got %.3f [%.3f %.3f], recorded %.3f [%.3f %.3f]\n', ...
                BB.ratioName(p,q), got, RECORDED(r,:));
            okAll = false;
        end
    end
end
if ~okAll
    error('figdyn:ratioMismatch', ...
        'Persisted ladder resamples do not reproduce the recorded pairwise ratio CIs.');
end
fprintf('       all six pairwise ratio CIs reproduce to printed precision\n');

%% ===================================================================
%  SHARED DERIVED QUANTITIES
%  ===================================================================
% step-3 fits (locked; reproduced only to label the legend)
x23 = K.v0_cm_s.^(2/3); y = K.d_final_cm; N = numel(y);
c23 = [ones(N,1) x23] \ y;
nfree = fminbnd(@(n) local_rss_pow(n, K.v0_cm_s, y), 0.3, 1.2);
cfree = [ones(N,1) K.v0_cm_s.^nfree] \ y;
SStot = sum((y - mean(y)).^2);
R2_23   = 1 - sum((y - [ones(N,1) x23]*c23).^2) / SStot;
R2_free = 1 - sum((y - [ones(N,1) K.v0_cm_s.^nfree]*cfree).^2) / SStot;
% KD prediction score (a prediction, not a fit)
dKDtrial = arrayfun(@(v) local_predict_d(v, k0, mass/c0, mass, grav), K.v0_cm_s);
R2_KD = 1 - sum((y - dKDtrial).^2) / SStot;
fprintf('\n--- step-3 (locked) ---\n');
fprintf('  d0 = %.4f | a = %.5f | n = %.4f\n', c23(1), c23(2), nfree);
fprintf('  R^2: fixed-2/3 %.4f | free-n %.4f | KD prediction score %.4f\n', R2_23, R2_free, R2_KD);
% v0 -> 0 limit. 2mg/k0 is the DRAG-FREE result (mgd = k d^2/2). With c0 the
% projectile still accelerates under gravity as it sinks, so drag dissipates
% energy even from rest and the exact root is lower. Asserting the exact root
% against 2mg/k0 would therefore be wrong; what is asserted instead is that the
% implementation recovers 2mg/k0 as c -> 0, which validates the closure.
d_at0     = local_predict_d(1e-9, k0, mass/c0,      mass, grav);
d_at0_nodrag = local_predict_d(1e-9, k0, mass/(c0/1000), mass, grav);
fprintf('  KD d(v0->0) = %.4f cm with c0 | %.4f cm as c->0 | drag-free 2mg/k0 = %.4f cm\n', ...
    d_at0, d_at0_nodrag, 2*mass*grav/k0);
fprintf('      the %.1f%% gap is drag dissipated during the sink from rest, not an error\n', ...
    100*(d_at0 - 2*mass*grav/k0)/(2*mass*grav/k0));
local_assert3sf('KD d(v0=0), drag-free limit', d_at0_nodrag, 2*mass*grav/k0);

%% ===================================================================
%  FIGURE 1 -- fig_parameter_ladder (1x2, medians + resample error bars)
%  ===================================================================
fprintf('\n--- fig_parameter_ladder ---\n');
Ccell = readtable(BB.cellCsv);  Ccell.model = string(Ccell.model);
okC = isfinite(Ccell.k_fit) & isfinite(Ccell.d1_fit);
cCell = mass ./ Ccell.d1_fit;

kMed = nan(1,3); cMed = nan(1,3); kEB = nan(3,2); cEB = nan(3,2);
for g = 1:3
    s = okC & Ccell.model == MDL(g);
    kMed(g) = median(Ccell.k_fit(s));
    cMed(g) = median(cCell(s));
    kEB(g,:) = local_prctile(BB.medBoot.k(:,g), [2.5 97.5]);
    cEB(g,:) = local_prctile(BB.medBoot.c(:,g), [2.5 97.5]);
    fprintf('  %-8s k = %.4e [%.4e %.4e] | c = %6.2f [%6.2f %6.2f]\n', ...
        MDL(g), kMed(g), kEB(g,:), cMed(g), cEB(g,:));
end
fprintf('  cell-median vs pooled M0: k %+.1f%% .. %+.1f%% | c %+.1f%% .. %+.1f%%\n', ...
    100*(min(kMed)-k0)/k0, 100*(max(kMed)-k0)/k0, ...
    100*(min(cMed)-c0)/c0, 100*(max(cMed)-c0)/c0);

fig1 = figure('Color','w','Units','inches','Position',[1 1 6.75 2.85],'PaperPositionMode','auto');
subplot(1,2,1); hold on; grid on; box on;
for g = 1:3
    errorbar(g, kMed(g)/1e5, (kMed(g)-kEB(g,1))/1e5, (kEB(g,2)-kMed(g))/1e5, ...
        char(MARK(g)), 'Color', COL(g,:), 'MarkerFaceColor', COL(g,:), ...
        'MarkerSize', 7, 'LineWidth', 1.1, 'CapSize', 6);
end
if SHOW_M0_REF
    yline(k0/1e5,'--','Color',[0.5 0.5 0.5]);
    patch([0.5 3.5 3.5 0.5],[K0CI(1) K0CI(1) K0CI(2) K0CI(2)]/1e5,[0.9 0.9 0.9], ...
        'EdgeColor','none','FaceAlpha',0.4);
end
set(gca,'XTick',1:3,'XTickLabel',cellstr(MDL)); xlim([0.5 3.5]);
ylim(local_pad([min(kEB(:,1)) max(kEB(:,2))]/1e5, 0.15));
ylabel('k (10^5 g/s^2)'); title('(a)','FontWeight','normal');

subplot(1,2,2); hold on; grid on; box on;
for g = 1:3
    errorbar(g, cMed(g), cMed(g)-cEB(g,1), cEB(g,2)-cMed(g), ...
        char(MARK(g)), 'Color', COL(g,:), 'MarkerFaceColor', COL(g,:), ...
        'MarkerSize', 7, 'LineWidth', 1.1, 'CapSize', 6);
end
if SHOW_M0_REF
    yline(c0,'--','Color',[0.5 0.5 0.5]);
end
set(gca,'XTick',1:3,'XTickLabel',cellstr(MDL)); xlim([0.5 3.5]);
ylim(local_pad([min(cEB(:,1)) max(cEB(:,2))], 0.15));
ylabel('c = m/d_1 (g/cm)'); title('(b)','FontWeight','normal');
p1 = local_save(fig1, OUTDIR, 'fig_parameter_ladder');

%% ===================================================================
%  FIGURE 2 -- fig_scaling_models (legend equations + R^2)
%  ===================================================================
fprintf('\n--- fig_scaling_models ---\n');
vg = linspace(min(K.v0_cm_s), max(K.v0_cm_s), 300).';
dKDcurve = arrayfun(@(v) local_predict_d(v, k0, mass/c0, mass, grav), vg);

fig2 = figure('Color','w','Units','inches','Position',[1 1 3.375 2.85],'PaperPositionMode','auto');
hold on; grid on; box on;
hG = gobjects(1,3);
for g = 1:3
    sd = model == MDL(g) & ~flagDisp;
    hG(g) = plot(K.v0_cm_s(sd), K.d_final_cm(sd), char(MARK(g)), 'MarkerSize', MSZ, ...
        'MarkerFaceColor', COL(g,:), 'MarkerEdgeColor', COL(g,:), 'LineStyle','none');
end
h1 = plot(vg, c23(1) + c23(2)*vg.^(2/3), 'k-',  'LineWidth', 1.4);
h2 = plot(vg, cfree(1) + cfree(2)*vg.^nfree, 'k--', 'LineWidth', 1.2);
h3 = plot(vg, dKDcurve, '-', 'Color', KDCOL, 'LineWidth', 1.6);
xlabel('$v_0$ (cm/s)','Interpreter','latex');
ylabel('$d$ (cm)','Interpreter','latex');
ylim([0.5 4.5]);
% The full equation entry overflows the single-column axes, so the third entry
% uses the citation form the brief permits. Noted in the report.
lg = legend([hG h1 h2 h3], [cellstr(MDL), { ...
    sprintf('$d = d_0 + a\\,v_0^{2/3}$,  $R^2=%.3f$', R2_23), ...
    sprintf('$d = d_0 + a\\,v_0^{%.3f}$,  $R^2=%.3f$', nfree, R2_free), ...
    sprintf('KD Eq. (1), no free parameters,  $R^2=%.3f$', R2_KD)}], ...
    'Interpreter','latex','Location','northwest','FontSize',6,'Box','off');
p2 = local_save(fig2, OUTDIR, 'fig_scaling_models');

%% ===================================================================
%  FIGURE 3 -- fig_tstop (single panel, KD Fig. 1b construction)
%  ===================================================================
fprintf('\n--- fig_tstop ---\n');

% characteristic scales. L_c is a single named parameter.
D_eq = 2*sqrt(A_BARE/pi);
L_c  = D_eq;                      % the literal analogue of KD's D_b
T_c  = sqrt(L_c/grav);  V_c = sqrt(L_c*grav);
fprintf('  L_c = D_eq = %.3f cm -> T_c = %.4f s, V_c = %.1f cm/s\n', L_c, T_c, V_c);
for alt = [sqrt(A_BARE), L0, DS, mass/c0]
    fprintf('    alt L_c = %.3f cm -> T_c = %.4f s, V_c = %.1f cm/s\n', ...
        alt, sqrt(alt/grav), sqrt(alt*grav));
end

% fallback-rule fraction (a_stop = NaN marks find_stop's fallback branch)
for g = 1:3
    s = model == MDL(g);
    fprintf('  %-8s t_stop fallback fired on %.1f%% of trials\n', ...
        MDL(g), 100*mean(~isfinite(K.a_stop_cm_s2(s))));
end

% ODE model curve, reusing step5_closure_diagnostics' local_integrate_model
vgrid = 5:5:300;
tODE = nan(size(vgrid));
for q = 1:numel(vgrid)
    [tt, vv] = local_integrate_model(vgrid(q), k0, mass/c0, mass, grav, 1e-5);
    j = find(vv <= 0, 1);
    if isempty(j), j = numel(vv); end
    tODE(q) = tt(j);
end
tHarm = pi*sqrt(mass/k0);
fprintf('  ODE t_stop at v0=5 cm/s = %.4f s | drag-free harmonic limit pi*sqrt(m/k0) = %.4f s\n', ...
    tODE(1), tHarm);

fig3 = figure('Color','w','Units','inches','Position',[1 1 3.375 2.85],'PaperPositionMode','auto');
hold on; grid on; box on;
yline(T_c,'--','Color',[0.45 0.45 0.45]);
xline(V_c,'--','Color',[0.45 0.45 0.45]);
hM = gobjects(1,3);
for g = 1:3
    sa = model == MDL(g);
    hh = unique(K.dropHeight_mm(sa));
    mv = nan(numel(hh),1); mt = nan(numel(hh),1); sv = nan(numel(hh),1); st = nan(numel(hh),1);
    for q = 1:numel(hh)
        s2 = sa & K.dropHeight_mm == hh(q);
        mv(q) = mean(K.v0_cm_s(s2),'omitnan');  sv(q) = std(K.v0_cm_s(s2),'omitnan');
        mt(q) = mean(K.t_stop_s(s2),'omitnan'); st(q) = std(K.t_stop_s(s2),'omitnan');
    end
    errorbar(mv, mt, st, st, sv, sv, char(MARK(g)), 'Color', COL(g,:), ...
        'MarkerFaceColor', COL(g,:), 'MarkerSize', 4, 'LineStyle','none', ...
        'LineWidth', 0.6, 'CapSize', 0);
    hM(g) = plot(NaN, NaN, char(MARK(g)), 'Color', COL(g,:), 'MarkerFaceColor', COL(g,:), ...
        'MarkerSize', 5, 'LineStyle','none');
end
hO = plot(vgrid, tODE, 'k:', 'LineWidth', 1.5);
text(V_c+6, 0.004, sprintf('V_c = %.0f cm/s', V_c), 'FontSize', 6, 'Color', [0.35 0.35 0.35]);
text(292, T_c+0.002, sprintf('T_c = %.3f s', T_c), 'FontSize', 6, 'Color', [0.35 0.35 0.35], ...
    'HorizontalAlignment','right');
xlim([0 300]); ylim([0 max(K.t_stop_s(~flagDisp))*1.08]);
xlabel('v_0 (cm/s)'); ylabel('t_{stop} (s)');
legend([hM hO], [cellstr(MDL), {'Eq. (1), (k, c) from trajectory fit'}], ...
    'Location','northeast','FontSize',6,'Box','off');
p3 = local_save(fig3, OUTDIR, 'fig_tstop');

%% ===================================================================
%  FIGURE 4 -- fig_force_decomposition (KD Fig. 3)
%  ===================================================================
fprintf('\n--- fig_force_decomposition ---\n');

% ---- bin extraction (step5_depth_term_study's rule, unchanged) ----------
% per trial per depth: mean v^2 and mean(a+g) over |z - z_i| <= 0.1 and v > V_MIN.
% v is displayed as sqrt(mean v^2) so the parabola and the inset stay consistent.
BP = struct('g',[],'b',[],'v2',[],'ag',[],'cell',[],'v0',[]);
gv=[]; bv=[]; v2v=[]; agv=[]; cv=[]; v0v=[];
[cellId, ~, ~] = findgroups(model, K.dropHeight_mm);
for i = 1:nTr
    j = find(tagsS == string(K.trialTag(i)), 1);
    z = S(j).z_cm(:); v = S(j).v_cm_s(:); ag = S(j).a_plus_g_cm_s2(:);
    gi = find(MDL == model(i));
    for b = 1:4
        m = abs(z - ZC(b)) <= ZHALF & v > V_MIN & isfinite(ag);
        if ~any(m), continue; end
        gv(end+1,1)=gi; bv(end+1,1)=b; v2v(end+1,1)=mean(v(m).^2); %#ok<AGROW>
        agv(end+1,1)=mean(ag(m)); cv(end+1,1)=cellId(i); v0v(end+1,1)=K.v0_cm_s(i); %#ok<AGROW>
    end
end
BP.g=gv; BP.b=bv; BP.v2=v2v; BP.ag=agv; BP.cell=cv; BP.v0=v0v;
nBP = numel(agv);
fprintf('  bin points: %d trial-level (geometry x depth groups: 12)\n', nBP);

% ---- Step 1: pooled common slope across all 12 groups -------------------
grp = (BP.g-1)*4 + BP.b;                       % 1..12
Ind = double(grp == 1:12);
Xc  = [Ind BP.v2];
bc  = Xc \ BP.ag;
slope_pool = bc(13);
c_a = mass * slope_pool;
d1_panel = mass / c_a;
rssC = sum((BP.ag - Xc*bc).^2);
Xf   = [Ind Ind.*BP.v2];
bf   = Xf \ BP.ag;
rssF = sum((BP.ag - Xf*bf).^2);
Fpool = ((rssC-rssF)/(24-13)) / (rssF/(nBP-24));
pPool = local_fpval(Fpool, 11, nBP-24);
% cluster CI over cells
rng(SEED);
uCell = unique(BP.cell);
caB = nan(NB_CELL,1);
for b = 1:NB_CELL
    pk = uCell(randi(numel(uCell), numel(uCell), 1));
    id = cell2mat(arrayfun(@(u) find(BP.cell==u), pk, 'UniformOutput', false));
    Xb = [double(grp(id) == 1:12) BP.v2(id)];
    cb = Xb \ BP.ag(id);
    caB(b) = mass * cb(13);
end
caCI = local_prctile(caB, [2.5 97.5]);
fprintf('  POOLED common slope -> c_a = %.3f g/cm [%.3f %.3f] (d1 = %.4f cm)\n', ...
    c_a, caCI, d1_panel);
fprintf('  pooled parallelism (12 slopes equal?): F = %.3f, p = %.4g\n', Fpool, pPool);
fprintf('  c0 = %.3f g/cm -> c_a is %+.1f%% vs c0; c_a CI contains c0: %d\n', ...
    c0, 100*(c_a-c0)/c0, caCI(1) <= c0 && caCI(2) >= c0);

% per-geometry slopes (acceleration-space geometry check)
cGeo = nan(1,3); pGeo = nan(1,3);
for g = 1:3
    s = BP.g == g;
    Ig = double(BP.b(s) == 1:4);
    Xg = [Ig BP.v2(s)];  bg = Xg \ BP.ag(s);
    cGeo(g) = mass*bg(5);
    rC = sum((BP.ag(s) - Xg*bg).^2);
    Xg2 = [Ig Ig.*BP.v2(s)];  bg2 = Xg2 \ BP.ag(s);
    rF = sum((BP.ag(s) - Xg2*bg2).^2);
    ng = sum(s);
    pGeo(g) = local_fpval(((rC-rF)/3)/(rF/(ng-8)), 3, ng-8);
end
fprintf('  per-geometry c_bins (corrected export): %s\n', mat2str(round(cGeo,2)));
fprintf('  per-geometry parallel p            : %s\n', mat2str(round(pGeo,4)));
fprintf('  pre-correction reference           : [14.1 12.3 13.9], p = [0.073 0.339 0.000]\n');

% ---- Step 2: intercepts at FIXED slope ---------------------------------
if strcmp(D1_SOURCE,'panel'), d1_use = d1_panel; else, d1_use = mass/c0; end
fprintf('  D1_SOURCE = %s -> d1 = %.4f cm\n', D1_SOURCE, d1_use);
Fint = nan(3,4); FintCI = nan(3,4,2); FintMed = nan(3,4); nTrb = nan(3,4); nCb = nan(3,4);
rng(SEED);
for g = 1:3
    for b = 1:4
        s = BP.g == g & BP.b == b;
        res = BP.ag(s) - BP.v2(s)/d1_use;
        Fint(g,b) = mean(res);  FintMed(g,b) = median(res);
        nTrb(g,b) = sum(s);     nCb(g,b) = numel(unique(BP.cell(s)));
        uc = unique(BP.cell(s));
        fb = nan(NB_CELL,1);
        for q = 1:NB_CELL
            pk = uc(randi(numel(uc), numel(uc), 1));
            id = cell2mat(arrayfun(@(u) find(BP.cell==u & s), pk, 'UniformOutput', false));
            fb(q) = mean(BP.ag(id) - BP.v2(id)/d1_use);
        end
        FintCI(g,b,:) = local_prctile(fb, [2.5 97.5]);
    end
end
fprintf('\n  F(z_i)/m  [cluster CI]   (n trials / n cells)\n');
for g = 1:3
    for b = 1:4
        flagM = '';
        if abs(FintMed(g,b)-Fint(g,b))/abs(Fint(g,b)) > 0.10, flagM = '  <- median differs >10%'; end
        fprintf('    %-8s z=%.1f : %8.1f [%8.1f %8.1f]  (%3d / %2d)  median %8.1f%s\n', ...
            MDL(g), ZC(b), Fint(g,b), FintCI(g,b,1), FintCI(g,b,2), ...
            nTrb(g,b), nCb(g,b), FintMed(g,b), flagM);
    end
end

% outlier count (kept in the fits)
nOut = zeros(1,3);
for g = 1:3
    for b = 1:4
        s = BP.g==g & BP.b==b;
        iqrv = local_iqr(BP.ag(s));
        nOut(g) = nOut(g) + sum(abs(BP.ag(s)-median(BP.ag(s))) > 3*iqrv);
    end
end
fprintf('  trial-level outliers (>3x IQR, retained): %s\n', mat2str(nOut));

% ---- FIGURE 4 layout ---------------------------------------------------
DEPCOL = [0.78 0.86 0.94; 0.52 0.71 0.86; 0.26 0.52 0.75; 0.03 0.32 0.55];
fig4 = figure('Color','w','Units','inches','Position',[1 1 6.75 5.6],'PaperPositionMode','auto');

for g = 1:3
    ax = subplot(2,3,g); hold on; grid on; box on;
    hD = gobjects(1,4);
    for b = 1:4
        s = BP.g==g & BP.b==b;
        % one point per height cell: cell medians of the trial-level values
        uc = unique(BP.cell(s));
        cvv = nan(numel(uc),1); cag = nan(numel(uc),1);
        for q = 1:numel(uc)
            m2 = s & BP.cell==uc(q);
            cvv(q) = sqrt(median(BP.v2(m2)));  cag(q) = median(BP.ag(m2));
        end
        hD(b) = plot(cvv, cag, char(MARK(g)), 'MarkerSize', 3, ...
            'MarkerFaceColor', DEPCOL(b,:), 'MarkerEdgeColor', DEPCOL(b,:), 'LineStyle','none');
        vv = linspace(0, max(cvv)*1.05, 60);
        plot(vv, Fint(g,b) + vv.^2/d1_use, 'k:', 'LineWidth', 0.9);
        plot(0, Fint(g,b), char(MARK(g)), 'MarkerSize', 5, ...
            'MarkerFaceColor','w', 'MarkerEdgeColor', DEPCOL(b,:), 'LineWidth', 1.0);
    end
    xlabel('v (cm/s)');
    if g==1, ylabel('a+g (cm/s^2)'); end
    title(MDL(g),'FontWeight','normal');
    xlim([0 max(sqrt(BP.v2))*1.05]);
    % tighten to this panel's own data: a shared 0-2e4 range buries the
    % intercept spread (185-730 cm/s^2) that distinguishes the four depths
    sg = BP.g == g;
    ylim([min(0, min(BP.ag(sg))) max(BP.ag(sg))*1.03]);
    if g==3
        legend(hD, arrayfun(@(b) sprintf('z = %.1f cm', ZC(b)), 1:4, 'UniformOutput', false), ...
            'Location','southeast','FontSize',5,'Box','off');
    end
    % inset: same points and fits against v^2 -> parallel lines
    pos = get(ax,'Position');
    % upper-left: the data rises left-to-right, so that quadrant is free and
    % the inset's own labels do not land on the main panel's curves
    axi = axes('Position',[pos(1)+0.13*pos(3), pos(2)+0.62*pos(4), 0.33*pos(3), 0.31*pos(4)]);
    hold(axi,'on'); box(axi,'on');
    % ZOOMED to the low-v^2 end. Against v^2 the four fits are parallel lines
    % offset by F(z_i)/m, so the offset is constant in y; at full scale it is
    % a fraction of a percent of the range and the lines superimpose. This
    % window is where the separation is actually readable.
    V2MAX = 1.0e4;
    for b = 1:4
        s = BP.g==g & BP.b==b;
        uc = unique(BP.cell(s));
        cv2 = nan(numel(uc),1); cag = nan(numel(uc),1);
        for q = 1:numel(uc)
            m2 = s & BP.cell==uc(q);
            cv2(q) = median(BP.v2(m2)); cag(q) = median(BP.ag(m2));
        end
        in = cv2 <= V2MAX;
        plot(axi, cv2(in), cag(in), '.', 'MarkerSize', 4, 'Color', DEPCOL(b,:));
        vv2 = linspace(0, V2MAX, 20);
        plot(axi, vv2, Fint(g,b) + vv2/d1_use, ':', 'Color', DEPCOL(b,:), 'LineWidth', 0.9);
    end
    ylo = min(Fint(g,:)) - 400;  yhi = max(Fint(g,:)) + V2MAX/d1_use + 400;
    xlim(axi, [0 V2MAX]); ylim(axi, [ylo yhi]);
    set(axi,'FontSize',5); xlabel(axi,'v^2 (zoom)','FontSize',5); ylabel(axi,'a+g','FontSize',5);
end

% ---- 4B: collapse ------------------------------------------------------
ax2 = subplot(2,1,2); hold on; grid on; box on;
zg = (0:GRID_STEP:3.5).';
uC = unique(cellId);
curves = nan(numel(zg), numel(uC)); cV0 = nan(numel(uC),1); cGeoId = nan(numel(uC),1);
maxDepth = nan(1,3);
for q = 1:numel(uC)
    idx = find(cellId == uC(q));
    cV0(q) = median(K.v0_cm_s(idx));
    cGeoId(q) = find(MDL == model(idx(1)));
    M = nan(numel(zg), numel(idx));
    for r2 = 1:numel(idx)
        j = find(tagsS == string(K.trialTag(idx(r2))), 1);
        z = S(j).z_cm(:); v = S(j).v_cm_s(:); ag = S(j).a_plus_g_cm_s2(:);
        m = isfinite(ag) & isfinite(v);
        yv = ag - v.^2/d1_use;
        [zz, uu] = local_monotone_pair(z(m), yv(m));
        if numel(zz) < 2, continue; end
        M(:,r2) = interp1(zz, uu, zg, 'linear', NaN);
    end
    % support gate: >= 50% of the cell's trials, and never past median d_final
    sup = sum(isfinite(M), 2);
    med = median(M, 2, 'omitnan');
    med(sup < ceil(numel(idx)/2)) = NaN;
    med(zg > median(K.d_final_cm(idx))) = NaN;
    curves(:,q) = med;
    mx = max(zg(isfinite(med)));
    if ~isempty(mx), maxDepth(cGeoId(q)) = max(maxDepth(cGeoId(q)), mx); end
end
% ---- curve selection (DISPLAY ONLY -- stated on the panel and in the report)
% Two filters, both explicit and deterministic; nothing is chosen by eye.
%   (1) drop any cell whose curve STARTS below -1000 cm/s^2
%   (2) of the survivors keep the 5 with the lowest RMSE against the locked
%       k0*z/m reference, ranked -- not selected by appearance
% This shows 5 of 54 cells. It is a legibility choice, NOT evidence about how
% well the fleet collapses; the full 54-curve version is what the statistics
% in the report are computed on.
startVal = nan(numel(uC),1);
rmseRef  = nan(numel(uC),1);
for q = 1:numel(uC)
    f = find(isfinite(curves(:,q)), 1);
    if isempty(f), continue; end
    startVal(q) = curves(f,q);
    o = isfinite(curves(:,q));
    rmseRef(q) = sqrt(mean((curves(o,q) - k0*zg(o)/mass).^2));
end
elig = isfinite(startVal) & startVal >= -1000 & isfinite(rmseRef);
[~, ord] = sort(rmseRef);
ord = ord(elig(ord));
nShow = min(5, numel(ord));
sel5 = ord(1:nShow);
fprintf('  4B display: %d of %d cells eligible (start >= -1000); showing %d best-RMSE to k0 z/m\n', ...
    sum(elig), numel(uC), nShow);

cmap = parula(256);
v0lo = min(cV0); v0hi = max(cV0);
hCv = gobjects(1,nShow); lblCv = cell(1,nShow);
for q = 1:nShow
    j5 = sel5(q);
    ci = round(1 + 255*(cV0(j5)-v0lo)/(v0hi-v0lo));
    hCv(q) = plot(zg, curves(:,j5), '-', 'Color', cmap(ci,:), 'LineWidth', 1.1);
    lblCv{q} = sprintf('%s, v_0 = %.0f cm/s', MDL(cGeoId(j5)), cV0(j5));
end
% candidate F(z)/m forms, no refitting
zl = linspace(0.05, 3.5, 200);
hc1 = plot(zl, k0*zl/mass, 'k:', 'LineWidth', 1.4);
hc2 = plot(zl, grav*(1 + (3*(zl/L0).^2 - 1).*exp(-2*zl/d1_use)), 'k--', 'LineWidth', 1.2);
% intercepts from 4A, with CI and depth-window bars
OFF = [-0.015 0 0.015];
for g = 1:3
    for b = 1:4
        xx = ZC(b) + OFF(g);
        % Bars are drawn plain -- a halo on the long horizontal depth-window bar
        % punches white gaps through the curve field. Only the marker gets a
        % halo, which is what the Default blue needs against the cool colormap.
        plot([xx xx], squeeze(FintCI(g,b,:)), '-', 'Color', COL(g,:), 'LineWidth', 1.2);
        plot([xx-ZHALF xx+ZHALF], [Fint(g,b) Fint(g,b)], '-', 'Color', COL(g,:), 'LineWidth', 0.9);
        plot(xx, Fint(g,b), char(MARK(g)), 'MarkerSize', 7, 'MarkerFaceColor','w', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.4);
        plot(xx, Fint(g,b), char(MARK(g)), 'MarkerSize', 6, 'MarkerFaceColor','w', ...
            'MarkerEdgeColor', COL(g,:), 'LineWidth', 1.4);
    end
end
xlabel('z (cm)'); ylabel('a+g - v^2/d_1 (cm/s^2)');
xlim([0 3.5]);
shown = curves(:, sel5);
yl = [min(0, min(shown(:),[],'omitnan')) ...
      max([max(shown(:),[],'omitnan'), max(hc2.YData), max(FintCI(:))])];
ylim(local_pad(yl, 0.10));
legend([hc1 hc2 hCv], [{'k z / m, k from trajectory fit', ...
    sprintf('g\\{1+[3(z/d_0)^2-1]e^{-2z/d_1}\\}, d_0 = %.3f cm', L0)}, lblCv], ...
    'Location','northwest','FontSize',5.5,'Box','off');
% the selection rule is stated on the panel so the subset cannot be mistaken
% for the fleet
text(0.985, 0.04, sprintf(['%d of %d cell curves shown: start \\geq -1000 cm/s^2, ' ...
    'lowest RMSE to k z/m'], nShow, numel(uC)), 'Units','normalized', ...
    'HorizontalAlignment','right', 'FontSize', 5.5, 'Color', [0.35 0.35 0.35]);
p4 = local_save(fig4, OUTDIR, 'fig_force_decomposition');

fprintf('  4B max plotted depth per geometry: %s cm\n', mat2str(round(maxDepth,2)));

% ---- collapse-quality diagnostic (report only) -------------------------
fprintf('\n  collapse-quality Spearman(y(z), cell median v0):\n');
fprintf('    %6s %14s %14s\n', 'z', 'd1 = m/c_a', 'd1 = m/c0');
for zq = ZC
    [~, iz] = min(abs(zg - zq));
    ya = curves(iz,:).';
    % recompute the same grid row under the trajectory d1 (arithmetic shift only)
    yb = nan(numel(uC),1);
    for q = 1:numel(uC)
        idx = find(cellId == uC(q));
        vv = nan(numel(idx),1);
        for r2 = 1:numel(idx)
            j = find(tagsS == string(K.trialTag(idx(r2))), 1);
            z = S(j).z_cm(:); v = S(j).v_cm_s(:);
            m = abs(z - zq) <= GRID_STEP & isfinite(v);
            if any(m), vv(r2) = mean(v(m).^2); end
        end
        yb(q) = ya(q) + median(vv,'omitnan')*(1/d1_use - c0/mass);
    end
    fprintf('    %6.1f %14.3f %14.3f\n', zq, ...
        local_spearman(ya, cV0), local_spearman(yb, cV0));
end

%% ===================================================================
%  OUTPUT
%  ===================================================================
fprintf('\n=== WRITTEN ===\n');
for p = [p1 p2 p3 p4], fprintf('  %s\n', p); end

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function paths = local_save(fig, outDir, stem)
% Vector PDF + 300 dpi PNG, as fig_scaling. Overwrites deliberately on rerun.
    pdfPath = fullfile(outDir, [stem '.pdf']);
    pngPath = fullfile(outDir, [stem '.png']);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 300);
    paths = [string(pdfPath) string(pngPath)];
end

function [t, v, z] = local_integrate_model(v0, k, d1, mass, grav, dt)
% VERBATIM from scripts/step5_closure_diagnostics.m (local_integrate_model).
% RK4 on dv/dt = grav - k*z/mass - v^2/d1, dz/dt = v.
    f = @(zz, vv) grav - k*zz/mass - vv.^2/d1;
    nMax = 200000;
    t = nan(nMax,1); v = nan(nMax,1); z = nan(nMax,1);
    t(1) = 0; v(1) = v0; z(1) = 0;
    for i = 1:nMax-1
        vi = v(i); zi = z(i);
        k1v = f(zi, vi);                   k1z = vi;
        k2v = f(zi+dt/2*k1z, vi+dt/2*k1v); k2z = vi + dt/2*k1v;
        k3v = f(zi+dt/2*k2z, vi+dt/2*k2v); k3z = vi + dt/2*k2v;
        k4v = f(zi+dt*k3z,   vi+dt*k3v);   k4z = vi + dt*k3v;
        v(i+1) = vi + dt/6*(k1v + 2*k2v + 2*k3v + k4v);
        z(i+1) = zi + dt/6*(k1z + 2*k2z + 2*k3z + k4z);
        t(i+1) = t(i) + dt;
        if v(i+1) <= 0, break; end
    end
    ok = isfinite(t); t = t(ok); v = v(ok); z = z(ok);
end

function d_pred = local_predict_d(v0, k, d1, mass, grav)
% Root of kd_speed2_model(z) = 0 -- the Lambert-W-equivalent closure, as in
% step5_closure_diagnostics.
    f = @(z) kd_speed2_model(z, v0, k, d1, mass, grav);
    zs = linspace(1e-9, 20, 3000);
    fs = arrayfun(f, zs);
    j  = find(fs(1:end-1) > 0 & fs(2:end) <= 0, 1);
    if isempty(j), d_pred = NaN; return; end
    try, d_pred = fzero(f, [zs(j) zs(j+1)]); catch, d_pred = NaN; end
end

function [zz, uu] = local_monotone_pair(z, u)
% Strictly increasing descent segment of an arbitrary paired series, the same
% construction step5_fleet_fits uses so interp1 gets unique samples.
    zz = []; uu = [];
    ok = isfinite(z) & isfinite(u);
    z = z(ok); u = u(ok);
    if numel(z) < 2, return; end
    [~, i2] = max(z);
    z = z(1:i2); u = u(1:i2);
    keep = false(numel(z),1); last = -inf;
    for q = 1:numel(z)
        if z(q) > last, keep(q) = true; last = z(q); end
    end
    zz = z(keep); uu = u(keep);
end

function r = local_rss_pow(n, v0, y)
    X = [ones(numel(y),1) v0.^n];
    r = sum((y - X*(X\y)).^2);
end

function local_assert3sf(name, got, want)
% Agreement with a locked value at 3-significant-figure precision, judged as a
% relative difference below 5e-3 (the published constants are quoted at mixed
% precision, so a literal digit comparison is not usable).
    rel = abs(got-want)/abs(want);
    if rel > 5e-3
        error('figdyn:assert3sf', '%s: got %.6g, locked %.6g (relative %.2e > 5e-3)', ...
            name, got, want, rel);
    end
end

function p = local_fpval(F, df1, df2)
% Upper-tail F probability via the regularised incomplete beta (base MATLAB).
    if ~isfinite(F) || F <= 0, p = 1; return; end
    p = betainc(df2/(df2 + df1*F), df2/2, df1/2);
end

function idx = randsample_local(n, k)
    p = randperm(n); idx = p(1:k);
end

function lims = local_pad(v, frac)
    v = double(v); d = diff(v);
    if d == 0, d = max(abs(v(1)), 1); end
    lims = [v(1)-frac*d, v(2)+frac*d];
end

function v = local_iqr(x)
    q = local_prctile(x, [25 75]); v = q(2) - q(1);
end

function y = local_prctile(x, p)
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0, y = nan(size(p)); return; end
    if n == 1, y = repmat(x, size(p)); return; end
    q = (0.5:1:(n-0.5)) / n * 100;
    y = interp1(q, x, p, 'linear', 'extrap');
    y = min(max(y, x(1)), x(end));
end

function r = local_spearman(x, y)
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3, r = NaN; return; end
    rx = local_tiedrank(x); ry = local_tiedrank(y);
    ax = rx - mean(rx); ay = ry - mean(ry);
    r = sum(ax.*ay) / sqrt(sum(ax.^2) * sum(ay.^2));
end

function r = local_tiedrank(x)
    n = numel(x);
    [xs, i] = sort(x(:));
    ranks = (1:n).'; rs = zeros(n,1); q = 1;
    while q <= n
        j = q;
        while j < n && xs(j+1) == xs(q), j = j + 1; end
        rs(q:j) = mean(ranks(q:j)); q = j + 1;
    end
    r = zeros(n,1); r(i) = rs;
end
