function [FIT, MEANS] = depth_scaling_forms(root, varargin)
%DEPTH_SCALING_FORMS  Competing published forms, head to head. READ-ONLY.
%
%   Companion to depth_scaling_literature.m. That script tests the Uehara /
%   Newhall / Ambroso branch. This one asks a different question: given the
%   SAME per-height means, which published functional form actually describes
%   each condition best?
%
%   THE CANDIDATES
%     F1  Ambroso et al. (2005); Katsuragi & Durian (2007)
%           d = (d0^2 H)^(1/3),  H = h_true + d      [1 parameter]
%         d0 = sqrt(d^3/H) is a measurable length, not a fitted constant.
%
%     F2  de Bruyn & Walsh; Goldman & Umbanhowar (2008); cited in KD (2007)
%           d = d0 + alpha*v0                        [2 parameters]
%         Goldman & Umbanhowar find this for spheres penetrating more than
%         about a radius. de Bruyn & Walsh model it as a Bingham fluid with a
%         yield stress, which permits a NEGATIVE intercept.
%
%     F3  Newhall & Durian Eq.(4) reduced for fixed geometry
%           d = K*v0^n, free n                       [2 parameters]
%
%   F1 is handicapped -- one parameter against two. If it still wins, that is
%   a strong result. All three are mu-free.
%
%   ALSO TESTED
%     - F-test of each velocity form against a CONSTANT. If velocity
%       dependence is not established, the fitted parameters mean nothing and
%       the flag says so rather than leaving the reader to notice.
%     - alpha ~ rho_g^(-1/2), predicted by Goldman & Umbanhowar for fixed
%       projectile. Uses finalised phi and rho_particle. No mu required.
%     - t_stop vs v0. KD report the counterintuitive signature that stopping
%       time DECREASES with impact speed, so deeper penetration takes less
%       time. Parameter-free, and t_stop_s is already in the scalars CSV.
%
%   v0 RELIABILITY is reported per condition as median v0/v_freefall. F2 and
%   F3 both use v0 as the predictor, so read their fits against that number.
%
%   USAGE
%       [FIT, MEANS] = depth_scaling_forms('D:\ME_GRANULAB\JerboaImpact');
%
%   OPTIONS
%       'OutDir'   default <root>/03_RESULTS/_batch_logs
%       'MinRep'   min trials per height to form a mean (default 3)
%       'RefCond'  reference for the alpha density check (default "GB/full")

opt.OutDir  = fullfile(root,'03_RESULTS','_batch_logs');
opt.MinRep  = 3;
opt.RefCond = "GB/full";
opt.RefV0   = 200;      % cm/s, common reference for the t0 comparison
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

G = 980;
FOOT = '\pm7.92\circ foot';

% finalised phi (5 preparations); supersedes the phi written into the CSVs
SUB = struct( ...
  'name',  {"GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"}, ...
  'phi',   {   0.624,      0.643,        0.280,          0.402   }, ...
  'rho_p', {   2.50,       2.50,         2.35,           2.35    });

fprintf('\n=== depth_scaling_forms | %s ===\n', FOOT);

% ------------------------------------------------------------- 1. load
D = dir(fullfile(root,'03_RESULTS','**','*_kin_scalars.csv')); D = D(~[D.isdir]);
if isempty(D), error('depth_scaling_forms:noFiles','No *_kin_scalars.csv found.'); end

n = numel(D);
tag=strings(n,1); cond=strings(n,1);
h=nan(n,1); v0=nan(n,1); d=nan(n,1); ts=nan(n,1);
for i = 1:n
    T = readtable(fullfile(D(i).folder,D(i).name));
    if height(T)<1, continue; end
    tag(i)  = string(erase(D(i).name,'_kin_scalars.csv'));
    cond(i) = local_str(T,'condition');
    h(i)    = local_num(T,'dropHeight_true_mm');
    v0(i)   = abs(local_num(T,'v0_cm_s'));
    d(i)    = local_num(T,'d_final_cm');
    ts(i)   = local_num(T,'t_stop_s');
end
keep = tag~="" & isfinite(h) & isfinite(v0) & isfinite(d) & v0>0 & d>0;
DATA = table(tag(keep),cond(keep),h(keep)/10,v0(keep),d(keep),ts(keep), ...
    'VariableNames',{'trialTag','condition','h_cm','v0_cm_s','d_cm','t_stop_s'});
DATA.H_cm = DATA.h_cm + DATA.d_cm;
DATA.v0_ff = sqrt(2*G*DATA.h_cm);
DATA.v0_ratio = DATA.v0_cm_s ./ DATA.v0_ff;
DATA.rho_g = nan(height(DATA),1);
for s = 1:numel(SUB)
    DATA.rho_g(DATA.condition==SUB(s).name) = SUB(s).phi*SUB(s).rho_p;
end
fprintf('trials: %d\n', height(DATA));

% -------------------------------------------------- 2. per-height means
condOrder = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
present = unique(DATA.condition,'stable');
conds = [condOrder(ismember(condOrder,present)), ...
         reshape(present(~ismember(present,condOrder)),1,[])];
mrows = {};
for c = 1:numel(conds)
    S = DATA(DATA.condition==conds(c),:);
    hr = round(S.h_cm,2); hs = unique(hr);
    for k = 1:numel(hs)
        g = S(hr==hs(k),:);
        if height(g) < opt.MinRep, continue; end
        mrows{end+1} = table(conds(c), hs(k), height(g), ...
            mean(g.d_cm), std(g.d_cm), mean(g.H_cm), mean(g.v0_cm_s), ...
            mean(g.t_stop_s), std(g.t_stop_s), g.rho_g(1), ...
            'VariableNames',{'condition','h_cm','nRep','d_mean','d_sd', ...
                'H_mean','v0_mean','tstop_mean','tstop_sd','rho_g'}); %#ok<AGROW>
    end
end
MEANS = vertcat(mrows{:});

% ------------------------------------------------------------- 3. fits
rows = {};
for c = 1:numel(conds)
    M = MEANS(MEANS.condition==conds(c),:);
    S = DATA(DATA.condition==conds(c),:);
    nH = height(M);
    if nH < 4
        fprintf('  %s: only %d height groups, skipped\n', conds(c), nH); continue
    end
    dm = M.d_mean; vm = M.v0_mean; Hm = M.H_mean;
    SStot = sum((dm-mean(dm)).^2);

    % F1 -- Ambroso, ONE parameter
    d0A   = mean(sqrt(dm.^3 ./ Hm));
    predA = (d0A^2 * Hm).^(1/3);
    rmseA = sqrt(mean((dm-predA).^2));
    R2A   = 1 - sum((dm-predA).^2)/SStot;

    % F2 -- linear in v0, TWO parameters
    L = local_ols(vm, dm);
    F2stat = ((SStot - sum((dm-(L.b(1)+L.b(2)*vm)).^2))/1) / ...
             (sum((dm-(L.b(1)+L.b(2)*vm)).^2)/(nH-2));

    % F3 -- free power law in v0, TWO parameters
    P = local_ols(log10(vm), log10(dm));
    predP = 10.^(P.b(1) + P.b(2)*log10(vm));
    rmseP = sqrt(mean((dm-predP).^2));

    % stopping time trend, in PHYSICAL units. No normalisation, no plateau
    % fit: KD themselves plot t_stop [s] vs v0 [cm/s] and merely DRAW the two
    % projectile-size scales as reference lines.
    Tst = local_ols(vm, M.tstop_mean);

    % t_stop evaluated at a COMMON reference speed, so conditions with
    % different v0 ranges are compared at the same place on the curve
    t0ref   = Tst.b(1) + Tst.b(2)*opt.RefV0;
    t0refSE = sqrt(Tst.s2*(1/Tst.n + (opt.RefV0-Tst.xbar)^2/Tst.Sxx));
    extrap  = opt.RefV0 < min(vm) || opt.RefV0 > max(vm);

    % Ambroso length scale, and the KD velocity/time scales built from it
    vStar = sqrt(d0A*980);  tStar = sqrt(d0A/980);

    rows{end+1} = table(conds(c), height(S), nH, S.rho_g(1), ...
        median(S.v0_ratio), ...
        d0A, rmseA, R2A, ...
        L.b(1), L.ci(1,1), L.ci(1,2), L.b(2), L.ci(2,1), L.ci(2,2), L.R2, L.rmse, ...
        F2stat, local_fcrit(nH-2), F2stat > local_fcrit(nH-2), ...
        P.b(2), P.ci(2,1), P.ci(2,2), rmseP, ...
        Tst.b(2), Tst.ci(2,1), Tst.ci(2,2), ...
        t0ref, t0ref-Tst.tc*t0refSE, t0ref+Tst.tc*t0refSE, extrap, ...
        vStar, tStar, min(vm)/vStar, ...
        local_best(rmseA, L.rmse, rmseP), ...
        'VariableNames',{'condition','n_trials','n_heights','rho_g', ...
            'v0_ratio_median', ...
            'F1_d0_cm','F1_rmse','F1_R2', ...
            'F2_d0_cm','F2_d0_lo','F2_d0_hi','F2_alpha_s','F2_a_lo','F2_a_hi', ...
            'F2_R2','F2_rmse','F_vs_const','F_crit05','vel_dep_established', ...
            'F3_n','F3_n_lo','F3_n_hi','F3_rmse', ...
            'dtstop_dv0','dt_lo','dt_hi', ...
            't0_at_ref_s','t0_ref_lo','t0_ref_hi','t0_extrapolated', ...
            'vStar_cm_s','tStar_s','v0min_over_vStar','best_form'}); %#ok<AGROW>
end
FIT = vertcat(rows{:});
nc = height(FIT);

% ---------------------------------------------- 4. alpha density check
iRef = find(FIT.condition==opt.RefCond,1);
FIT.F2_alpha_pred = nan(nc,1);
if ~isempty(iRef)
    FIT.F2_alpha_pred = FIT.F2_alpha_s(iRef) * sqrt(FIT.rho_g(iRef)./FIT.rho_g);
end
FIT.alpha_obs_over_pred = FIT.F2_alpha_s ./ FIT.F2_alpha_pred;

% Goldman & Umbanhowar Eq.(8): t0 ~ (rho_s/rho_g)^(1/4) sqrt(R/g).
% Fixed foot -> t0 ~ rho_g^(-1/4). Evaluated at the COMMON reference speed so
% the comparison does not depend on each condition's v0 range.
FIT.t0_pred_rho14 = nan(nc,1);
if ~isempty(iRef)
    FIT.t0_pred_rho14 = FIT.t0_at_ref_s(iRef) * (FIT.rho_g(iRef)./FIT.rho_g).^(1/4);
end
FIT.t0_obs_over_pred = FIT.t0_at_ref_s ./ FIT.t0_pred_rho14;

% ------------------------------------------------------------ 5. console
fprintf('\n--- Head to head on identical per-height means (lower rmse wins) ---\n');
fprintf('  %-16s %-22s %-22s %-22s %s\n','condition','F1 Ambroso (1 par)', ...
        'F2 linear v0 (2 par)','F3 power v0 (2 par)','best');
for c = 1:nc
    fprintf('  %-16s rmse %.3f R2 %+.2f   rmse %.3f R2 %+.2f   rmse %.3f n=%.2f     %s\n', ...
        FIT.condition(c), FIT.F1_rmse(c), FIT.F1_R2(c), ...
        FIT.F2_rmse(c), FIT.F2_R2(c), FIT.F3_rmse(c), FIT.F3_n(c), FIT.best_form(c));
end

fprintf('\n--- F2 parameters, and is velocity dependence established at all? ---\n');
for c = 1:nc
    fprintf('  %-16s d0 = %+6.3f [%+6.3f,%+6.3f] cm   alpha = %.5f [%.5f,%.5f] s\n', ...
        FIT.condition(c), FIT.F2_d0_cm(c), FIT.F2_d0_lo(c), FIT.F2_d0_hi(c), ...
        FIT.F2_alpha_s(c), FIT.F2_a_lo(c), FIT.F2_a_hi(c));
    fprintf('  %-16s F = %6.2f vs F_crit(0.05) = %5.2f  ->  %s\n', '', ...
        FIT.F_vs_const(c), FIT.F_crit05(c), ...
        local_tern(FIT.vel_dep_established(c), 'established', ...
                   'NOT established -- do not interpret d0 or alpha'));
end

fprintf('\n--- Goldman & Umbanhowar: alpha ~ rho_g^(-1/2)  (ref %s) ---\n', opt.RefCond);
for c = 1:nc
    fprintf('  %-16s alpha obs %.5f   pred %.5f   obs/pred = %.2f %s\n', ...
        FIT.condition(c), FIT.F2_alpha_s(c), FIT.F2_alpha_pred(c), ...
        FIT.alpha_obs_over_pred(c), ...
        local_tern(abs(FIT.alpha_obs_over_pred(c)-1) > 0.3, '  <-- departs', ''));
end

fprintf('\n--- t_stop vs v0, physical units, single slope (no plateau fit) ---\n');
for c = 1:nc
    sig = 'flat / not resolved';
    if FIT.dt_hi(c) < 0, sig = 'DECREASES (matches KD)';
    elseif FIT.dt_lo(c) > 0, sig = 'increases (opposite)'; end
    fprintf('  %-16s dt_stop/dv0 = %+9.2e [%+9.2e,%+9.2e]   %s\n', ...
        FIT.condition(c), FIT.dtstop_dv0(c), FIT.dt_lo(c), FIT.dt_hi(c), sig);
end

fprintf('\n--- KD scales from d0, and where the data sit relative to them ---\n');
for c = 1:nc
    fprintf('  %-16s d0 = %.3f cm  ->  v* = %5.1f cm/s, t* = %5.1f ms   lowest v0 = %.1f x v*\n', ...
        FIT.condition(c), FIT.F1_d0_cm(c), FIT.vStar_cm_s(c), 1000*FIT.tStar_s(c), ...
        FIT.v0min_over_vStar(c));
end
fprintf(['  KD find t_stop rising steeply below v0 ~ v*, then levelling off. Every\n' ...
         '  condition here starts well above v*, so that low-speed rise is OUTSIDE\n' ...
         '  the measured range and is not tested by these data.\n']);

fprintf('\n--- G&U Eq.(8): t0 ~ rho_g^(-1/4), evaluated at v0 = %g cm/s ---\n', opt.RefV0);
for c = 1:nc
    fprintf('  %-16s t0 obs %5.1f ms [%5.1f,%5.1f]   pred %5.1f ms   obs/pred %.2f%s%s\n', ...
        FIT.condition(c), 1000*FIT.t0_at_ref_s(c), 1000*FIT.t0_ref_lo(c), ...
        1000*FIT.t0_ref_hi(c), 1000*FIT.t0_pred_rho14(c), FIT.t0_obs_over_pred(c), ...
        local_tern(abs(FIT.t0_obs_over_pred(c)-1) > 0.3, '  <-- departs', ''), ...
        local_tern(FIT.t0_extrapolated(c), '  [EXTRAPOLATED past this v0 range]', ''));
end

fprintf('\n--- v0 reliability (F2 and F3 depend on v0; F1 does not) ---\n');
for c = 1:nc
    fprintf('  %-16s median v0/v_freefall = %.2f\n', FIT.condition(c), FIT.v0_ratio_median(c));
end

% ------------------------------------------------------------- 6. figure
col = containers.Map({'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {[0 .45 .74],[.30 .75 .93],[.85 .33 .10],[.64 .08 .18]});
mk = containers.Map({'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {'o','s','^','d'});
isGB = @(s) startsWith(s,"GB");

fig = figure('Color','w','Position',[60 60 1180 470]);
tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% A: F2 linear form
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
hA = gobjects(nc,1); lA = strings(nc,1);
for c = 1:nc
    k = char(FIT.condition(c)); cc = col(k);
    S = DATA(DATA.condition==FIT.condition(c),:);
    M = MEANS(MEANS.condition==FIT.condition(c),:);
    plot(ax1,S.v0_cm_s,S.d_cm,mk(k),'LineStyle','none','MarkerSize',3.5, ...
        'MarkerEdgeColor',cc,'MarkerFaceColor','none','HandleVisibility','off');
    errorbar(ax1,M.v0_mean,M.d_mean,M.d_sd,'LineStyle','none','Color',cc, ...
        'LineWidth',1.1,'CapSize',3,'HandleVisibility','off');
    fc = local_tern(isGB(FIT.condition(c)),cc,[1 1 1]);
    hA(c) = plot(ax1,M.v0_mean,M.d_mean,mk(k),'LineStyle','none','MarkerSize',7, ...
        'MarkerEdgeColor',cc,'MarkerFaceColor',fc,'LineWidth',1.2);
    xx = linspace(min(M.v0_mean),max(M.v0_mean),30);
    plot(ax1,xx,FIT.F2_d0_cm(c)+FIT.F2_alpha_s(c)*xx, ...
        local_tern(isGB(FIT.condition(c)),'-','--'),'Color',cc,'LineWidth',1.6, ...
        'HandleVisibility','off');
    lA(c) = sprintf('%s  d_0=%.2f, \\alpha=%.4f%s',FIT.condition(c), ...
        FIT.F2_d0_cm(c),FIT.F2_alpha_s(c), ...
        local_tern(FIT.vel_dep_established(c),'',' (n.s.)'));
end
xlabel(ax1,'impact velocity  v_0  (cm s^{-1})'); ylabel(ax1,'penetration depth  d  (cm)');
title(ax1,'A.  de Bruyn-Walsh / Goldman-Umbanhowar:  d = d_0 + \alphav_0');
legend(ax1,hA,cellstr(lA),'Location','northwest','Box','off','FontSize',8);

% B: stopping time, PHYSICAL units, with the KD reference scales drawn on
% (KD plot t_stop [s] vs v0 [cm/s] and mark the two projectile-size scales as
% lines -- they do not rescale the axes, and neither do we).
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
hB = gobjects(nc,1); lB = strings(nc,1);
for c = 1:nc
    k = char(FIT.condition(c)); cc = col(k);
    M = MEANS(MEANS.condition==FIT.condition(c),:);
    errorbar(ax2,M.v0_mean,M.tstop_mean*1000,M.tstop_sd*1000,'LineStyle','none', ...
        'Color',cc,'LineWidth',1.1,'CapSize',3,'HandleVisibility','off');
    fc = local_tern(isGB(FIT.condition(c)),cc,[1 1 1]);
    hB(c) = plot(ax2,M.v0_mean,M.tstop_mean*1000,mk(k),'LineStyle','none', ...
        'MarkerSize',7,'MarkerEdgeColor',cc,'MarkerFaceColor',fc,'LineWidth',1.2);
    xx = linspace(min(M.v0_mean),max(M.v0_mean),30);
    b0 = mean(M.tstop_mean) - FIT.dtstop_dv0(c)*mean(M.v0_mean);
    plot(ax2,xx,(b0+FIT.dtstop_dv0(c)*xx)*1000, ...
        local_tern(isGB(FIT.condition(c)),'-','--'),'Color',cc,'LineWidth',1.6, ...
        'HandleVisibility','off');
    lB(c) = sprintf('%s   t_0(%g) = %.0f ms%s', FIT.condition(c), opt.RefV0, ...
        1000*FIT.t0_at_ref_s(c), local_tern(FIT.t0_extrapolated(c),' *',''));
end
% KD reference scales, built from the reference condition's d0
if ~isempty(iRef)
    yline(ax2, 1000*FIT.tStar_s(iRef), 'k:', ...
        sprintf('(d_0/g)^{1/2} = %.0f ms', 1000*FIT.tStar_s(iRef)), ...
        'FontSize',8,'LabelHorizontalAlignment','right');
    xline(ax2, FIT.vStar_cm_s(iRef), 'k--', ...
        sprintf('(d_0g)^{1/2} = %.0f cm/s', FIT.vStar_cm_s(iRef)), ...
        'FontSize',8,'LabelOrientation','horizontal');
end
xline(ax2, opt.RefV0, 'Color',[.6 .6 .6], 'LineStyle','-', ...
    'Label',sprintf('t_0 ref %g', opt.RefV0), 'FontSize',8, ...
    'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xlim(ax2,[0 max(MEANS.v0_mean)*1.08]);
xlabel(ax2,'impact velocity  v_0  (cm s^{-1})');
ylabel(ax2,'stopping time  t_{stop}  (ms)');
title(ax2,'B.  t_{stop} vs v_0  (KD Fig. 1b construction)');
legend(ax2,hB,cellstr(lB),'Location','northeast','Box','off','FontSize',8);

title(tl,sprintf('Competing published forms  |  %s  |  %d trials, %d height means', ...
    FOOT,height(DATA),height(MEANS)),'FontWeight','bold');

% -------------------------------------------------------------- 7. save
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
st = char(datetime('now','Format','yyyyMMdd_HHmmss'));
p = @(s) fullfile(opt.OutDir,sprintf(s,st));
exportgraphics(fig,p('depth_scaling_forms_%s.png'),'Resolution',200);
savefig(fig,p('depth_scaling_forms_%s.fig'));
writetable(FIT,  p('depth_scaling_forms_fits_%s.csv'));
writetable(MEANS,p('depth_scaling_forms_means_%s.csv'));
fprintf('\nwrote:\n  %s\n  %s\n  %s\n  %s\n\n', ...
    p('depth_scaling_forms_%s.png'),p('depth_scaling_forms_%s.fig'), ...
    p('depth_scaling_forms_fits_%s.csv'),p('depth_scaling_forms_means_%s.csv'));
end

% ------------------------------------------------------------------ helpers
function S = local_ols(x,y)
x=x(:); y=y(:); ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
n=numel(x); X=[ones(n,1),x];
if n<3, S.b=[NaN;NaN]; S.se=[NaN;NaN]; S.ci=nan(2,2); S.R2=NaN; S.rmse=NaN;
    S.n=n; S.xbar=NaN; S.Sxx=NaN; S.s2=NaN; S.tc=NaN; return; end
b=X\y; r=y-X*b; s2=(r'*r)/(n-2); se=sqrt(diag(s2*((X'*X)\eye(2))));
tc=local_t95(n-2);
S.b=b; S.se=se; S.ci=[b-tc*se, b+tc*se];
S.R2=1-(r'*r)/sum((y-mean(y)).^2); S.rmse=sqrt(s2);
S.n=n; S.xbar=mean(x); S.Sxx=sum((x-mean(x)).^2); S.s2=s2; S.tc=tc;
end
function t=local_t95(dof), try, t=tinv(0.975,dof); catch, t=1.96; end, end
function f=local_fcrit(dof2)
% F(1, dof2) upper 5% critical value
tbl=[1 161;2 18.5;3 10.13;4 7.71;5 6.61;6 5.99;7 5.59;8 5.32;9 5.12;10 4.96; ...
     12 4.75;14 4.60;16 4.49;18 4.41;20 4.35;25 4.24;30 4.17];
if dof2>=30, f=4.0; else, f=interp1(tbl(:,1),tbl(:,2),max(dof2,1),'linear','extrap'); end
end
function s=local_best(a,b,c)
[~,i]=min([a b c]); n=["F1 Ambroso","F2 linear","F3 power"]; s=n(i);
end
function v=local_num(T,name)
v=NaN;
if ismember(name,T.Properties.VariableNames)
    x=T.(name)(1);
    if isnumeric(x)||islogical(x), v=double(x); else, v=str2double(string(x)); end
end
end
function s=local_str(T,name)
s=""; if ismember(name,T.Properties.VariableNames), s=strtrim(string(T.(name)(1))); end
end
function s=local_tern(c,a,b), if c, s=a; else, s=b; end, end
