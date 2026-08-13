function [FIT, MEANS, DATA] = ambroso_collapse(root, varargin)
%AMBROSO_COLLAPSE  Penetration collapse onto one curve, driven by MEASURED v0.
%
%   Ambroso, Santore, Abate & Durian (2005), Phys. Rev. E 71, 051305, Fig. 4;
%   same form in Katsuragi & Durian (2007), Nat. Phys. 3, 420, Fig. 2a.
%
%   WHY v0 AND NOT DROP DISTANCE
%   The papers write the law in the total drop distance H = h + d, which
%   presumes a frictionless fall so that the release height alone sets the
%   energy delivered. The carriage here runs on a rail whose friction has not
%   been measured, so h is a control setting rather than a measured input.
%   Substituting the free-fall relation H = v0^2/2g puts the same law entirely
%   in terms of the MEASURED impact speed:
%
%       d  = ( d0^2 v0^2 / 2g )^(1/3)          [equivalently d = K v0^(2/3)]
%       d0 = sqrt(2g) * d^(3/2) / v0
%       d / d0 = ( v0 / v* )^(2/3),   v* = sqrt(2 g d0)
%       d / v0^(2/3) = ( d0^2 / 2g )^(1/3) = constant if the law holds
%
%   Nothing but the driving variable has changed; the exponent in v0 is 2/3
%   because 1/3 in distance is 2/3 in speed.
%
%   ONE CONSEQUENCE, STATED PLAINLY
%   In the original form d0 is exactly the penetration at zero drop height,
%   because the +d inside H survives as v0 -> 0. Dropping to pure v0 removes
%   that term, so here d0 is a characteristic length fitted from d and v0
%   rather than the h = 0 depth itself. The two agree closely when the
%   penetration is small next to the fall, which holds for most of this data,
%   but the h = 0 measurement now tests the law approximately rather than
%   exactly. Pass 'KeepDropTerm' true to restore d0's exact meaning by adding
%   the measured penetration back, v0^2/2g + d, still with no reliance on h.
%
%   USAGE
%       [FIT, MEANS, DATA] = ambroso_collapse('D:\ME_GRANULAB\JerboaImpact');
%
%   OPTIONS
%       'OutDir'       default <root>/03_RESULTS/_batch_logs
%       'MinRep'       min trials per drop height to form a mean (default 3)
%       'PoolGB'       also fit GB/full + GB/shallow pooled (default true)
%       'Conditions'   restrict to these conditions ([] = all present)
%       'Select'       list of trialTags to restrict to. Use this to drive the
%                      fit from a curated set, e.g. the cleaned Default GB/full
%                      dataset returned by load_default_gb.
%       'KeepDropTerm' add the measured penetration to the fall distance
%                      (default false); see the note above
%       'ScaleA'       panel A axes, 'log' (default) or 'linear'. Log axes
%                      show the shape of the law directly -- depth rising and
%                      flattening -- rather than straightening it into a line.
%                      Panel B stays logarithmic: a collapse onto a single
%                      power law is conventionally read on log axes.

opt.OutDir       = fullfile(root,'03_RESULTS','_batch_logs');
opt.MinRep       = 3;
opt.PoolGB       = true;
opt.Conditions   = [];
opt.KeepDropTerm = false;
opt.Select       = [];   % list of trialTags to restrict to ([] = all found)
opt.ScaleA       = 'log';      % panel A axes: 'log' (default) or 'linear'
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

G = 980;
FOOT = '\pm7.92\circ foot';
fprintf('\n=== ambroso_collapse | %s ===\n', FOOT);

% ------------------------------------------------------------------ load
F = dir(fullfile(root,'03_RESULTS','**','*_kin_scalars.csv')); F = F(~[F.isdir]);
if isempty(F), error('ambroso_collapse:noFiles','No *_kin_scalars.csv found.'); end
n = numel(F);
tag=strings(n,1); cond=strings(n,1); d=nan(n,1); v0=nan(n,1); grpLab=nan(n,1);
for i = 1:n
    T = readtable(fullfile(F(i).folder,F(i).name));
    if height(T)<1, continue; end
    tag(i)  = string(erase(F(i).name,'_kin_scalars.csv'));
    cond(i) = local_str(T,'condition');
    d(i)    = local_num(T,'d_final_cm');
    v0(i)   = abs(local_num(T,'v0_cm_s'));
    % Drop setting, used ONLY to identify which trials are repeats of each
    % other. It never enters a fit -- every fitted quantity uses measured v0.
    grpLab(i) = local_num(T,'dropHeight_mm');
end
ok = tag~="" & isfinite(d) & d>0 & isfinite(v0) & v0>0 & isfinite(grpLab);
if ~isempty(opt.Select)
    sel = ismember(tag, string(opt.Select));
    fprintf('Select: %d of %d requested tags found\n', ...
            sum(ismember(string(opt.Select), tag)), numel(opt.Select));
    ok = ok & sel;
end
DATA = table(tag(ok),cond(ok),v0(ok),d(ok),grpLab(ok), ...
    'VariableNames',{'trialTag','condition','v0_cm_s','d_cm','repGroup'});
if ~isempty(opt.Conditions)
    DATA = DATA(ismember(DATA.condition, string(opt.Conditions)), :);
    if isempty(DATA), error('ambroso_collapse:noneKept','No trials match Conditions.'); end
end

% Fall distance implied by the MEASURED impact speed. No h anywhere.
if opt.KeepDropTerm
    fall = DATA.v0_cm_s.^2/(2*G) + DATA.d_cm;
else
    fall = DATA.v0_cm_s.^2/(2*G);
end
DATA.d0_cm = sqrt(DATA.d_cm.^3 ./ fall);        % = sqrt(2g) d^(3/2)/v0
DATA.vStar = sqrt(2*G*DATA.d0_cm);
fprintf('trials: %d   (driving variable: measured v_0)\n', height(DATA));

% -------------------------------------------------------- per-v0 means
% Replicates are grouped by their drop setting so each mean averages ~10
% repeats, as the literature figures do. Grouping by rounded v0 instead splits
% those repeats across bins and starves conditions with few settings.
order = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
pres  = unique(DATA.condition,'stable');
conds = [order(ismember(order,pres)), reshape(pres(~ismember(pres,order)),1,[])];
mrows = {};
for c = 1:numel(conds)
    S = DATA(DATA.condition==conds(c),:);
    grp = findgroups(S.repGroup);
    for k = 1:max(grp)
        g = S(grp==k,:);
        if height(g) < opt.MinRep, continue; end
        mrows{end+1} = table(conds(c), height(g), ...
            mean(g.d_cm), std(g.d_cm), mean(g.v0_cm_s), ...
            mean(g.d0_cm), std(g.d0_cm), ...
            mean(g.d_cm)/mean(g.v0_cm_s)^(2/3), ...
            std(g.d_cm)/mean(g.v0_cm_s)^(2/3), ...
            'VariableNames',{'condition','nRep','d_mean','d_sd','v0_mean', ...
                'd0_mean','d0_sd','comp','comp_sd'}); %#ok<AGROW>
    end
end
MEANS = vertcat(mrows{:});

% ---------------------------------------------------------------- fits
rows = {};
for c = 1:numel(conds)
    M = MEANS(MEANS.condition==conds(c),:);
    S = DATA(DATA.condition==conds(c),:);
    if height(M) < 3
        fprintf('  SKIPPED %s: only %d group(s) with >= %d repeats\n', ...
                conds(c), height(M), opt.MinRep);
        continue
    end
    rows{end+1} = local_assess(conds(c), height(S), M, G); %#ok<AGROW>
end
FIT = vertcat(rows{:});
if opt.PoolGB && all(ismember(["GB/full","GB/shallow"], conds))
    Mg = MEANS(startsWith(MEANS.condition,"GB"),:);
    Sg = DATA(startsWith(DATA.condition,"GB"),:);
    FIT = [FIT; local_assess("GB pooled", height(Sg), Mg, G)];
end

% -------------------------------------------------------------- console
fprintf('\n--- d0, the single length (from d and v_0 only) ---\n');
for c = 1:height(FIT)
    fprintf('  %-16s d0 = %5.3f +/- %5.3f cm   v* = %5.1f cm/s   (%2d speeds, %3d trials)\n', ...
        FIT.condition(c), FIT.d0_cm(c), FIT.d0_sd(c), FIT.vStar(c), ...
        FIT.n_speeds(c), FIT.n_trials(c));
end
fprintf('\n--- Collapse quality: deviation from d/d0 = (v_0/v*)^(2/3) ---\n');
for c = 1:height(FIT)
    fprintf('  %-16s v_0/v* spans %5.2f - %5.2f (%4.1fx)   rms deviation %5.1f%%\n', ...
        FIT.condition(c), FIT.x_min(c), FIT.x_max(c), FIT.x_max(c)/FIT.x_min(c), ...
        FIT.rms_dev_pct(c));
end
fprintf('\n--- Compensated test: d/v_0^{2/3} is constant iff the law holds ---\n');
for c = 1:height(FIT)
    fprintf('  %-16s slope on log10(v_0) = %+7.4f +/- %6.4f   %s\n', ...
        FIT.condition(c), FIT.comp_slope(c), FIT.comp_slope_se(c), ...
        local_tern(abs(FIT.comp_slope(c)) < 2*FIT.comp_slope_se(c), ...
                   'FLAT', 'SLOPED -- law does not hold here'));
end
fprintf('\n--- Free exponent, cross-check (target 2/3 = 0.667) ---\n');
for c = 1:height(FIT)
    fprintf('  %-16s n = %5.3f [%5.3f, %5.3f]%s\n', FIT.condition(c), ...
        FIT.exponent(c), FIT.exp_lo(c), FIT.exp_hi(c), ...
        local_tern(FIT.exp_incl_target(c),'   2/3 in CI',''));
end
fprintf('\n');

% --------------------------------------------------------------- figure
col = containers.Map({'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {[0 .45 .74],[.30 .75 .93],[.85 .33 .10],[.64 .08 .18]});
mk  = containers.Map({'GB/full','GB/shallow','CHIN/as_poured','CHIN/dense'}, ...
    {'o','s','^','d'});
isGB = @(s) startsWith(s,"GB");

fig = figure('Color','w','Position',[60 60 1180 500]);
tl  = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% Panel A -- physical units
ax1 = nexttile(tl); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
if strcmpi(opt.ScaleA,'log'), set(ax1,'XScale','log','YScale','log'); end
hA = gobjects(numel(conds),1); lA = strings(numel(conds),1);
for c = 1:numel(conds)
    k = char(conds(c)); cc = col(k);
    i = find(FIT.condition==conds(c),1); if isempty(i), continue; end
    d0 = FIT.d0_cm(i);
    S = DATA(DATA.condition==conds(c),:);
    M = MEANS(MEANS.condition==conds(c),:);
    errorbar(ax1,M.v0_mean,M.d_mean,M.d_sd,'LineStyle','none','Color',cc, ...
        'LineWidth',1.1,'CapSize',3,'HandleVisibility','off');
    fc = local_tern(isGB(conds(c)),cc,[1 1 1]);
    hA(c) = plot(ax1,M.v0_mean,M.d_mean,mk(k),'LineStyle','none','MarkerSize',7, ...
        'MarkerEdgeColor',cc,'MarkerFaceColor',fc,'LineWidth',1.2);
    if strcmpi(opt.ScaleA,'log')
        xx = logspace(log10(min(M.v0_mean)*0.9),log10(max(M.v0_mean)*1.1),60);
    else
        xx = linspace(0, max(DATA.v0_cm_s)*1.05, 120);
    end
    plot(ax1,xx,(d0^2*xx.^2/(2*G)).^(1/3),local_tern(isGB(conds(c)),'-','--'), ...
        'Color',cc,'LineWidth',1.5,'HandleVisibility','off');
    lA(c) = sprintf('%s   n = %.3f \\pm %.3f', conds(c), FIT.exponent(i), ...
                    (FIT.exp_hi(i)-FIT.exp_lo(i))/2);
end
if ~strcmpi(opt.ScaleA,'log')
    xlim(ax1,[0 max(DATA.v0_cm_s)*1.05]);
    ylim(ax1,[0 max(DATA.d_cm)*1.10]);
end
% Reference slope 2/3, from Newhall & Durian (2003) Eq.(4) -- the only
% published form that is an explicit power law in v0. (KD Fig. 2b also plots
% depth against v0, but fits the linear d = d0 + alpha|v0| there instead.)
iG = find(startsWith(FIT.condition,"GB"),1); if isempty(iG), iG = 1; end
MGb = MEANS(MEANS.condition==FIT.condition(iG),:);
if ~isempty(MGb)
    xg = exp(mean(log(MGb.v0_mean))); yg = exp(mean(log(MGb.d_mean)));
    xr = [min(MEANS.v0_mean)*0.9 max(MEANS.v0_mean)*1.1];
    hRefA = plot(ax1, xr, yg*(xr/xg).^(2/3), 'k--', 'LineWidth',1.3);
else
    hRefA = gobjects(1);
end
xlabel(ax1,'impact velocity  v_0   (cm s^{-1})');
ylabel(ax1,'penetration depth  d   (cm)');
title(ax1,'A.  Katsuragi & Durian (2007):  d = (d_0^2v_0^2/2g)^{1/3}');
legend(ax1,[hA(isgraphics(hA)); hRefA], ...
    [cellstr(lA(isgraphics(hA))); {'slope 2/3 (Newhall Eq. 4)'}], ...
    'Location','southeast','Box','off','FontSize',8);

% Panel B -- the collapse
ax2 = nexttile(tl); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log','YScale','log');
hB = gobjects(numel(conds),1); lB = strings(numel(conds),1);
for c = 1:numel(conds)
    k = char(conds(c)); cc = col(k);
    i = find(FIT.condition==conds(c),1); if isempty(i), continue; end
    d0 = FIT.d0_cm(i); vS = FIT.vStar(i);
    S = DATA(DATA.condition==conds(c),:);
    M = MEANS(MEANS.condition==conds(c),:);
    errorbar(ax2,M.v0_mean/vS,M.d_mean/d0,M.d_sd/d0,'LineStyle','none', ...
        'Color',cc,'LineWidth',1.1,'CapSize',3,'HandleVisibility','off');
    fc = local_tern(isGB(conds(c)),cc,[1 1 1]);
    hB(c) = plot(ax2,M.v0_mean/vS,M.d_mean/d0,mk(k),'LineStyle','none', ...
        'MarkerSize',7,'MarkerEdgeColor',cc,'MarkerFaceColor',fc,'LineWidth',1.2);
    lB(c) = sprintf('%s   rms %.0f%%', conds(c), FIT.rms_dev_pct(i));
end
xx = logspace(log10(min(FIT.x_min)*0.85),log10(max(FIT.x_max)*1.15),100);
hRef = plot(ax2,xx,xx.^(2/3),'k-','LineWidth',1.8);
xlabel(ax2,'v_0 / v^*        ( v^* = (2g d_0)^{1/2} )');
ylabel(ax2,'d / d_0');
title(ax2,'B.  Ambroso et al. (2005):  d/d_0 = (v_0/v^*)^{2/3}');
legend(ax2,[hB(isgraphics(hB)); hRef], ...
    [cellstr(lB(isgraphics(hB))); {'(v_0/v^*)^{2/3}'}], ...
    'Location','southeast','Box','off','FontSize',8);

title(tl,sprintf(['Ambroso et al. (2005) / Katsuragi & Durian (2007)  |  %s  |  ' ...
    '%d trials, %d speed means'], FOOT, height(DATA), height(MEANS)), ...
    'FontWeight','bold','FontSize',9);

% ----------------------------------------------------------------- save
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
st = char(datetime('now','Format','yyyyMMdd_HHmmss'));
p = @(s) fullfile(opt.OutDir,sprintf(s,st));
exportgraphics(fig,p('ambroso_collapse_%s.png'),'Resolution',200);
savefig(fig,p('ambroso_collapse_%s.fig'));
writetable(FIT,  p('ambroso_collapse_fits_%s.csv'));
writetable(MEANS,p('ambroso_collapse_means_%s.csv'));
fprintf('wrote:\n  %s\n  %s\n  %s\n  %s\n\n', ...
    p('ambroso_collapse_%s.png'),p('ambroso_collapse_%s.fig'), ...
    p('ambroso_collapse_fits_%s.csv'),p('ambroso_collapse_means_%s.csv'));
end

% ------------------------------------------------------------------ helpers
function R = local_assess(name, nTrials, M, G)
d0 = mean(M.d0_mean);
vS = sqrt(2*G*d0);
x  = M.v0_mean/vS;  y = M.d_mean/d0;
dev = 100*(y - x.^(2/3)) ./ x.^(2/3);
C = local_ols(log10(M.v0_mean), M.comp);            % compensated flatness
P = local_ols(log10(M.v0_mean), log10(M.d_mean));   % free exponent
R = table(string(name), nTrials, height(M), d0, std(M.d0_mean), vS, ...
    min(x), max(x), sqrt(mean(dev.^2)), mean(dev), ...
    C.b(2), C.se(2), abs(C.b(2)) < 2*C.se(2), ...
    P.b(2), P.ci(2,1), P.ci(2,2), (P.ci(2,1)<=2/3)&&(P.ci(2,2)>=2/3), ...
    'VariableNames',{'condition','n_trials','n_speeds','d0_cm','d0_sd','vStar', ...
        'x_min','x_max','rms_dev_pct','mean_dev_pct', ...
        'comp_slope','comp_slope_se','comp_flat', ...
        'exponent','exp_lo','exp_hi','exp_incl_target'});
end

function S = local_ols(x,y)
x=x(:); y=y(:); ok=isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
n=numel(x); X=[ones(n,1),x];
if n<3, S.b=[NaN;NaN]; S.se=[NaN;NaN]; S.ci=nan(2,2); return; end
b=X\y; r=y-X*b; s2=(r'*r)/(n-2); se=sqrt(diag(s2*((X'*X)\eye(2))));
try, tc=tinv(0.975,n-2); catch, tc=1.96; end
S.b=b; S.se=se; S.ci=[b-tc*se, b+tc*se];
end

function v = local_num(T,name)
v=NaN;
if ismember(name,T.Properties.VariableNames)
    x=T.(name)(1);
    if isnumeric(x)||islogical(x), v=double(x); else, v=str2double(string(x)); end
end
end
function s = local_str(T,name)
s=""; if ismember(name,T.Properties.VariableNames), s=strtrim(string(T.(name)(1))); end
end
function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
