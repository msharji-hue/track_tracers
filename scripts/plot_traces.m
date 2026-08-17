function OUT = plot_traces(root, varargin)
%PLOT_TRACES  Representative z(t) and v(t), 2x1 per condition. READ-ONLY.
%
%   One figure per condition: z(t) on top, v(t) below, sharing a time axis
%   with t = 0 at impact. Time limits are common to EVERY figure so panels can
%   be laid side by side on a slide and read directly against each other.
%
%   SIGN CONVENTION. Positive UP, negative down -- matching Katsuragi & Durian
%   (2007) Fig. 1, where z runs 0 to -8 cm and v runs -400 to 0 cm/s. So z is
%   negative below the impact point and v is negative while the foot falls,
%   and no axis reversal is needed: down on the page is down in the bed.
%   Orientation is derived from each trace at run time, not assumed, so a
%   flipped axis upstream cannot silently invert a panel.
%
%   HOW THE THREE TRIALS ARE CHOSEN -- and why it matters
%   Each drop height holds 8-11 replicates scattering ~19% in depth, so
%   "nearest v0 to the target" draws an essentially RANDOM member of that
%   group; it produced a 165 mm trial shallower than a 65 mm one. Instead:
%   collect every trial within a tolerance of the target speed, then take the
%   one whose penetration depth is closest to the MEDIAN of that group. The
%   result is representative rather than an arbitrary draw. No trial may serve
%   two velocity levels.
%
%   WHY THE DEFAULT TARGETS ARE 90 / 120 / 190 cm/s
%   An earlier set (110/145/185) put the middle target in a gap: GB/shallow and
%   CHIN/as_poured each had exactly ONE trial within tolerance of 145 cm/s, and
%   that lone trial happened to be shallow, so those panels showed the medium
%   speed penetrating less than the low speed. That was a sampling artefact of
%   the target choice, not a property of the substrate. With 90/120/190 at a
%   +/-25 cm/s window every condition has >=8 trials per level (dense: 31/18/4)
%   and the MEDIAN depth rises with speed in all four:
%       GB/full        1.69 -> 1.86 -> 2.51 cm
%       GB/shallow     1.28 -> 1.84 -> 2.26 cm
%       CHIN/as_poured 1.59 -> 1.63 -> 2.41 cm
%       CHIN/dense     0.28 -> 0.76 -> 0.84 cm
%   Contrast drops from 2.5x to 2.1x, which is the price of sampling all four.
%
%   MONOTONIC SELECTION -- and its honest limits
%   Even with monotonic group medians, a single draw can still come out
%   ordered wrongly, because depth scatters ~19% within a level. With
%   'Monotonic' true the routine searches combinations of candidate trials and
%   keeps the one that (a) has depth increasing with level and (b) sits closest
%   to each level's median depth. These are therefore REPRESENTATIVE traces
%   chosen for clarity, not a random sample: the figure should be captioned as
%   such, and the real trial-to-trial spread belongs in the fit figures where
%   it is quantified. If no ordered combination exists the routine says so and
%   falls back to median-representative selection rather than failing quietly.
%
%   MATCHING VARIABLE -- default 'v0', but read this
%   'v0' matches impact speed across conditions, so every legend reads the
%   same rounded value and the panels answer "at the same impact speed, how do
%   substrates differ?". The cost: the v0 window common to all four conditions
%   is only 101-190 cm/s, a 1.7x span.
%
%   'height' matches TRUE drop height instead -- the controlled variable, exact
%   and free of impact-detection error. Heights 25/65/165 mm are common to all
%   four conditions and deliver a 4x span in v0. On GB/full that makes the
%   Katsuragi & Durian (2007) stopping-time signature visible in three traces:
%   t_stop falls 41.8 -> 35.3 -> 18.8 ms as impact speed rises. v0-matching
%   cannot show this, because its speed range is too narrow. Use 'height' when
%   the point is how the response varies with drop energy; use 'v0' when the
%   point is a like-for-like substrate comparison.
%
%   'height' is also the only fair option for CHIN/dense, whose v0 does not
%   track drop height -- matching that condition by v0 can pair a 285 mm drop
%   against a 25 mm one.
%
%   STOP MARKERS use the same symbol per condition as the depth-scaling
%   figures: o GB/full, s GB/shallow, ^ CHIN/as_poured, d CHIN/dense; filled
%   for glass beads, open for chinchilla.
%
%   SMOOTHING is cosmetic and display-only. A short moving mean, not
%   Savitzky-Golay or a polynomial: this pipeline has already had a polynomial
%   smoother removed for erasing the stop discontinuity and SG retired for the
%   same reason. At ~2778 fps a 9-frame window is 3.2 ms against a 19-42 ms
%   stopping event. Pass 'ShowRaw' true to check that for yourself.
%
%   USAGE
%       OUT = plot_traces('F:\ME_GRANULAB\JerboaImpact');
%       OUT = plot_traces(root, 'MatchBy','height');            % see above
%       OUT = plot_traces(root, 'ShowRaw', true);
%
%   OPTIONS
%       'MatchBy'     'v0' (default) or 'height'
%       'Targets'     [] uses [110 145 185] cm/s for v0, [25 65 165] mm for height
%       'Conditions'  default all four
%       'TolV0'       candidate window for v0 matching, cm/s (default 25)
%       'ClipVAtZero' truncate each v(t) curve at its first zero crossing
%                     after impact, so no negative (rebound / post-stop noise)
%                     velocities are drawn. Pre-impact data is kept. Default true.
%       'TargetWeight' weight on staying near the target speed relative to
%                     being depth-representative (default 2). With 0 the
%                     chosen trial can drift the full tolerance from target --
%                     that produced a curve labelled 190 cm/s from a 174 cm/s
%                     trial, and two 'different' levels only 7 cm/s apart.
%       'LabelMode'   'actual' (default) prints the measured v0 rounded to the
%                     nearest 5; 'target' prints the nominal target. 'actual'
%                     is the default because a label must not claim a speed the
%                     trial does not have.
%       'Monotonic'   require depth to increase with velocity level (default
%                     true). See the note below before turning this off.
%       'SmoothN'     moving-mean window in frames (default 9)
%       'ShowRaw'     overlay unsmoothed traces (default false)
%       'PadBefore'   ms of approach shown before impact (default 5). Applies
%                     to both panels so their time axes stay aligned.
%       'ClampV'      hold the velocity axis at v >= 0 (default true). Any
%                     rebound after the stop is sub-millimetre, so the small
%                     negative excursion it produces is display noise here.
%       'OutDir'      default <root>/03_RESULTS/_batch_logs

opt.MatchBy    = 'v0';
opt.Monotonic  = true;
opt.TargetWeight = 2;      % how hard to pull the chosen trial toward its target v0
opt.LabelMode  = 'actual'; % 'actual' (rounded measured v0) or 'target'
opt.ClipVAtZero = true;
opt.Targets    = [];
opt.Conditions = ["GB/full","GB/shallow","CHIN/as_poured","CHIN/dense"];
opt.TolV0      = 25;
opt.SmoothN    = 9;
opt.ShowRaw    = false;
opt.ClampV     = true;
opt.PadBefore  = 5;    % ms of approach shown before impact
opt.PadAfter   = 0.30;
opt.OutDir     = fullfile(root,'03_RESULTS','_batch_logs');
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

byV0 = strcmpi(opt.MatchBy,'v0');
if isempty(opt.Targets)
    if byV0, opt.Targets = [90 120 190]; else, opt.Targets = [25 65 165]; end
end
nT = numel(opt.Targets);

% colour is tied to LEVEL, identical in every figure
lvlCol = [0.00 0.45 0.74;   % low  - blue
          0.93 0.69 0.13;   % mid  - amber
          0.80 0.15 0.15];  % high - red
if nT > 3, lvlCol = lines(nT); end
lvlName = ["low","medium","high"];

fprintf('\n=== plot_traces  (READ-ONLY) ===\n');
fprintf('matching on %s, targets %s %s\n', opt.MatchBy, ...
    strjoin(compose('%g',opt.Targets(:)'), ' / '), local_tern(byV0,'cm/s','mm'));

% ---------------------------------------------------------------- index
F = dir(fullfile(root,'03_RESULTS','**','*_kin.mat')); F = F(~[F.isdir]);
if isempty(F), error('plot_traces:noFiles','No *_kin.mat under %s', root); end
n = numel(F);
tag=strings(n,1); cnd=strings(n,1); v0=nan(n,1); hh=nan(n,1); dd=nan(n,1);
ts=nan(n,1); pth=strings(n,1);
% Some *_kin.mat in the results tree are from other experiments and carry no
% meta/kin (e.g. the *_long files). Silence load's per-file warning and skip
% them, then report how many were skipped rather than hiding it.
ws = warning('off','MATLAB:load:variableNotFound');
cleanupWarn = onCleanup(@() warning(ws));
nSkip = 0; skipped = strings(0);
for i = 1:n
    p = fullfile(F(i).folder, F(i).name);
    S = load(p,'meta','kin');
    if ~isfield(S,'meta') || ~isfield(S,'kin')
        nSkip = nSkip + 1;
        skipped(end+1) = string(erase(F(i).name,'_kin.mat')); %#ok<AGROW>
        continue
    end
    tag(i) = string(S.meta.trialTag);
    cnd(i) = string(S.meta.material) + "/" + string(S.meta.container);
    v0(i)  = abs(S.kin.v0_cm_s);
    dd(i)  = S.kin.d_final_cm;
    ts(i)  = S.kin.t_stop_s * 1000;
    hh(i) = S.meta.dropHeight_mm;
    pth(i) = p;
end
ok = tag~="" & isfinite(v0) & isfinite(dd) & isfinite(hh);
if nSkip > 0
    fprintf('skipped %d file(s) with no meta/kin: %s\n', nSkip, strjoin(skipped, ', '));
end
fprintf('usable trials: %d\n', sum(ok));

% ------------------------------------------------- representative trials
conds = opt.Conditions(ismember(opt.Conditions, unique(cnd(ok))));
rows = {};
for c = 1:numel(conds)
    m = ok & cnd == conds(c);

    % candidate pool per level, and that pool's median depth
    pool = cell(nT,1); medD = nan(nT,1);
    for k = 1:nT
        avail = find(m);
        if byV0
            near = avail(abs(v0(avail)-opt.Targets(k)) <= opt.TolV0);
        else
            near = avail(abs(hh(avail)-opt.Targets(k)) <= 1);
        end
        if isempty(near)                   % target outside range: nearest one
            key = local_tern(byV0, v0, hh);
            [~,j] = min(abs(key(avail)-opt.Targets(k)));
            near = avail(j);
        end
        % keep the 25 closest to the level median so the search stays small
        medD(k) = median(dd(near));
        [~,srt] = sort(abs(dd(near)-medD(k)));
        pool{k} = near(srt(1:min(25,numel(srt))));
    end

    % Stopping time: Katsuragi & Durian report t_stop FALLING as impact speed
    % rises. Enforce that ordering only where this condition's own group
    % medians already show it -- otherwise the selection would be hunting for
    % atypical trials to display a trend the condition does not have. For
    % CHIN/as_poured the median t_stop actually RISES with speed and the
    % fitted trend is not significant, so the constraint is skipped there and
    % the panel shows the behaviour as measured.
    medT = arrayfun(@(k) median(ts(pool{k})), 1:nT);
    wantTs = all(diff(medT) < 0);

    % choose one trial per level: distinct, depth increasing with level,
    % t_stop decreasing where supported, each close to its level's median
    % depth AND close to its target speed
    best = []; bestCost = inf; nOrdered = 0;
    idxGrid = local_combos(cellfun(@numel,pool));
    for r = 1:size(idxGrid,1)
        pick = arrayfun(@(k) pool{k}(idxGrid(r,k)), 1:nT);
        if numel(unique(pick)) < nT, continue; end
        if opt.Monotonic
            if ~all(diff(dd(pick)) > 0), continue; end
            if wantTs && ~all(diff(ts(pick)) < 0), continue; end
            nOrdered = nOrdered + 1;
        end
        key = local_tern(byV0, v0, hh);
        cost = sum(abs(dd(pick)-medD(:)) ./ medD(:)) ...
             + opt.TargetWeight * sum(abs(key(pick)-opt.Targets(:)) ./ opt.Targets(:));
        if cost < bestCost, bestCost = cost; best = pick; end
    end
    if wantTs
        fprintf('  %s: t_stop ordering enforced (group medians support it)\n', conds(c));
    else
        fprintf('  %s: t_stop ordering NOT enforced (group medians do not decrease)\n', conds(c));
    end
    monoUsed = true;
    if isempty(best)                       % no ordered set exists -- say so
        monoUsed = false;
        fprintf(['  NOTE: %s has no depth-ordered combination at these targets;\n' ...
                 '        falling back to median-representative selection.\n'], conds(c));
        best = arrayfun(@(k) pool{k}(1), 1:nT);
    end

    for k = 1:nT
        idx = best(k);
        rows{end+1} = table(conds(c), k, lvlName(min(k,3)), opt.Targets(k), ...
            tag(idx), v0(idx), hh(idx), dd(idx), numel(pool{k}), monoUsed, ...
            pth(idx), ...
            'VariableNames',{'condition','level','levelName','target', ...
                'trialTag','v0_cm_s','h_true_mm','d_final_cm','nCandidates', ...
                'depthOrdered','path'}); %#ok<AGROW>
    end
end
OUT.picked = vertcat(rows{:});

fprintf('\n--- representative trials (closest to the median depth of each group) ---\n');
for c = 1:numel(conds)
    T = OUT.picked(OUT.picked.condition==conds(c),:);
    fprintf('  %s\n', conds(c));
    for k = 1:height(T)
        fprintf('    %-6s %-22s v0 = %5.1f cm/s   h = %3.0f mm   d = %.2f cm   (of %d candidates)\n', ...
            T.levelName(k), T.trialTag(k), T.v0_cm_s(k), T.h_true_mm(k), ...
            T.d_final_cm(k), T.nCandidates(k));
    end
end
if any(conds=="CHIN/dense") && byV0
    fprintf(['\n  NOTE: CHIN/dense matched by v0. Its v0 does not track drop height,\n' ...
             '  so its three traces may come from very different falls. Check the\n' ...
             '  h column above; use ''MatchBy'',''height'' for a fair dense panel.\n']);
end

% -------------------------------------------------- common time window
tEnd = 0;
for i = 1:height(OUT.picked)
    S = load(OUT.picked.path(i),'kin');   %#ok<*LOAD>
    t = (S.kin.t_s - S.kin.t_s(S.kin.impact_index))*1000;
    tEnd = max(tEnd, t(min(S.kin.stopFrame,numel(t))));
end
tEnd = tEnd*(1+opt.PadAfter);
fprintf('\ncommon time axis: %.0f to %.0f ms\n', -opt.PadBefore, tEnd);

% --------------------------------------------------------------- figures
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
st = char(datetime('now','Format','yyyyMMdd_HHmmss'));
OUT.files = strings(0);

for c = 1:numel(conds)
    T = OUT.picked(OUT.picked.condition==conds(c),:);
    key = char(conds(c));

    fig = figure('Color','w','Position',[80 80 620 700]);
    tl  = tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');

    axZ = nexttile(tl); hold(axZ,'on'); grid(axZ,'on'); box(axZ,'on');
    axV = nexttile(tl); hold(axV,'on'); grid(axV,'on'); box(axV,'on');
    hZ = gobjects(height(T),1); lbl = strings(height(T),1);

    for k = 1:height(T)
        S = load(T.path(k),'kin'); kin = S.kin;
        iImp = kin.impact_index; iStp = min(kin.stopFrame, numel(kin.t_s));
        t = (kin.t_s - kin.t_s(iImp))*1000;

        % Sign convention applied as a VARIABLE TRANSFORMATION; the axis is
        % left as a normal MATLAB axis. z = 0 at the surface, z < 0 above it,
        % z > 0 into the bed. Where the saved kinematics store penetration as
        % negative this is simply z_plot = -z_cm. Direction is read from each
        % trace rather than assumed, so an upstream sign change cannot
        % silently invert a panel.
        z = kin.z - kin.z(iImp);
        if median(z(iImp:iStp),'omitnan') < 0, z = -z; end
        v = kin.v;
        if median(v(max(1,iImp-5):iImp),'omitnan') < 0, v = -v; end

        zs = smoothdata(z,'movmean',opt.SmoothN);
        vs = smoothdata(v,'movmean',opt.SmoothN);
        w  = t >= -opt.PadBefore & t <= tEnd;   % z keeps the approach
        wv = t >= 0 & t <= tEnd;                % v starts at impact
        cc = lvlCol(min(k,size(lvlCol,1)),:);

        % Velocity window: keep everything before impact, then stop at the
        % first point where the smoothed speed reaches zero. Beyond that the
        % trace is rebound and post-stop noise, which appears as negative
        % values under a positive-down convention. Defined BEFORE any plotting
        % so the raw overlay uses the same window as the smoothed curve.
        wv = w;
        if opt.ClipVAtZero
            idxAll = (1:numel(vs))';
            iz = find(vs(:) <= 0 & idxAll > iImp, 1);
            if ~isempty(iz), wv = w & (idxAll <= iz); end
        end

        if opt.ShowRaw
            plot(axZ,t(w),z(w),'-','Color',[cc 0.20],'LineWidth',0.5,'HandleVisibility','off');
            plot(axV,t(wv),v(wv),'-','Color',[cc 0.20],'LineWidth',0.5,'HandleVisibility','off');
        end
        hZ(k) = plot(axZ,t(w),zs(w),'-','Color',cc,'LineWidth',1.8);
        plot(axV,t(wv),vs(wv),'-','Color',cc,'LineWidth',1.8);

        % legend shows the ROUNDED target so labels match across conditions
        if byV0
            if strcmpi(opt.LabelMode,'target')
                lbl(k) = sprintf('v_0 \\approx %g cm s^{-1}', T.target(k));
            else
                lbl(k) = sprintf('v_0 \\approx %g cm s^{-1}', 5*round(T.v0_cm_s(k)/5));
            end
        else
            lbl(k) = sprintf('h = %g mm  (v_0 = %.0f cm s^{-1})', T.target(k), T.v0_cm_s(k));
        end
    end

    for ax = [axZ axV]
        if opt.PadBefore > 0
            xline(ax,0,'k-','LineWidth',0.9,'HandleVisibility','off');
        end
        xlim(ax,[-opt.PadBefore tEnd]);
    end
    yline(axV,0,'k:','HandleVisibility','off');
    if opt.ClipVAtZero
        yl = ylim(axV); ylim(axV,[0 yl(2)]);
    end
    if opt.ClampV
        yl = ylim(axV);
        ylim(axV,[0 yl(2)]);
    end

    ylabel(axZ,'z / penetration depth  (cm)');
    ylabel(axV,'v  (cm s^{-1})');
    xlabel(axV,'time from impact  t  (ms)');
    set(axZ,'XTickLabel',[]);
    legend(axZ,hZ,cellstr(lbl),'Location','southwest','Box','off','FontSize',9);
    title(tl, sprintf('%s   -   \\pm7.92\\circ foot', conds(c)), ...
        'FontWeight','bold','Interpreter','tex');
    subtitle(tl, sprintf('t = 0 at impact; moving mean %d frames (display only)', ...
        opt.SmoothN), 'FontSize',8);

    fn = fullfile(opt.OutDir, sprintf('traces_%s_%s.png', ...
        matlab.lang.makeValidName(key), st));
    exportgraphics(fig, fn, 'Resolution', 200);
    savefig(fig, strrep(fn,'.png','.fig'));
    OUT.files(end+1) = string(fn); %#ok<AGROW>
end

csvp = fullfile(opt.OutDir, sprintf('traces_selected_%s.csv', st));
writetable(OUT.picked(:,1:10), csvp);   % path column omitted
fprintf('\nwrote %d figures + %s\n\n', numel(OUT.files), csvp);
end

% ------------------------------------------------------------------ helpers
function G = local_combos(sizes)
%LOCAL_COMBOS  All index combinations across pools, capped for speed.
sizes = min(sizes(:)', 25);
v = arrayfun(@(n) 1:n, sizes, 'UniformOutput', false);
[g{1:numel(sizes)}] = ndgrid(v{:});
G = cell2mat(cellfun(@(x) x(:), g, 'UniformOutput', false));
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
