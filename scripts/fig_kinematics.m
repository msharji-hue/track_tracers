function R = fig_kinematics(varargin)
%FIG_KINEMATICS  Step-2 kinematics figure: depth, velocity, net acceleration.
%
%   3x3. Columns are the foot models (Default, Tight, Wide); rows are depth,
%   velocity and net acceleration against time from impact. READ-ONLY: reads
%   the kinematics tree and writes only the figure files.
%
%   USAGE
%       fig_kinematics
%       fig_kinematics('Root', 'D:\ME_GRANULAB\JerboaImpact')
%       R = fig_kinematics('AllTrials', true, 'Save', false);
%
%   OPTIONS (name-value)
%       'Root'       results root (parent of 03_RESULTS); '' prompts
%       'Models'     column order (default ["Default","Tight","Wide"])
%       'AllTrials'  plot every kept trial (default false: one representative
%                    per model x height, see below)
%       'PreCapMs'   pre-impact context cap, ms (default 20)
%       'OutDir'     default <Root>/03_RESULTS/_figures
%       'Save'       write PDF + PNG (default true)
%       'Show'       display the figure (default true)
%
%   FIELD NAMES are read from kd_kinematics' save block, not assumed:
%       kin.t_s        time, ALREADY re-zeroed at impact (kd_kinematics step 5:
%                      t_s = t - t(impact_index)), so no further shift is made
%       kin.z          depth, alias of kin.depthRod_cm, zeroed at impact
%       kin.v          smoothed velocity. The local variable inside
%                      kd_kinematics is v_smooth; it is NOT saved under that
%                      name, and reading kin.v_smooth returns nothing.
%       kin.a_plus_g   net acceleration, the pipeline's own field
%       kin.impact_index, kin.stopFrame
%   Table columns (trialTag, model, dropHeight_mm, v0_cm_s, isZeroDrop,
%   kinPath) come from load_kinematics_set.
%
%   SIGN CONVENTION is pipeline-native and is NOT mirrored to KD 2007: depth
%   is positive INTO the bed (rod_displacement, signConvention +1).
%
%   ROW 3 PLOTS kin.a_plus_g EXACTLY AS THE PIPELINE COMPUTES IT, namely
%   a_plus_g = -a - g (kd_kinematics step 4). It is deliberately not
%   recomputed or re-signed here: the pipeline's convention is what every
%   number in docs/QUANTITIES.md uses, and a figure that quietly adopted a
%   different one would disagree with all of them.
%
%   !! TWO POINTS TO CONFIRM BEFORE PUBLISHING ROW 3 !!
%   Flagged rather than resolved, because neither can be checked without
%   running against real data:
%
%   (a) THERE IS NO FREE-FALL BASELINE TO SEE. kd_kinematics masks a and
%       a_plus_g to NaN outside [impact_index, stopFrame] (step 7), so row 3
%       is empty before impact however large PreCapMs is. KD 2007 Fig. 1 shows
%       a+g flat near zero during the fall; that segment does not exist in this
%       pipeline's output. Rows 1 and 2 do show the pre-impact context.
%
%   (b) FREE FALL WOULD NOT READ ZERO UNDER THIS SIGN. With depth positive into
%       the bed a free-falling rod has a = +g, so -a - g evaluates to -2g, not
%       0. Because of (a) that region is never drawn, so nothing in this figure
%       depends on it. But if row 3 is meant to be read the KD way -- free fall
%       at zero, resistance positive -- the convention needs checking against
%       kd_kinematics first. Confirm on a known trial; do not adjust this
%       script on a guess.
%
%   TIME WINDOW per trial: every available pre-impact frame up to PreCapMs
%   before impact, through stopFrame, and nothing after the stop.
%
%   ZERO-DROP TRIALS are drawn as flat lines at zero, one per model,
%   de-emphasized. load_kinematics_set quarantines their impact-derived
%   quantities -- released from contact, so there is no impact, no stop and no
%   depth -- and plotting their raw traces would present a release transient as
%   though it were an impact. A line at zero is what the protocol asserts.
%
%   COLOUR encodes v0 on one perceptually-uniform map shared by every panel,
%   with a single colorbar. Model is carried by the columns, so colouring by
%   model would spend the colour channel on information the layout already
%   gives.

% ── options ──────────────────────────────────────────────────────────────
opt.Root      = '';
opt.Models    = ["Default","Tight","Wide"];
opt.AllTrials = false;
opt.PreCapMs  = 20;
opt.OutDir    = '';
opt.Save      = true;
opt.Show      = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(thisDir), 'src'));

root = resolve_output_root(opt.Root);
if isempty(root), fprintf('Cancelled.\n'); return; end
if isempty(opt.OutDir)
    opt.OutDir = fullfile(root, '03_RESULTS', '_figures');
end
opt.Models = string(opt.Models);
nCol = numel(opt.Models);

fprintf('\n=== fig_kinematics ===\n');
fprintf('  Root     : %s\n', root);
fprintf('  Models   : %s\n', strjoin(cellstr(opt.Models), ', '));
fprintf('  Trials   : %s\n', local_tern(opt.AllTrials, 'all kept', ...
        'one per (model, height): v0 closest to the height mean'));
fprintf('  PreCapMs : %g ms of pre-impact context\n', opt.PreCapMs);
fprintf('  OutDir   : %s\n', opt.OutDir);

% ── data ─────────────────────────────────────────────────────────────────
% Manual exclusions and the zero-drop quarantine are applied by the loader, so
% every row returned is already a kept trial.
K = load_kinematics_set(root);
K = K(ismember(K.model, opt.Models), :);
if isempty(K)
    error('fig_kinematics:noTrials', 'No trials for model(s) %s under %s.', ...
          strjoin(cellstr(opt.Models), ', '), root);
end

[sel, nGroups] = local_select(K, opt);

fprintf('\n--- trials per panel ---\n');
if ~opt.AllTrials
    fprintf('  %d (model, height) group(s) -> one representative each\n', nGroups);
end
for c = 1:nCol
    nReal = sum(sel.model == opt.Models(c) & ~sel.isZeroDrop);
    nZero = sum(sel.model == opt.Models(c) &  sel.isZeroDrop);
    fprintf('  %-8s %3d curve(s)%s\n', opt.Models(c), nReal, ...
        local_tern(nZero > 0, sprintf('  + %d zero-drop (flat at 0)', nZero), ''));
end

% ── load the series ──────────────────────────────────────────────────────
% load_kinematics_set's 'series' option returns only z/v/t, not a_plus_g, so
% each trial's _kin.mat is read directly here.
nS      = height(sel);
S       = repmat(struct('t_ms',[], 'z',[], 'v',[], 'ag',[]), nS, 1);
usable  = false(nS,1);
for i = 1:nS
    if sel.isZeroDrop(i), continue; end     % flat lines, no series needed
    [S(i), usable(i)] = local_load_series(sel.kinPath(i), opt.PreCapMs);
end
nSkipped = sum(~usable & ~sel.isZeroDrop);
if nSkipped > 0
    fprintf('  %d trial(s) skipped: unreadable _kin.mat or empty time window\n', ...
            nSkipped);
end
if ~any(usable)
    error('fig_kinematics:nothingToPlot', ...
        ['No usable series. Every selected trial had an unreadable _kin.mat ' ...
         'or no frames within [-%g ms, stopFrame].'], opt.PreCapMs);
end

% ── shared colour scale ──────────────────────────────────────────────────
v0all = sel.v0_cm_s(usable);
CLIM  = [min(v0all), max(v0all)];
if diff(CLIM) <= 0, CLIM = CLIM + [-1 1]; end   % single-v0 edge case
CMAP  = parula(256);

% ── figure ───────────────────────────────────────────────────────────────
fig = figure('Color','w','Units','inches','Position',[1 1 7.0 6.4], ...
             'Visible', local_tern(opt.Show,'on','off'));
tl  = tiledlayout(fig, 3, nCol, 'Padding','compact', 'TileSpacing','compact');
colormap(fig, CMAP);

ROWS = { 'z',  'depth (cm)'      ; ...
         'v',  'v (cm/s)'        ; ...
         'ag', 'a + g (cm/s^2)'  };

AX = gobjects(3, nCol);
for r = 1:3
    for c = 1:nCol
        ax = nexttile(tl, (r-1)*nCol + c);
        hold(ax,'on'); box(ax,'on');
        AX(r,c) = ax;

        inCol = find(sel.model == opt.Models(c)).';

        % Zero-drop first so real traces draw over it. Depth and velocity are
        % flat at zero by protocol; row 3 is left empty, since a_plus_g is
        % undefined without an impact.
        if r <= 2 && any(sel.isZeroDrop(inCol))
            yline(ax, 0, '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.0);
        end

        for i = inCol
            if ~usable(i), continue; end
            y = S(i).(ROWS{r,1});
            if isempty(y) || all(isnan(y)), continue; end
            plot(ax, S(i).t_ms, y, '-', 'LineWidth', 0.8, ...
                 'Color', local_v0_color(sel.v0_cm_s(i), CLIM, CMAP));
        end

        % Sets the colorbar's scale; clim() is R2022a+, caxis works everywhere.
        caxis(ax, CLIM);

        set(ax, 'FontSize', 8, 'LineWidth', 0.5, 'Layer', 'top');
        if r == 1, title(ax, opt.Models(c), 'FontSize', 9); end
        if r == 3, xlabel(ax, 't - t_{impact} (ms)', 'FontSize', 8); end
        if c == 1, ylabel(ax, ROWS{r,2}, 'FontSize', 8); end
    end
end

% ── identical limits across the columns of each row ──────────────────────
% Applied after everything is drawn. The columns exist to be compared, and
% per-panel autoscaling would show three different y scales as though they
% were one.
for r = 1:3
    xl = [inf, -inf];  yl = [inf, -inf];
    for c = 1:nCol
        x = xlim(AX(r,c));  y = ylim(AX(r,c));
        xl = [min(xl(1), x(1)), max(xl(2), x(2))];
        yl = [min(yl(1), y(1)), max(yl(2), y(2))];
    end
    for c = 1:nCol
        xlim(AX(r,c), xl);  ylim(AX(r,c), yl);
        if c > 1, set(AX(r,c), 'YTickLabel', []); end
    end
end

cb = colorbar(AX(1,nCol));
cb.Layout.Tile   = 'east';
cb.Label.String  = 'v_0 (cm/s)';
cb.Label.FontSize = 8;
cb.FontSize      = 8;

title(tl, sprintf(['Impact kinematics: depth positive into the bed, ' ...
                   't from impact, %g ms pre-impact context'], opt.PreCapMs), ...
      'FontSize', 9);

% ── export ───────────────────────────────────────────────────────────────
written = strings(0,1);
if opt.Save
    if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
    pdfPath = fullfile(opt.OutDir, 'fig1_kinematics.pdf');
    pngPath = fullfile(opt.OutDir, 'fig1_kinematics.png');
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 600);
    written = [string(pdfPath); string(pngPath)];
    fprintf('\n--- written ---\n');
    fprintf('  %s\n', written);
end

if ~opt.Show
    close(fig);
    fig = gobjects(1);
end
fprintf('\n');

R = struct();
R.selected = sel;
R.series   = S;
R.usable   = usable;
R.clim     = CLIM;
R.written  = written;
R.figure   = fig;
end

% ═════════════════════════════════════════════════════════════════════════
function [sel, nGroups] = local_select(K, opt)
%LOCAL_SELECT  One representative trial per (model, height), or everything.
%
%   The representative is the trial whose v0 is closest to the MEAN v0 of its
%   (model, height) group -- the most typical trial by the variable actually
%   plotted, rather than whichever happened to be first on disk.
%
%   Zero-drop rows are thinned to one per model either way: they are identical
%   flat lines by construction, so drawing several would darken a feature that
%   carries no extra information.

    drop = K(~K.isZeroDrop, :);
    nGroups = 0;

    if opt.AllTrials
        keepTags = drop.trialTag;
    else
        g = findgroups(drop.model, drop.dropHeight_mm);
        nGroups = max(g);
        keepTags = strings(nGroups, 1);
        for k = 1:nGroups
            inGrp = find(g == k);
            v0s   = drop.v0_cm_s(inGrp);
            [~, best] = min(abs(v0s - mean(v0s)));
            keepTags(k) = drop.trialTag(inGrp(best));
        end
    end

    rows = ismember(K.trialTag, keepTags) | K.isZeroDrop;
    sel  = K(rows, :);

    % Thin the zero-drop rows to one per model.
    kill = false(height(sel), 1);
    for m = 1:numel(opt.Models)
        z = find(sel.model == opt.Models(m) & sel.isZeroDrop);
        if numel(z) > 1, kill(z(2:end)) = true; end
    end
    sel = sel(~kill, :);
end

function [s, ok] = local_load_series(kinPath, preCapMs)
%LOCAL_LOAD_SERIES  Per-trial series, windowed to [-preCapMs, stopFrame].
%
%   kin.t_s is ALREADY relative to impact (kd_kinematics step 5), so it is used
%   as-is with no further shift. The window keeps every available pre-impact
%   frame up to the cap and ends AT stopFrame: after the stop the traces are
%   either masked or physically meaningless.

    s  = struct('t_ms',[], 'z',[], 'v',[], 'ag',[]);
    ok = false;
    if ~isfile(kinPath), return; end
    try
        L = load(kinPath, 'kin');
    catch
        return
    end
    if ~isfield(L, 'kin'), return; end
    kin = L.kin;
    if ~all(isfield(kin, {'t_s','z','v','a_plus_g','stopFrame'})), return; end

    n   = numel(kin.t_s);
    idx = (1:n).' <= kin.stopFrame & kin.t_s(:) >= -preCapMs/1000;
    if ~any(idx), return; end

    s.t_ms = kin.t_s(idx) * 1000;
    s.z    = kin.z(idx);
    s.v    = kin.v(idx);
    s.ag   = kin.a_plus_g(idx);   % NaN before impact by design; see header (a)
    ok = true;
end

function col = local_v0_color(v0, clim, cmap)
%LOCAL_V0_COLOR  Map one v0 onto the shared colormap.
    f = (v0 - clim(1)) / (clim(2) - clim(1));
    f = min(max(f, 0), 1);
    col = cmap(1 + round(f * (size(cmap,1) - 1)), :);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
