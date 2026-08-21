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
%       R = fig_kinematics('Style', 'trials');        % QA: every trial
%
%   OPTIONS (name-value)
%       'Root'           results root (parent of 03_RESULTS); '' prompts
%       'Models'         column order (default ["Default","Tight","Wide"])
%       'Style'          'mean' (default) one ensemble curve per height
%                        'trials'         every kept trial, for QA
%       'PreCapMs'       pre-impact context cap, ms (default 20)
%       'GridMs'         ensemble time step, ms (default 0.2)
%       'SmoothAccelMs'  display movmean window on a+g, ms (default 1.5)
%       'MinReplicates'  ensemble ends below this many trials (default 3)
%       'OutDir'         default <Root>/03_RESULTS/_figures
%       'Save'           write PDF + PNG (default true)
%       'Show'           display the figure (default true)
%
%   FIELD NAMES are read from kd_kinematics' save block, not assumed:
%       kin.t_s        time, ALREADY re-zeroed at impact (kd_kinematics step 5:
%                      t_s = t - t(impact_index)), so no further shift is made
%       kin.z          depth, alias of kin.depthRod_cm, zeroed at impact
%       kin.v          smoothed velocity. The local variable inside
%                      kd_kinematics is v_smooth; it is NOT saved under that
%                      name, and reading kin.v_smooth returns nothing.
%       kin.stopFrame
%   Table columns (trialTag, model, dropHeight_mm, v0_cm_s, isZeroDrop,
%   kinPath) come from load_kinematics_set.
%
%   ENSEMBLE MEANS (Style 'mean'). For each (model, height) every kept trial is
%   interpolated onto one uniform grid relative to impact, then averaged across
%   replicates point by point, using only the trials that still have data
%   there. Pre-impact context and stop times differ per trial, so coverage
%   falls off at both ends; a curve TERMINATES where fewer than MinReplicates
%   trials remain. Shorter trials are never padded with their final value,
%   which would invent a plateau that no trial measured.
%
%   h = 0 ANCHOR. Zero-drop trials go through the IDENTICAL pipeline as every
%   other height for depth and velocity -- same grid, same averaging, same
%   termination rule -- from their measured traces. load_kinematics_set
%   quarantines the zero-drop SCALARS (d_final, t_stop, a_stop are NaN), but
%   the kin TRACES are intact and are what this figure uses. Those trials have
%   no meaningful stop, so their window is the full available span rather than
%   [-PreCapMs, t_stop].
%
%   ROW 3 uses net_accel(kin), which computes the net grain acceleration
%   as g - a from the stored raw acceleration trace. In this convention,
%   free fall is 0, rest is +g, and deceleration is positive.

%   SIGN CONVENTION is pipeline-native and is NOT mirrored to KD 2007: depth is
%   positive INTO the bed (rod_displacement, signConvention +1).
%
%   ROW 3 uses net_accel(kin), which computes the net grain acceleration
%   as g - a from the stored raw acceleration trace.
%
%   This is implemented exactly as specified and flagged rather than silently
%   reconciled, because the fix belongs in kd_kinematics or in the caption, not
%   in a plotting script. Confirm against a known trial before this figure goes
%   in a paper.
%
%   DISPLAY LAYER. Everything below affects the drawing only. Nothing is
%   written back, and no fit ever sees it:
%     - v and a+g are clamped at max(value, 0). Depth is unclamped.
%     - a+g gets a light movmean (SmoothAccelMs) AFTER averaging and BEFORE
%       clamping. Depth and v get none: replicate averaging is enough.
%     - no raw markers, no scatter, no shaded bands.

% ── options ──────────────────────────────────────────────────────────────
opt.Root          = '';
opt.Models        = ["Default","Tight","Wide"];
opt.Style         = 'mean';
opt.PreCapMs      = 20;
opt.GridMs        = 0.2;
opt.SmoothAccelMs = 1.5;
opt.MinReplicates = 3;
opt.OutDir        = '';
opt.Save          = true;
opt.Show          = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

opt.Style  = lower(char(opt.Style));
if ~ismember(opt.Style, {'mean','trials'})
    error('fig_kinematics:badStyle', ...
          'Style must be ''mean'' or ''trials'', got "%s".', opt.Style);
end
opt.Models = string(opt.Models);
nCol = numel(opt.Models);

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(thisDir), 'src'));

root = resolve_output_root(opt.Root);
if isempty(root), fprintf('Cancelled.\n'); return; end
if isempty(opt.OutDir)
    opt.OutDir = fullfile(root, '03_RESULTS', '_figures');
end

fprintf('\n=== fig_kinematics ===\n');
fprintf('  Root          : %s\n', root);
fprintf('  Models        : %s\n', strjoin(cellstr(opt.Models), ', '));
fprintf('  Style         : %s\n', opt.Style);
fprintf('  PreCapMs      : %g ms\n', opt.PreCapMs);
fprintf('  GridMs        : %g ms\n', opt.GridMs);
fprintf('  SmoothAccelMs : %g ms (a+g display smooth)\n', opt.SmoothAccelMs);
fprintf('  MinReplicates : %d\n', opt.MinReplicates);
fprintf('  OutDir        : %s\n', opt.OutDir);

% ── data ─────────────────────────────────────────────────────────────────
% Manual exclusions and the zero-drop scalar quarantine are applied by the
% loader; every row returned is a kept trial.
K = load_kinematics_set(root);
K = K(ismember(K.model, opt.Models), :);
if isempty(K)
    error('fig_kinematics:noTrials', 'No trials for model(s) %s under %s.', ...
          strjoin(cellstr(opt.Models), ', '), root);
end

switch opt.Style
    case 'mean',   C = local_build_ensembles(K, opt);
    case 'trials', C = local_build_trials(K, opt);
end
if isempty(C)
    error('fig_kinematics:nothingToPlot', ...
        ['No curve survived. Every group was either unreadable or fell below ' ...
         'MinReplicates = %d across the whole grid.'], opt.MinReplicates);
end

% ── provenance ───────────────────────────────────────────────────────────
fprintf('\n--- curves per panel ---\n');
for c = 1:nCol
    inCol = [C.model] == opt.Models(c);
    if ~any(inCol)
        fprintf('  %-8s none\n', opt.Models(c));
        continue
    end
    reps = [C(inCol).nMax];
    fprintf('  %-8s %2d height curve(s), replicates %d-%d per height\n', ...
            opt.Models(c), sum(inCol), min(reps), max(reps));
    z0 = C(inCol);
    z0 = z0([z0.isZeroDrop]);
    if isempty(z0)
        fprintf('           (no h = 0 anchor: no zero-drop trials for this model)\n');
    elseif all(arrayfun(@(s) all(isnan(s.z)), z0))
        fprintf(['           WARNING: h = 0 anchor is empty -- fewer than %d ' ...
                 'zero-drop replicates\n'], opt.MinReplicates);
    end
end
fprintf(['\n  display layer: per-height ensemble means; a+g movmean %g ms; ' ...
         'v and a+g\n  clamped at 0; quantitative fits use unmodified ' ...
         'pipeline kinematics.\n'], opt.SmoothAccelMs);

% ── colour: mean v0 per height, on turbo over [0, max] ───────────────────
CLIM = [0, max([C.meanV0])];
if diff(CLIM) <= 0, CLIM = [0, 1]; end
CMAP = local_colormap(256);

% ── figure ───────────────────────────────────────────────────────────────
fig = figure('Color','w','Units','inches','Position',[1 1 7.0 6.4], ...
             'Visible', local_tern(opt.Show,'on','off'));
tl  = tiledlayout(fig, 3, nCol, 'Padding','compact', 'TileSpacing','compact');

ROWS = { 'z',  'depth (cm)'     ; ...
         'v',  'v (cm/s)'       ; ...
         'ag', 'a + g (cm/s^2)' };

AX = gobjects(3, nCol);
for r = 1:3
    for c = 1:nCol
        ax = nexttile(tl, (r-1)*nCol + c);
        hold(ax,'on'); box(ax,'on');
        AX(r,c) = ax;

        % Thin light zero-line on depth and velocity.
        if r <= 2
            yline(ax, 0, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
        end

        for k = find([C.model] == opt.Models(c))
            y = C(k).(ROWS{r,1});
            if isempty(y) || all(isnan(y)), continue; end
            plot(ax, C(k).t_ms, y, '-', 'LineWidth', 0.9, ...
                 'Color', local_v0_color(C(k).meanV0, CLIM, CMAP));

            % 'h = 0' label at the right end of the anchor, depth and a+g only.
            % Omitted in v, where it would collide with the axis.
            if C(k).isZeroDrop && ismember(r, [1 3])
                last = find(~isnan(y), 1, 'last');
                if ~isempty(last)
                    text(ax, C(k).t_ms(last), y(last), ' h = 0', ...
                         'FontSize', 7, 'VerticalAlignment', 'middle', ...
                         'HorizontalAlignment', 'left');
                end
            end
        end

        set(ax, 'FontSize', 8, 'LineWidth', 0.5, 'Layer', 'top');
        if r == 1
            title(ax, opt.Models(c), 'FontSize', 9, 'Interpreter', 'none');
        end
        if r == 3
            xlabel(ax, 't - t_impact (ms)', 'FontSize', 8, 'Interpreter', 'none');
        end
        if c == 1
            ylabel(ax, ROWS{r,2}, 'FontSize', 8, 'Interpreter', 'none');
        end
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
R.curves  = C;
R.clim    = CLIM;
R.written = written;
R.figure  = fig;
end

% ═════════════════════════════════════════════════════════════════════════
function C = local_build_ensembles(K, opt)
%LOCAL_BUILD_ENSEMBLES  One mean curve per (model, height).
%
%   Every height, h = 0 included, goes through this same path: interpolate each
%   replicate onto a shared grid, average where the replicates actually have
%   data, and stop where coverage drops below MinReplicates.

    G_CM_S2 = 980;
    C = struct('model',{}, 'height',{}, 'meanV0',{}, 'isZeroDrop',{}, ...
               'nMin',{}, 'nMax',{}, 't_ms',{}, 'z',{}, 'v',{}, 'ag',{});

    g = findgroups(K.model, K.dropHeight_mm);
    for k = 1:max(g)
        rows = find(g == k);
        isZD = all(K.isZeroDrop(rows));

        % Load every replicate's windowed trace.
        T = cell(numel(rows),1);
        for j = 1:numel(rows)
            T{j} = local_load_trace(K.kinPath(rows(j)), opt.PreCapMs, isZD);
        end
        T = T(~cellfun(@isempty, T));
        if isempty(T), continue; end

        % Shared grid: from the pre-impact cap out to the longest replicate.
        tEnd = max(cellfun(@(s) s.t_ms(end), T));
        if tEnd <= -opt.PreCapMs, continue; end
        grid = (-opt.PreCapMs : opt.GridMs : tEnd).';

        [zm, nz] = local_mean_on_grid(T, 'z',  grid, opt.MinReplicates);
        [vm, nv] = local_mean_on_grid(T, 'v',  grid, opt.MinReplicates);
        [am, na] = local_mean_on_grid(T, 'ag', grid, opt.MinReplicates);

        % ── display layer ────────────────────────────────────────────────
        % a+g: light smooth AFTER averaging, BEFORE clamping, so the clamp acts
        % on the curve that will actually be drawn.
        w = max(1, round(opt.SmoothAccelMs / opt.GridMs));
        am = movmean(am, w, 'omitnan');

        if isZD
            % a+g is NaN outside the [impact, stop] window, and for a zero-drop
            % trial that window is degenerate -- there is no impact and no
            % stop -- so the pipeline has nothing to report here.
            %
            % DERIVATION (as specified): the measured v is flat, so
            %     a = dv/dt = 0   ->   a + g = g
            % and the anchor is the constant g over the span the depth curve
            % covers.

            am = nan(size(grid));
            am(~isnan(zm)) = G_CM_S2;
            na = nz;
        end

        % v and a+g clamped at zero for display; depth left alone.
        vm = max(vm, 0);
        am = max(am, 0);

        v0s = K.v0_cm_s(rows);
        nAll = [nz; nv; na];
        C(end+1) = struct( ...                                    %#ok<AGROW>
            'model',      K.model(rows(1)), ...
            'height',     K.dropHeight_mm(rows(1)), ...
            'meanV0',     mean(v0s, 'omitnan'), ...
            'isZeroDrop', isZD, ...
            'nMin',       min(nAll(nAll > 0), [], 'omitnan'), ...
            'nMax',       numel(T), ...
            't_ms',       grid, 'z', zm, 'v', vm, 'ag', am);
    end
end

function C = local_build_trials(K, opt)
%LOCAL_BUILD_TRIALS  Every kept trial, unaveraged. QA view.
%   Same display clamps as the ensemble path so the two are comparable, but no
%   averaging and no smoothing: this exists to show what the means came from.
    C = struct('model',{}, 'height',{}, 'meanV0',{}, 'isZeroDrop',{}, ...
               'nMin',{}, 'nMax',{}, 't_ms',{}, 'z',{}, 'v',{}, 'ag',{});
    for i = 1:height(K)
        isZD = K.isZeroDrop(i);
        s = local_load_trace(K.kinPath(i), opt.PreCapMs, isZD);
        if isempty(s), continue; end
        C(end+1) = struct( ...                                    %#ok<AGROW>
            'model',      K.model(i), ...
            'height',     K.dropHeight_mm(i), ...
            'meanV0',     K.v0_cm_s(i), ...
            'isZeroDrop', isZD, ...
            'nMin',       1, 'nMax', 1, ...
            't_ms',       s.t_ms, 'z', s.z, ...
            'v',          max(s.v, 0), 'ag', max(s.ag, 0));
    end
end

function s = local_load_trace(kinPath, preCapMs, isZeroDrop)
%LOCAL_LOAD_TRACE  One trial's windowed series, or [] if unusable.
%
%   kin.t_s is ALREADY relative to impact (kd_kinematics step 5), so it is used
%   as-is. The window keeps every available pre-impact frame up to the cap and
%   ends AT stopFrame -- after the stop the traces are masked or meaningless.
%
%   Zero-drop trials have no meaningful stop, so their window is the full
%   available span instead.

    s = [];
    if ~isfile(kinPath), return; end
    try
        L = load(kinPath, 'kin');
    catch
        return
    end
    if ~isfield(L, 'kin'), return; end
    kin = L.kin;
    if ~all(isfield(kin, {'t_s','z','v','a','stopFrame'})), return; end

    n = numel(kin.t_s);
    if isZeroDrop
        idx = true(n,1);
    else
        idx = (1:n).' <= kin.stopFrame;
    end
    idx = idx & kin.t_s(:) >= -preCapMs/1000;
    if ~any(idx), return; end

    s = struct();
    s.t_ms = kin.t_s(idx) * 1000;
    s.z    = kin.z(idx);
    s.v    = kin.v(idx);
    ag     = net_accel(kin);
    s.ag   = ag(idx);
end

function [m, n] = local_mean_on_grid(T, field, grid, minRep)
%LOCAL_MEAN_ON_GRID  Replicate mean at each grid point, with a coverage floor.
%
%   Each replicate is interpolated onto the grid with NaN outside its own span,
%   so a trial contributes only where it actually has data. Points covered by
%   fewer than minRep replicates are NaN: the curve terminates rather than
%   trailing off into a one-trial extrapolation, and no trial is ever padded
%   with its final value.

    Y = nan(numel(grid), numel(T));
    for j = 1:numel(T)
        t = T{j}.t_ms;
        y = T{j}.(field);
        ok = isfinite(t) & isfinite(y);
        if nnz(ok) < 2, continue; end
        % 'linear' with NaN fill: outside a replicate's span it contributes
        % nothing, which is what makes the coverage count meaningful.
        Y(:,j) = interp1(t(ok), y(ok), grid, 'linear', NaN);
    end
    n = sum(~isnan(Y), 2);
    m = mean(Y, 2, 'omitnan');
    m(n < minRep) = NaN;
end

function cmap = local_colormap(nLevels)
%LOCAL_COLORMAP  Perceptually ordered rainbow; turbo, or parula where absent.
    if exist('turbo', 'file') == 2 || exist('turbo', 'builtin') == 5
        cmap = turbo(nLevels);
    else
        % turbo is R2020b+. parula keeps the figure readable on older MATLAB
        % rather than erroring on a colormap name.
        warning('fig_kinematics:noTurbo', ...
            'turbo() unavailable; falling back to parula.');
        cmap = parula(nLevels);
    end
end

function col = local_v0_color(v0, clim, cmap)
%LOCAL_V0_COLOR  Map one v0 onto the shared colormap.
    f = (v0 - clim(1)) / (clim(2) - clim(1));
    if ~isfinite(f), f = 0; end
    f = min(max(f, 0), 1);
    col = cmap(1 + round(f * (size(cmap,1) - 1)), :);
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
