function R = fig_kinematics(varargin)
%FIG_KINEMATICS  Kinematics figure: depth, velocity, net acceleration.
%
%   One panel row per quantity (depth, v, a+g) against time from impact.
%   READ-ONLY: reads the kinematics tree and writes only figure files.
%
%   USAGE
%       fig_kinematics
%       fig_kinematics('Root', 'D:\ME_GRANULAB\JerboaImpact')
%       R = fig_kinematics('Layout', 'grid3x3', 'Save', false);
%       R = fig_kinematics('Style', 'trials');        % QA: every trial
%
%   OPTIONS (name-value)
%       'Root'              results root (parent of 03_RESULTS); '' prompts
%       'Models'            order (default ["Default","Tight","Wide"])
%       'PreCapMs'          pre-impact context cap, ms (default 20)
%       'GridMs'            ensemble time step, ms (default 0.2)
%       'MinReplicates'     ensemble ends below this many trials (default 3)
%       'Average'           'median' (default) | 'mean', across replicates
%       'SmoothMs'          struct('z',1.0,'v',2.0,'ag',3.0): display movmean
%                           window per row, ms, applied AFTER aggregation
%       'IntrusionThreshCm' zero-drop depth-range threshold for an individual
%                           overlay (default 0.1 cm)
%       'Layout'            'per-model' (default): three figures, one per
%                           model, each a 3x1 vertical stack sharing time
%                           'grid3x3': one figure, rows x models, reproducing
%                           the earlier combined layout
%       'Style'             'mean' (default) ensemble curves
%                           'trials' every kept trial individually, for QA
%       'OutDir'            default <Root>/03_RESULTS/_figures
%       'Save'              write PDF + PNG (default true)
%       'Show'              display the figure(s) (default true)
%
%   FIELD NAMES are read from kd_kinematics' save block, not assumed:
%       kin.t_s        time, ALREADY re-zeroed at impact (kd_kinematics step 5:
%                      t_s = t - t(impact_index)), so no further shift is made
%       kin.z          depth, alias of kin.depthRod_cm, zeroed at impact.
%                      Masked to NaN only after stopFrame + postCapFrames, so
%                      the pre-impact segment IS present (unlike kin.a).
%       kin.v          smoothed velocity (the local variable inside
%                      kd_kinematics is v_smooth; it is NOT saved under that
%                      name -- kin.v_smooth returns nothing)
%       kin.a          raw per-frame acceleration, masked to NaN outside
%                      [impact_index, stopFrame]. a+g is RECOMPUTED from this
%                      via src/net_accel.m as g - a; kin.a_plus_g is never
%                      read (stored _kin.mat files written before the 2026-08
%                      sign fix carry the old formula in that column, but the
%                      raw a trace beside it was never affected).
%       kin.impact_index, kin.stopFrame
%   _kin.mat also stores calib alongside kin -- net_accel(kin, calib) takes g
%   from the trial's own calibration rather than assuming 980.
%
%   Table columns (trialTag, model, dropHeight_mm, v0_cm_s, isZeroDrop,
%   kinPath) come from load_kinematics_set.
%
%   SIGN CONVENTION is pipeline-native, NOT mirrored to KD 2007: depth is
%   positive INTO the bed. With the 2026-08 fix, free fall sits at a+g = 0, a
%   resting rod at +g, and deceleration is large positive -- matching KD 2007
%   Fig. 1.
%
%   ══════════════════ DISPLAY LAYER -- NONE OF THIS REACHES A FIT ══════════
%
%   ENSEMBLE (Style 'mean'). For each (model, height) -- h = 0 included, with
%   NO special-case logic -- every trial is interpolated onto one uniform grid
%   relative to impact ('GridMs'), using each trial's OWN [-PreCapMs, t_stop]
%   window (kin.stopFrame read from that trial's own file, never the
%   quarantined table scalar). The pointwise 'Average' (median by default) is
%   taken across replicates at each grid point, using only trials with data
%   there. A curve TERMINATES where fewer than 'MinReplicates' trials remain
%   -- never padded with a shorter trial's final value.
%
%   SMOOTHING is a movmean applied AFTER aggregation, window per row
%   ('SmoothMs'). The pre-smoothing NaN mask (points below MinReplicates
%   coverage) is RE-APPLIED afterwards, so the window blends real neighbours
%   without extending a curve past where the replicates actually end -- and so
%   a plotted line never draws a connector down to the axis at the curve's
%   start or end (MATLAB breaks a line at NaN).
%
%   ZERO-DROP TRIALS are not special-cased in the ensemble: their traces are
%   real, kd_kinematics runs on them like any other trial, and the median
%   represents the majority behaviour without inventing a separate model. A
%   trial is ALSO overlaid individually (depth and v panels only, thin line,
%   the h=0 colour, standard row smoothing) when its depth range over the
%   WHOLE recorded trace exceeds 'IntrusionThreshCm' -- a threshold rule, not
%   a hardcoded tag list. The threshold is evaluated on the full trace because
%   it is a fact about the trial, independent of the display window.
%
%   a+g (row c) is plotted on a LOG y-axis, matching KD 2007 Fig. 1c. Values
%   <= 0 are masked to NaN -- a log axis cannot show them, and clamping would
%   misrepresent rather than omit. h = 0 draws NO curve in this row: the
%   dashed reference line at g is its quasi-static anchor. v is NOT clamped.
%
%   Base MATLAB only.

% ── options ──────────────────────────────────────────────────────────────
opt.Root               = '';
opt.Models             = ["Default","Tight","Wide"];
opt.PreCapMs           = 20;
opt.GridMs             = 0.2;
opt.MinReplicates      = 3;
opt.Average            = 'median';
opt.SmoothMs           = struct('z',1.0, 'v',2.0, 'ag',3.0);
opt.IntrusionThreshCm  = 0.1;
opt.Layout             = 'per-model';
opt.Style              = 'mean';
opt.OutDir             = '';
opt.Save               = true;
opt.Show               = true;
for i = 1:2:numel(varargin), opt.(varargin{i}) = varargin{i+1}; end

opt.Average = lower(char(opt.Average));
if ~ismember(opt.Average, {'median','mean'})
    error('fig_kinematics:badAverage', ...
          'Average must be ''median'' or ''mean'', got "%s".', opt.Average);
end
opt.Layout = lower(char(opt.Layout));
if ~ismember(opt.Layout, {'per-model','grid3x3'})
    error('fig_kinematics:badLayout', ...
          'Layout must be ''per-model'' or ''grid3x3'', got "%s".', opt.Layout);
end
opt.Style = lower(char(opt.Style));
if ~ismember(opt.Style, {'mean','trials'})
    error('fig_kinematics:badStyle', ...
          'Style must be ''mean'' or ''trials'', got "%s".', opt.Style);
end
if ~isstruct(opt.SmoothMs) || ~all(isfield(opt.SmoothMs, {'z','v','ag'}))
    error('fig_kinematics:badSmoothMs', ...
        'SmoothMs must be a struct with fields z, v, ag (ms).');
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
fprintf('  Root              : %s\n', root);
fprintf('  Models            : %s\n', strjoin(cellstr(opt.Models), ', '));
fprintf('  Layout            : %s\n', opt.Layout);
fprintf('  Style             : %s\n', opt.Style);
fprintf('  PreCapMs          : %g ms\n', opt.PreCapMs);
fprintf('  GridMs            : %g ms\n', opt.GridMs);
fprintf('  MinReplicates     : %d\n', opt.MinReplicates);
fprintf('  Average           : %s\n', opt.Average);
fprintf('  SmoothMs          : depth %g ms, v %g ms, a+g %g ms\n', ...
        opt.SmoothMs.z, opt.SmoothMs.v, opt.SmoothMs.ag);
fprintf('  IntrusionThreshCm : %g cm\n', opt.IntrusionThreshCm);
fprintf('  OutDir            : %s\n', opt.OutDir);

% ── data ─────────────────────────────────────────────────────────────────
% Manual exclusions and the zero-drop SCALAR quarantine are applied by the
% loader; every row is a kept trial. The scalar quarantine does not touch the
% traces this figure reads, and windowing here always uses each trial's own
% stored kin.stopFrame, never the table's (quarantined) t_stop_s column.
K = load_kinematics_set(root);
K = K(ismember(K.model, opt.Models), :);
if isempty(K)
    error('fig_kinematics:noTrials', 'No trials for model(s) %s under %s.', ...
          strjoin(cellstr(opt.Models), ', '), root);
end

switch opt.Style
    case 'mean'
        [C, O] = local_build_ensembles(K, opt);
    case 'trials'
        C = local_build_trials(K, opt);
        O = local_empty_overlay();
end
if isempty(C)
    error('fig_kinematics:nothingToPlot', ...
        ['No curve survived. Every group was unreadable, or fell below ' ...
         'MinReplicates = %d across its whole grid.'], opt.MinReplicates);
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
    if ~isempty(O)
        nOut = sum([O.model] == opt.Models(c));
        if nOut > 0
            fprintf(['           %d zero-drop trial(s) overlaid ' ...
                     'individually (depth range > %g cm)\n'], ...
                    nOut, opt.IntrusionThreshCm);
        end
    end
end
fprintf(['\n  display layer: pointwise %ss over replicates; movmean ' ...
         'smoothing per row; a+g on\n  log axis with values <= 0 masked; ' ...
         'quantitative fits use unmodified pipeline kinematics.\n'], ...
        opt.Average);

% ── shared colour scale: v0 on turbo over [0, max] ────────────────────────
allV0 = [C.meanV0];
CLIM  = [0, max(allV0)];
if diff(CLIM) <= 0, CLIM = [0, 1]; end
CMAP  = local_colormap(256);

% ── shared axis limits, computed once across every model and row ─────────
LIMS = local_shared_limits(C, O);

% ── figures ──────────────────────────────────────────────────────────────
figs    = struct();
written = strings(0,1);

if strcmp(opt.Layout, 'per-model')
    for c = 1:nCol
        f = local_draw_one_model(opt.Models(c), C, O, CLIM, CMAP, LIMS, opt);
        figs.(matlab.lang.makeValidName(char(opt.Models(c)))) = f;
        if opt.Save
            stem = sprintf('fig1_kinematics_%s', lower(char(opt.Models(c))));
            written = [written; local_export(f, opt.OutDir, stem)]; %#ok<AGROW>
        end
    end
else % grid3x3
    fg = local_draw_grid3x3(opt.Models, C, O, CLIM, CMAP, LIMS, opt);
    figs.grid3x3 = fg;
    if opt.Save
        written = [written; ...
            local_export(fg, opt.OutDir, 'fig1_kinematics_grid3x3')]; %#ok<AGROW>
    end
end

if opt.Save
    fprintf('\n--- written ---\n');
    fprintf('  %s\n', written);
end

if ~opt.Show
    fn = fieldnames(figs);
    for i = 1:numel(fn), close(figs.(fn{i})); figs.(fn{i}) = gobjects(1); end
end
fprintf('\n');

R = struct();
R.curves   = C;
R.overlays = O;
R.clim     = CLIM;
R.limits   = LIMS;
R.written  = written;
R.figures  = figs;
end

% ═════════════════════════════════════════════════════════════════════════
%  DATA: loading and windowing
% ═════════════════════════════════════════════════════════════════════════
function s = local_load_kin(kinPath)
%LOCAL_LOAD_KIN  One trial's kin + calib, with a+g recomputed via net_accel.
%   Returns [] if the file is missing, unreadable, or missing a required
%   field. Everything downstream slices from here, so the file is opened once
%   per trial regardless of how many things are computed from it.
    s = [];
    if ~isfile(kinPath), return; end
    try
        L = load(kinPath, 'kin', 'calib');
    catch
        return
    end
    if ~isfield(L, 'kin') || ~all(isfield(L.kin, ...
            {'t_s','z','v','a','impact_index','stopFrame'}))
        return
    end
    kin = L.kin;
    calib = [];
    if isfield(L, 'calib'), calib = L.calib; end

    s = struct();
    s.t_ms  = kin.t_s(:) * 1000;
    s.z     = kin.z(:);
    s.v     = kin.v(:);
    s.ag    = net_accel(kin, calib);           % g - a, from the raw trace
    s.stopFrame = kin.stopFrame;
    n = numel(s.t_ms);
    if s.stopFrame >= 1 && s.stopFrame <= n
        s.stop_ms = s.t_ms(s.stopFrame);
    else
        s.stop_ms = NaN;
    end
end

function w = local_window(s, preCapMs)
%LOCAL_WINDOW  [-preCapMs, stopFrame] slice of a loaded trace, as its own
%   struct with the same field names, ready to interpolate or plot directly.
    n = numel(s.t_ms);
    idx = (1:n).' <= s.stopFrame & s.t_ms >= -preCapMs;
    w = struct('t_ms', s.t_ms(idx), 'z', s.z(idx), 'v', s.v(idx), ...
               'ag', s.ag(idx), 'stop_ms', s.stop_ms);
end

% ═════════════════════════════════════════════════════════════════════════
%  ENSEMBLE CONSTRUCTION (Style 'mean')
% ═════════════════════════════════════════════════════════════════════════
function [C, O] = local_build_ensembles(K, opt)
%LOCAL_BUILD_ENSEMBLES  One aggregate curve per (model, height); h = 0 goes
%   through the identical path. O collects the individually-overlaid
%   high-intrusion zero-drop trials (depth/v only).

    C = local_empty_curve();
    O = local_empty_overlay();

    g = findgroups(K.model, K.dropHeight_mm);
    for k = 1:max(g)
        rows = find(g == k);
        isZD = all(K.isZeroDrop(rows));

        loaded = cell(numel(rows),1);
        for j = 1:numel(rows)
            loaded{j} = local_load_kin(K.kinPath(rows(j)));
        end
        okLoad = ~cellfun(@isempty, loaded);
        loaded = loaded(okLoad);
        rowsOk = rows(okLoad);
        if isempty(loaded), continue; end

        W = cellfun(@(s) local_window(s, opt.PreCapMs), loaded, ...
                    'UniformOutput', false);

        tEnd = max(cellfun(@(w) local_last_or(w.t_ms, -opt.PreCapMs), W));
        if tEnd <= -opt.PreCapMs, continue; end
        grid = (-opt.PreCapMs : opt.GridMs : tEnd).';

        [zRaw, nz] = local_aggregate_on_grid(W, 'z',  grid, opt.MinReplicates, opt.Average);
        [vRaw, nv] = local_aggregate_on_grid(W, 'v',  grid, opt.MinReplicates, opt.Average);
        [aRaw, na] = local_aggregate_on_grid(W, 'ag', grid, opt.MinReplicates, opt.Average);

        medStopMs = median(cellfun(@(w) w.stop_ms, W), 'omitnan');

        zC = local_smooth_row(zRaw, opt.GridMs, opt.SmoothMs.z,  medStopMs, 'depth');
        vC = local_smooth_row(vRaw, opt.GridMs, opt.SmoothMs.v,  medStopMs, 'v');
        aC = local_smooth_row(aRaw, opt.GridMs, opt.SmoothMs.ag, medStopMs, 'a+g');

        v0s  = K.v0_cm_s(rowsOk);
        nAll = [nz; nv; na];
        C(end+1) = struct( ...                                    %#ok<AGROW>
            'model',        K.model(rowsOk(1)), ...
            'height',       K.dropHeight_mm(rowsOk(1)), ...
            'meanV0',       mean(v0s, 'omitnan'), ...
            'isZeroDrop',   isZD, ...
            'nMin',         min(nAll(nAll > 0), [], 'omitnan'), ...
            'nMax',         numel(W), ...
            't_ms', grid, 'z', zC, 'v', vC, 'ag', aC);

        % ── individual zero-drop outliers (depth/v only) ──────────────────
        if isZD
            for j = 1:numel(loaded)
                sFull = loaded{j};
                depthRange = max(sFull.z, [], 'omitnan') - min(sFull.z, [], 'omitnan');
                if ~(depthRange > opt.IntrusionThreshCm), continue; end

                w = W{j};
                zO = local_smooth_row(w.z, [], opt.SmoothMs.z, medStopMs, ...
                                       'depth (overlay)', w.t_ms);
                vO = local_smooth_row(w.v, [], opt.SmoothMs.v, medStopMs, ...
                                       'v (overlay)', w.t_ms);
                O(end+1) = struct('model', K.model(rowsOk(j)), ...
                    'meanV0', 0, 't_ms', w.t_ms, 'z', zO, 'v', vO); %#ok<AGROW>
            end
        end
    end
end

function C = local_build_trials(K, opt)
%LOCAL_BUILD_TRIALS  Every kept trial as its OWN curve, unaveraged. QA view.
%   Same struct shape as the ensemble curves, so the drawing code is shared:
%   nMin = nMax = 1, and each trial is smoothed individually with the same
%   row windows so the two styles are visually comparable.
    C = local_empty_curve();
    for i = 1:height(K)
        s = local_load_kin(K.kinPath(i));
        if isempty(s), continue; end
        w = local_window(s, opt.PreCapMs);
        if isempty(w.t_ms), continue; end

        zC = local_smooth_row(w.z,  [], opt.SmoothMs.z,  w.stop_ms, 'depth', w.t_ms);
        vC = local_smooth_row(w.v,  [], opt.SmoothMs.v,  w.stop_ms, 'v',     w.t_ms);
        aC = local_smooth_row(w.ag, [], opt.SmoothMs.ag, w.stop_ms, 'a+g',   w.t_ms);

        C(end+1) = struct( ...                                    %#ok<AGROW>
            'model',      K.model(i), ...
            'height',     K.dropHeight_mm(i), ...
            'meanV0',     K.v0_cm_s(i), ...
            'isZeroDrop', K.isZeroDrop(i), ...
            'nMin', 1, 'nMax', 1, ...
            't_ms', w.t_ms, 'z', zC, 'v', vC, 'ag', aC);
    end
end

function C = local_empty_curve()
    C = struct('model',{}, 'height',{}, 'meanV0',{}, 'isZeroDrop',{}, ...
               'nMin',{}, 'nMax',{}, 't_ms',{}, 'z',{}, 'v',{}, 'ag',{});
end

function O = local_empty_overlay()
    O = struct('model',{}, 'meanV0',{}, 't_ms',{}, 'z',{}, 'v',{});
end

function v = local_last_or(x, dflt)
    if isempty(x), v = dflt; else, v = x(end); end
end

% ═════════════════════════════════════════════════════════════════════════
%  AGGREGATION AND SMOOTHING
% ═════════════════════════════════════════════════════════════════════════
function [m, n] = local_aggregate_on_grid(W, field, grid, minRep, avg)
%LOCAL_AGGREGATE_ON_GRID  Pointwise median/mean at each grid point, with a
%   coverage floor.
%
%   Each replicate is interpolated onto the grid with NaN outside its own
%   window, so a trial contributes only where it actually has data. Points
%   covered by fewer than minRep replicates are NaN: the curve terminates
%   rather than trailing into a one- or two-trial extrapolation, and no trial
%   is ever padded with its final value.

    Y = nan(numel(grid), numel(W));
    for j = 1:numel(W)
        t = W{j}.t_ms;
        y = W{j}.(field);
        ok = isfinite(t) & isfinite(y);
        if nnz(ok) < 2, continue; end
        Y(:,j) = interp1(t(ok), y(ok), grid, 'linear', NaN);
    end
    n = sum(~isnan(Y), 2);
    if strcmp(avg, 'median')
        m = median(Y, 2, 'omitnan');
    else
        m = mean(Y, 2, 'omitnan');
    end
    m(n < minRep) = NaN;
end

function y = local_smooth_row(yRaw, gridMs, windowMs, medianStopMs, rowLabel, tMs)
%LOCAL_SMOOTH_ROW  movmean AFTER aggregation, then re-mask to the original
%   support so smoothing cannot extend a curve past where the data ends.
%
%   gridMs is the uniform sample spacing when known (ensemble curves); pass []
%   and a t_ms vector instead for a trial's own (also uniform, but different)
%   sampling, and the spacing is taken from the median step of tMs.
    if isempty(gridMs)
        if numel(tMs) < 2
            y = yRaw; return
        end
        gridMs = median(diff(tMs), 'omitnan');
    end
    if ~(gridMs > 0) || ~isfinite(gridMs)
        y = yRaw; return
    end

    w = max(1, round(windowMs / gridMs));
    y = movmean(yRaw, w, 'omitnan');
    y(isnan(yRaw)) = NaN;

    if isfinite(medianStopMs) && medianStopMs > 0 && windowMs > 0.10*medianStopMs
        warning('fig_kinematics:overSmooth', ...
            ['%s smoothing window (%.2f ms) exceeds 10%% of the median stop ' ...
             'time (%.2f ms). The curve may be over-smoothed.'], ...
            rowLabel, windowMs, medianStopMs);
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  SHARED LIMITS AND COLOUR
% ═════════════════════════════════════════════════════════════════════════
function LIMS = local_shared_limits(C, O)
%LOCAL_SHARED_LIMITS  One set of x/y limits per row, shared across every
%   model and (for per-model layout) every figure -- computed once so the
%   panels stay comparable regardless of which model is drawn where.

    allT = [vertcat(C.t_ms); vertcat(O.t_ms)];
    XL = [min(allT), max(allT)];

    zAll = [vertcat(C.z); vertcat(O.z)];
    YL_z = local_pad_ylim(zAll);

    vAll = [vertcat(C.v); vertcat(O.v)];
    YL_v = local_pad_ylim(vAll);

    % a+g: only h > 0 curves carry this row.
    isZD = [C.isZeroDrop];
    agAll = vertcat(C(~isZD).ag);
    YL_ag = local_log_ylim(agAll);

    LIMS = struct('t', XL, 'z', YL_z, 'v', YL_v, 'ag_log', YL_ag);
end

function yl = local_pad_ylim(vals)
    vals = vals(isfinite(vals));
    if isempty(vals), yl = [-1 1]; return; end
    lo = min(vals); hi = max(vals); r = hi - lo;
    if r == 0, r = max(abs(hi), 1); end
    yl = [lo - 0.05*r, hi + 0.05*r];
end

function yl = local_log_ylim(vals)
%LOCAL_LOG_YLIM  Data-driven bounds that always cover the KD-like reference
%   range [1e2, 3e4] (matching Fig. 1c), extended further only if the actual
%   data exceeds it.
    vals = vals(isfinite(vals) & vals > 0);
    lo = 1e2; hi = 3e4;
    if ~isempty(vals)
        lo = min(lo, 10^floor(log10(min(vals))));
        hi = max(hi, 10^ceil(log10(max(vals))));
    end
    yl = [lo, hi];
end

function cmap = local_colormap(nLevels)
%LOCAL_COLORMAP  Perceptually ordered rainbow; turbo, or parula where absent.
    if exist('turbo', 'file') == 2 || exist('turbo', 'builtin') == 5
        cmap = turbo(nLevels);
    else
        warning('fig_kinematics:noTurbo', ...
            'turbo() unavailable; falling back to parula.');
        cmap = parula(nLevels);
    end
end

function col = local_v0_color(v0, clim, cmap)
    f = (v0 - clim(1)) / (clim(2) - clim(1));
    if ~isfinite(f), f = 0; end
    f = min(max(f, 0), 1);
    col = cmap(1 + round(f * (size(cmap,1) - 1)), :);
end

% ═════════════════════════════════════════════════════════════════════════
%  DRAWING
% ═════════════════════════════════════════════════════════════════════════
function fig = local_draw_one_model(model, C, O, CLIM, CMAP, LIMS, opt)
%LOCAL_DRAW_ONE_MODEL  3x1 vertical stack, single-column APS width.
    fig = figure('Color','w','Units','inches','Position',[1 1 3.4 7.0], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, 3, 1, 'Padding','compact', 'TileSpacing','compact');

    ax1 = nexttile(tl); local_draw_panel(ax1, 'z',  model, C, O, CLIM, CMAP, LIMS);
    ax2 = nexttile(tl); local_draw_panel(ax2, 'v',  model, C, O, CLIM, CMAP, LIMS);
    ax3 = nexttile(tl); local_draw_panel(ax3, 'ag', model, C, O, CLIM, CMAP, LIMS);

    xlabel(ax3, 't - t_impact (ms)', 'FontSize', 8, 'Interpreter', 'none');
    set([ax1 ax2], 'XTickLabel', []);
end

function fig = local_draw_grid3x3(models, C, O, CLIM, CMAP, LIMS, opt)
%LOCAL_DRAW_GRID3X3  One figure, rows x models -- the earlier combined layout.
    nCol = numel(models);
    fig = figure('Color','w','Units','inches','Position',[1 1 7.0 6.4], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl  = tiledlayout(fig, 3, nCol, 'Padding','compact', 'TileSpacing','compact');

    ROWS = {'z','v','ag'};
    AX = gobjects(3, nCol);
    for r = 1:3
        for c = 1:nCol
            ax = nexttile(tl, (r-1)*nCol + c);
            local_draw_panel(ax, ROWS{r}, models(c), C, O, CLIM, CMAP, LIMS);
            AX(r,c) = ax;
            if r < 3, set(ax, 'XTickLabel', []); end
            if c > 1, set(ax, 'YTickLabel', []); end
            if r == 3 && c == 1
                xlabel(ax, 't - t_impact (ms)', 'FontSize', 8, 'Interpreter', 'none');
            end
        end
    end
end

function local_draw_panel(ax, rowKey, model, C, O, CLIM, CMAP, LIMS)
%LOCAL_DRAW_PANEL  One (row, model) panel: depth, v, or a+g -- shared by both
%   layouts and both styles, since C/O are already in a common shape.

    hold(ax, 'on');

    % t = 0 reference, drawn FIRST so it sits behind the data.
    xline(ax, 0, '--', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.75, ...
          'HandleVisibility', 'off');
    if ismember(rowKey, {'z','v'})
        yline(ax, 0, '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end

    mkMap = containers.Map({'Default','Tight','Wide'}, {'o','s','^'});
    % Deliberately different from depth_scaling.m's Tight=^/Wide=s convention:
    % this figure's marker mapping is specified independently in its own
    % brief (Default circle, Tight square, Wide triangle) and is not meant to
    % be unified with that script's.

    inCol = find([C.model] == model);
    for k = inCol
        if strcmp(rowKey, 'ag') && C(k).isZeroDrop
            continue      % h = 0 has no a+g curve; the g-line is its anchor
        end
        y = C(k).(rowKey);
        if isempty(y) || all(isnan(y)), continue; end
        if strcmp(rowKey, 'ag')
            y(y <= 0) = NaN;   % log axis: mask non-positive values, never clamp
        end
        col = local_v0_color(C(k).meanV0, CLIM, CMAP);
        plot(ax, C(k).t_ms, y, '-', 'LineWidth', 1.5, 'Color', col);

        if strcmp(rowKey, 'z')
            last = find(~isnan(y), 1, 'last');
            if ~isempty(last)
                mk = 'd';
                if isKey(mkMap, char(model)), mk = mkMap(char(model)); end
                plot(ax, C(k).t_ms(last), y(last), mk, 'MarkerSize', 5, ...
                     'MarkerFaceColor', col, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
                     'LineWidth', 0.5);
            end
        end
    end

    if ismember(rowKey, {'z','v'})
        inColO = find([O.model] == model);
        for k = inColO
            y = O(k).(rowKey);
            if isempty(y) || all(isnan(y)), continue; end
            col = local_v0_color(O(k).meanV0, CLIM, CMAP);   % meanV0 = 0 -> h=0 colour
            plot(ax, O(k).t_ms, y, '-', 'LineWidth', 0.75, 'Color', col);
        end
    end

    if strcmp(rowKey, 'ag')
        set(ax, 'YScale', 'log');
        yline(ax, 980, '--', 'Color', 'k', 'LineWidth', 1.0, ...
              'HandleVisibility', 'off');   % KD's dashed line is g
        ylim(ax, LIMS.ag_log);
        ylabel(ax, 'a + g (cm/s^2)', 'FontSize', 8, 'Interpreter', 'none');
    elseif strcmp(rowKey, 'v')
        ylim(ax, LIMS.v);
        ylabel(ax, 'v (cm/s)', 'FontSize', 8, 'Interpreter', 'none');
    else
        ylim(ax, LIMS.z);
        ylabel(ax, 'depth (cm)', 'FontSize', 8, 'Interpreter', 'none');
    end
    xlim(ax, LIMS.t);

    box(ax, 'on');
    grid(ax, 'on');
    set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
            'FontSize', 8, 'LineWidth', 0.5, 'Layer', 'top');
end

% ═════════════════════════════════════════════════════════════════════════
function paths = local_export(fig, outDir, stem)
    if ~isfolder(outDir), mkdir(outDir); end
    pdfPath = fullfile(outDir, [stem '.pdf']);
    pngPath = fullfile(outDir, [stem '.png']);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 600);
    paths = [string(pdfPath); string(pngPath)];
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
