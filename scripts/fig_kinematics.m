function R = fig_kinematics(varargin)
%FIG_KINEMATICS  Kinematics figure: depth, velocity, net acceleration.
%
%   One panel row per quantity (depth, v, a+g) against time from impact.
%   READ-ONLY: reads the kinematics tree and writes figure files plus one
%   small CSV of per-height reportable scalars.
%
%   USAGE
%       fig_kinematics
%       fig_kinematics('Root', 'D:\ME_GRANULAB\JerboaImpact')
%       fig_kinematics('Diagnose', true)              % print the diagnostics
%       R = fig_kinematics('Layout', 'grid3x3', 'Save', false);
%       R = fig_kinematics('Style', 'trials');        % QA: every trial
%
%   OPTIONS (name-value)
%       'Root'              results root (parent of 03_RESULTS); '' prompts
%       'Models'            order (default ["Default","Tight","Wide"])
%       'PreCapMs'          pre-impact context cap, ms (default 20)
%       'GridMs'            ensemble time step, ms (default 0.2)
%       'MinReplicates'     minimum trials for a height to be PLOTTED AT ALL
%                           (default 3). It no longer decides where a curve
%                           ends -- see ENSEMBLE below.
%       'Average'           'median' (default) | 'mean', across replicates
%       'SmoothMs'          struct('z',1.0,'v',2.0,'ag',1.5): display movmean
%                           window per row, ms, applied AFTER aggregation
%       'AccelYMin'         a+g log-axis floor, cm/s^2 (default 1e2). KD used
%                           1e1; ours is set by this rig's noise floor.
%       'PreImpactAccelMs'  window for the display-only pre-impact
%                           acceleration fit, ms (default 1.0)
%       'ShowZeroDrop'      include zero-drop trials (default FALSE -- they
%                           are reserved for scalar measurements)
%       'IntrusionThreshCm' zero-drop depth-range threshold for an individual
%                           overlay (default 0.1 cm; only used when
%                           ShowZeroDrop is true)
%       'Layout'            'per-model' (default) | 'grid3x3'
%       'Style'             'mean' (default) | 'trials' (QA)
%       'Diagnose'          print per-height diagnostics (default false)
%       'OutDir'            default <Root>/03_RESULTS/_figures
%       'Save'              write PDF + PNG + CSV (default true)
%       'Show'              display the figure(s) (default true)
%
%   FIELD NAMES are read from kd_kinematics' save block, not assumed:
%       kin.t_s        time, ALREADY re-zeroed at impact (kd_kinematics step 5)
%       kin.z          depth (alias of kin.depthRod_cm), zeroed at impact
%       kin.v          smoothed velocity (the local name inside kd_kinematics
%                      is v_smooth; kin.v_smooth does not exist)
%       kin.a          raw acceleration, masked to NaN outside
%                      [impact_index, stopFrame]
%       kin.t_stop_s   the SAME stop scalar every analysis uses
%       kin.d_final_cm the measured final depth
%       kin.impact_index, kin.stopFrame
%   a+g is recomputed via src/net_accel.m as g - a; kin.a_plus_g is never read
%   (pre-2026-08 files hold the old formula there, the raw a beside it is
%   fine). _kin.mat also stores calib, so g comes from the trial's own
%   calibration.
%
%   ══════════════════ DISPLAY LAYER -- NONE OF THIS REACHES A FIT ══════════
%
%   ENSEMBLE. Per (model, height), each trial is resampled onto one shared
%   grid relative to impact, and every grid point AFTER that trial's detected
%   t_stop is filled with the trial's own measured REST STATE: v = 0, depth =
%   its d_final, a+g = NaN. The pointwise 'Average' (median by default) is
%   then taken over the FULL replicate set at every grid point, and the curve
%   terminates at that height's MEDIAN replicate t_stop -- the same endpoint
%   in all three panels.
%
%   This replaces an earlier "terminate below MinReplicates" rule, which was
%   survivor-biased: replicates stop at different times, so past the earliest
%   stop the aggregate was taken over only the longest-lasting trials, and the
%   curve was pulled toward them before ending early. Filling the rest state
%   keeps every replicate in the aggregate for its whole duration. A
%   consequence worth checking in the render: v now decays to ~0 at the
%   endpoint by construction of the median -- about half the replicates are
%   already at rest there -- rather than by being forced to zero.
%
%   a+g is deliberately NOT rest-extended. The pipeline does not compute
%   post-stop acceleration, and inventing a value (0, or g) would be
%   fabricating data rather than displaying it.
%
%   PRE-IMPACT a+g is computed HERE, display-only, from each trial's measured
%   velocity by the same local line-fit differentiation the pipeline uses
%   ('PreImpactAccelMs' window). kd_kinematics deliberately clamps its stored
%   a to [impact, stop]; this does not alter that, is never written to any kin
%   file, and never reaches a fit. During free fall these values sit near
%   zero, which is what puts the rising edge on the plot.
%
%   SMOOTHING is a movmean applied AFTER aggregation, window per row
%   ('SmoothMs'), NaN-aware and edge-shrinking: at the first post-impact
%   sample the window uses only the samples that exist, so the leading edge is
%   neither shifted nor truncated. The pre-smoothing NaN mask is re-applied
%   afterwards, which is also what stops a line drawing a connector to the
%   axis at a curve's start or end.
%
%   ZERO-DROP TRIALS are excluded by default: they are reserved for scalar
%   measurements and are not shown in these time histories. 'ShowZeroDrop'
%   restores them, including the individual-intrusion overlay rule.
%
%   Base MATLAB only.

% ── options ──────────────────────────────────────────────────────────────
opt.Root               = '';
opt.Models             = ["Default","Tight","Wide"];
opt.PreCapMs           = 20;
opt.GridMs             = 0.2;
opt.MinReplicates      = 3;
opt.Average            = 'median';
opt.SmoothMs           = struct('z',1.0, 'v',2.0, 'ag',1.5);
opt.AccelYMin          = 1e2;
opt.PreImpactAccelMs   = 1.0;
opt.ShowZeroDrop       = false;
opt.IntrusionThreshCm  = 0.1;
opt.Layout             = 'per-model';
opt.Style              = 'mean';
opt.Diagnose           = false;
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
fprintf('  Layout / Style    : %s / %s\n', opt.Layout, opt.Style);
fprintf('  PreCapMs          : %g ms\n', opt.PreCapMs);
fprintf('  GridMs            : %g ms\n', opt.GridMs);
fprintf('  MinReplicates     : %d (gates whether a height is plotted)\n', ...
        opt.MinReplicates);
fprintf('  Average           : %s\n', opt.Average);
fprintf('  SmoothMs          : depth %g ms, v %g ms, a+g %g ms\n', ...
        opt.SmoothMs.z, opt.SmoothMs.v, opt.SmoothMs.ag);
fprintf('  AccelYMin         : %g cm/s^2 (log-axis floor)\n', opt.AccelYMin);
fprintf('  PreImpactAccelMs  : %g ms (display-only pre-impact a fit)\n', ...
        opt.PreImpactAccelMs);
fprintf('  ShowZeroDrop      : %s\n', local_tern(opt.ShowZeroDrop,'yes','no'));
fprintf('  OutDir            : %s\n', opt.OutDir);

% ── data ─────────────────────────────────────────────────────────────────
% Manual exclusions and the zero-drop SCALAR quarantine are applied by the
% loader. Windowing here always uses each trial's own stored kin.t_stop_s,
% never the table's (quarantined) t_stop_s column.
K = load_kinematics_set(root);
K = K(ismember(K.model, opt.Models), :);
if isempty(K)
    error('fig_kinematics:noTrials', 'No trials for model(s) %s under %s.', ...
          strjoin(cellstr(opt.Models), ', '), root);
end

nZeroDropAvail = sum(K.isZeroDrop);
if ~opt.ShowZeroDrop
    % Reserved for scalar measurements; not shown in the time histories.
    K = K(~K.isZeroDrop, :);
    if isempty(K)
        error('fig_kinematics:onlyZeroDrop', ...
            ['Every remaining trial is zero-drop, and those are excluded from ' ...
             'the time-history figures. Pass ShowZeroDrop true to include them.']);
    end
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
        ['No curve survived. Every group was unreadable, or had fewer than ' ...
         'MinReplicates = %d trials.'], opt.MinReplicates);
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
    Ci = C(inCol);
    for q = 1:numel(Ci)
        fprintf('           h = %4g mm : median t_stop = %7.2f ms  (n = %d)\n', ...
                Ci(q).height, Ci(q).medStopMs, Ci(q).nMax);
    end
    if ~isempty(O)
        nOut = sum([O.model] == opt.Models(c));
        if nOut > 0
            fprintf(['           %d zero-drop trial(s) overlaid ' ...
                     'individually (depth range > %g cm)\n'], ...
                    nOut, opt.IntrusionThreshCm);
        end
    end
end
if ~opt.ShowZeroDrop && nZeroDropAvail > 0
    fprintf('  (%d zero-drop trial(s) excluded: reserved for scalar measurements)\n', ...
            nZeroDropAvail);
end
fprintf(['\n  display layer: pointwise %ss over replicates with post-stop rest\n' ...
         '  extension; movmean smoothing per row; a+g on log axis with values\n' ...
         '  <= AccelYMin masked; quantitative fits use unmodified pipeline\n' ...
         '  kinematics.\n'], opt.Average);

if opt.Diagnose
    local_diagnose(C, opt);
end

% ── shared colour scale: median v0 over the plotted (nonzero) heights ────
allV0 = [C.medV0];
CLIM  = [min(allV0), max(allV0)];
if ~(diff(CLIM) > 0), CLIM = CLIM(1) + [0, max(1, abs(CLIM(1)))]; end
CMAP  = local_colormap(256);

% ── shared axis limits, computed once across every model and row ─────────
LIMS = local_shared_limits(C, O, opt);

% ── figures ──────────────────────────────────────────────────────────────
figs    = struct();
written = strings(0,1);

if strcmp(opt.Layout, 'per-model')
    for c = 1:nCol
        if ~any([C.model] == opt.Models(c)), continue; end
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
    written = [written; string(local_write_stops_csv(C, opt))];
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
%  DATA: loading, pre-impact acceleration, gridding
% ═════════════════════════════════════════════════════════════════════════
function s = local_load_kin(kinPath, opt)
%LOCAL_LOAD_KIN  One trial's traces, with a+g recomputed and the pre-impact
%   segment filled in. Returns [] if unreadable or missing a required field.
    s = [];
    if ~isfile(kinPath), return; end
    try
        L = load(kinPath, 'kin', 'calib');
    catch
        return
    end
    need = {'t_s','z','v','a','impact_index','stopFrame','t_stop_s','d_final_cm'};
    if ~isfield(L, 'kin') || ~all(isfield(L.kin, need))
        return
    end
    kin = L.kin;
    calib = [];
    if isfield(L, 'calib'), calib = L.calib; end

    s = struct();
    s.t_ms  = kin.t_s(:) * 1000;
    s.z     = kin.z(:);
    s.v     = kin.v(:);
    s.ag    = net_accel(kin, calib);           % g - a; NaN outside [impact,stop]

    % Display-only pre-impact acceleration, from the measured velocity.
    s.ag = local_fill_preimpact_ag(s.ag, kin, calib, opt.PreImpactAccelMs);

    % The SAME stop scalar every analysis uses, and the measured rest depth.
    s.t_stop_ms  = kin.t_stop_s * 1000;
    s.d_final_cm = kin.d_final_cm;
end

function ag = local_fill_preimpact_ag(ag, kin, calib, windowMs)
%LOCAL_FILL_PREIMPACT_AG  Pre-impact net acceleration, display-only.
%
%   kd_kinematics DELIBERATELY clamps its stored a to [impact_index,
%   stopFrame], and nothing here alters that: this value is computed in the
%   figure, never written back to any _kin.mat, and never reaches a fit. It
%   exists so the a+g panel shows the free-fall approach rather than starting
%   abruptly at the impact peak.
%
%   Method matches the pipeline's: a local order-1 fit of the measured
%   velocity against time, over a short symmetric window. During free fall
%   a ~ +g, so a+g = g - a sits near zero -- below the plotted floor, which is
%   exactly where the rising edge comes from.
%
%   The fit window is clamped to end at impact_index-1 so it NEVER straddles
%   the impact discontinuity. Letting it reach across would smear the impact
%   spike backwards and manufacture an onset that the data does not show.

    g = 980;
    if isstruct(calib) && isfield(calib, 'g_cm_s2'), g = calib.g_cm_s2; end

    t_s = kin.t_s(:);
    v   = kin.v(:);
    imp = kin.impact_index;
    if ~isfinite(imp) || imp < 2, return; end

    dt_ms = median(diff(t_s), 'omitnan') * 1000;
    if ~(dt_ms > 0) || ~isfinite(dt_ms), return; end
    half = max(1, round((windowMs/2) / dt_ms));

    lastPre = min(imp - 1, numel(t_s));
    for i = 1:lastPre
        lo = max(1, i - half);
        hi = min(lastPre, i + half);        % never crosses into the impact
        xx = t_s(lo:hi);  yy = v(lo:hi);
        ok = isfinite(xx) & isfinite(yy);
        if nnz(ok) < 3, continue; end
        p = polyfit(xx(ok), yy(ok), 1);     % slope = a, cm/s^2 (t in seconds)
        ag(i) = g - p(1);
    end
end

function T = local_trial_on_grid(s, grid)
%LOCAL_TRIAL_ON_GRID  One trial resampled onto the shared grid, with its own
%   measured REST STATE filling every grid point after its detected t_stop.
%
%   Before/at t_stop: the measured traces, linearly interpolated, NaN outside
%   the trial's own recorded span. After t_stop: v = 0 and depth = d_final,
%   both measured facts about that trial. a+g stays NaN -- the pipeline does
%   not compute post-stop acceleration and this figure will not invent it.

    T = struct('z',  nan(size(grid)), ...
               'v',  nan(size(grid)), ...
               'ag', nan(size(grid)));

    pre  = grid <= s.t_stop_ms;
    post = ~pre;

    T.z(pre)  = local_interp(s.t_ms, s.z,  grid(pre));
    T.v(pre)  = local_interp(s.t_ms, s.v,  grid(pre));
    T.ag(pre) = local_interp(s.t_ms, s.ag, grid(pre));

    T.z(post)  = s.d_final_cm;
    T.v(post)  = 0;
    T.ag(post) = NaN;
end

function y = local_interp(t, x, q)
    ok = isfinite(t) & isfinite(x);
    if nnz(ok) < 2, y = nan(size(q)); return; end
    y = interp1(t(ok), x(ok), q, 'linear', NaN);
end

% ═════════════════════════════════════════════════════════════════════════
%  ENSEMBLE CONSTRUCTION (Style 'mean')
% ═════════════════════════════════════════════════════════════════════════
function [C, O] = local_build_ensembles(K, opt)
%LOCAL_BUILD_ENSEMBLES  One aggregate curve per (model, height).

    C = local_empty_curve();
    O = local_empty_overlay();

    g = findgroups(K.model, K.dropHeight_mm);
    for k = 1:max(g)
        rows = find(g == k);
        isZD = all(K.isZeroDrop(rows));

        loaded = cell(numel(rows),1);
        for j = 1:numel(rows)
            loaded{j} = local_load_kin(K.kinPath(rows(j)), opt);
        end
        okLoad = ~cellfun(@isempty, loaded);
        % A trial with no finite t_stop cannot be rest-extended -- the whole
        % grid would count as post-stop -- so it is dropped here rather than
        % silently contributing a flat d_final everywhere.
        okLoad(okLoad) = cellfun(@(s) isfinite(s.t_stop_ms), loaded(okLoad));
        loaded = loaded(okLoad);
        rowsOk = rows(okLoad);

        % MinReplicates now gates only whether the height is plotted at all.
        if numel(loaded) < opt.MinReplicates, continue; end

        tStops = cellfun(@(s) s.t_stop_ms, loaded);
        medStopMs = median(tStops);

        % Grid spans out to the LONGEST replicate; every trial is defined over
        % all of it thanks to the rest extension. The curve is truncated to
        % the median stop afterwards.
        grid = (-opt.PreCapMs : opt.GridMs : max(tStops)).';
        if numel(grid) < 2, continue; end

        nT = numel(loaded);
        Z  = nan(numel(grid), nT);
        V  = nan(numel(grid), nT);
        A  = nan(numel(grid), nT);
        for j = 1:nT
            Tj = local_trial_on_grid(loaded{j}, grid);
            Z(:,j) = Tj.z;  V(:,j) = Tj.v;  A(:,j) = Tj.ag;
        end

        [zRaw, nz] = local_aggregate(Z, opt.Average);
        [vRaw, nv] = local_aggregate(V, opt.Average);
        [aRaw, na] = local_aggregate(A, opt.Average);

        % MEASURED support: how many replicates are still before their own stop
        % at each grid point. This is what the retired MinReplicates rule saw,
        % and it is the quantity that made that rule survivor-biased -- the
        % rest-extended support below is full everywhere, so it would hide the
        % effect entirely.
        nMeas    = arrayfun(@(t) sum(tStops >= t), grid);
        oldEndMs = local_old_rule_end(nMeas, grid, opt.MinReplicates);

        zC = local_smooth_row(zRaw, opt.GridMs, opt.SmoothMs.z,  medStopMs, 'depth');
        vC = local_smooth_row(vRaw, opt.GridMs, opt.SmoothMs.v,  medStopMs, 'v');
        aC = local_smooth_row(aRaw, opt.GridMs, opt.SmoothMs.ag, medStopMs, 'a+g');

        % Unified endpoint: every panel ends at this height's median t_stop.
        past = grid > medStopMs;
        zC(past) = NaN;  vC(past) = NaN;  aC(past) = NaN;

        % At the endpoint every replicate contributes (rest-extended), while
        % only about half are still in their measured phase -- which is exactly
        % why the median decays to rest there by construction rather than by
        % being forced.
        supportAtEnd     = local_support_at(nz,    grid, medStopMs);
        supportMeasAtEnd = local_support_at(nMeas, grid, medStopMs);

        C(end+1) = struct( ...                                    %#ok<AGROW>
            'model',        K.model(rowsOk(1)), ...
            'height',       K.dropHeight_mm(rowsOk(1)), ...
            'medV0',        median(K.v0_cm_s(rowsOk), 'omitnan'), ...
            'isZeroDrop',   isZD, ...
            'nMax',         nT, ...
            'medStopMs',    medStopMs, ...
            'tStopMinMs',   min(tStops), ...
            'tStopMaxMs',   max(tStops), ...
            'supportAtEnd', supportAtEnd, ...
            'supportMeasAtEnd', supportMeasAtEnd, ...
            'oldRuleEndMs', oldEndMs, ...
            'agRaw',        aRaw, ...
            't_ms', grid, 'z', zC, 'v', vC, 'ag', aC);   %#ok<AGROW>

        % ── individual zero-drop overlays (only when explicitly enabled) ──
        if isZD && opt.ShowZeroDrop
            for j = 1:numel(loaded)
                sFull = loaded{j};
                depthRange = max(sFull.z, [], 'omitnan') - min(sFull.z, [], 'omitnan');
                if ~(depthRange > opt.IntrusionThreshCm), continue; end
                Tj = local_trial_on_grid(sFull, grid);
                zO = local_smooth_row(Tj.z, opt.GridMs, opt.SmoothMs.z, medStopMs, 'depth (overlay)');
                vO = local_smooth_row(Tj.v, opt.GridMs, opt.SmoothMs.v, medStopMs, 'v (overlay)');
                zO(past) = NaN;  vO(past) = NaN;
                O(end+1) = struct('model', K.model(rowsOk(j)), ...
                    'medV0', median(K.v0_cm_s(rowsOk), 'omitnan'), ...
                    't_ms', grid, 'z', zO, 'v', vO); %#ok<AGROW>
            end
        end
    end
end

function C = local_build_trials(K, opt)
%LOCAL_BUILD_TRIALS  Every kept trial as its own curve, unaveraged. QA view.
%   Same struct shape as the ensemble curves so the drawing code is shared.
%   No rest extension here: this view exists to show what each trial actually
%   recorded, ending at its own stop.
    C = local_empty_curve();
    for i = 1:height(K)
        s = local_load_kin(K.kinPath(i), opt);
        if isempty(s), continue; end

        idx = s.t_ms >= -opt.PreCapMs & s.t_ms <= s.t_stop_ms;
        if nnz(idx) < 2, continue; end
        t = s.t_ms(idx);

        zC = local_smooth_row(s.z(idx),  [], opt.SmoothMs.z,  s.t_stop_ms, 'depth', t);
        vC = local_smooth_row(s.v(idx),  [], opt.SmoothMs.v,  s.t_stop_ms, 'v',     t);
        aC = local_smooth_row(s.ag(idx), [], opt.SmoothMs.ag, s.t_stop_ms, 'a+g',   t);

        C(end+1) = struct( ...                                    %#ok<AGROW>
            'model',        K.model(i), ...
            'height',       K.dropHeight_mm(i), ...
            'medV0',        K.v0_cm_s(i), ...
            'isZeroDrop',   K.isZeroDrop(i), ...
            'nMax',         1, ...
            'medStopMs',    s.t_stop_ms, ...
            'tStopMinMs',   s.t_stop_ms, ...
            'tStopMaxMs',   s.t_stop_ms, ...
            'supportAtEnd', 1, ...
            'supportMeasAtEnd', 1, ...
            'oldRuleEndMs', NaN, ...
            'agRaw',        s.ag(idx), ...
            't_ms', t, 'z', zC, 'v', vC, 'ag', aC);
    end
end

function C = local_empty_curve()
    C = struct('model',{}, 'height',{}, 'medV0',{}, 'isZeroDrop',{}, ...
               'nMax',{}, 'medStopMs',{}, 'tStopMinMs',{}, 'tStopMaxMs',{}, ...
               'supportAtEnd',{}, 'supportMeasAtEnd',{}, 'oldRuleEndMs',{}, ...
               'agRaw',{}, ...
               't_ms',{}, 'z',{}, 'v',{}, 'ag',{});
end

function O = local_empty_overlay()
    O = struct('model',{}, 'medV0',{}, 't_ms',{}, 'z',{}, 'v',{});
end

% ═════════════════════════════════════════════════════════════════════════
%  AGGREGATION AND SMOOTHING
% ═════════════════════════════════════════════════════════════════════════
function [m, n] = local_aggregate(Y, avg)
%LOCAL_AGGREGATE  Pointwise median/mean down the replicate dimension.
%   No coverage floor is applied: with the rest extension every replicate
%   contributes at every grid point it can, and the curve's endpoint is set
%   by the median stop rather than by where support runs out.
    n = sum(~isnan(Y), 2);
    if strcmp(avg, 'median')
        m = median(Y, 2, 'omitnan');
    else
        m = mean(Y, 2, 'omitnan');
    end
    m(n == 0) = NaN;
end

function tEnd = local_old_rule_end(n, grid, minRep)
%LOCAL_OLD_RULE_END  Where the retired "terminate below MinReplicates" rule
%   would have ended this curve. Diagnostic only -- it documents how much of
%   the curve was previously lost to survivor bias.
    tEnd = NaN;
    postImpact = grid >= 0;
    idx = find(postImpact & n < minRep, 1, 'first');
    if ~isempty(idx), tEnd = grid(idx); end
end

function sAt = local_support_at(n, grid, tMs)
%LOCAL_SUPPORT_AT  Replicate count at the grid point nearest tMs.
    sAt = NaN;
    if isempty(grid) || ~isfinite(tMs), return; end
    [~, i] = min(abs(grid - tMs));
    sAt = n(i);
end

function y = local_smooth_row(yRaw, gridMs, windowMs, medianStopMs, rowLabel, tMs)
%LOCAL_SMOOTH_ROW  movmean AFTER aggregation, then re-mask to the original
%   support so smoothing cannot extend a curve past where the data ends.
%
%   movmean's default 'shrink' endpoint rule plus 'omitnan' is exactly what is
%   wanted at the impact edge: the window at the first post-impact sample
%   spans only the samples that exist, so the leading edge is neither shifted
%   later nor truncated. Re-applying the input NaN mask afterwards keeps the
%   curve's support identical to the data's.
    if isempty(gridMs)
        if nargin < 6 || numel(tMs) < 2, y = yRaw; return; end
        gridMs = median(diff(tMs), 'omitnan');
    end
    if ~(gridMs > 0) || ~isfinite(gridMs), y = yRaw; return; end

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
%  DIAGNOSTICS
% ═════════════════════════════════════════════════════════════════════════
function local_diagnose(C, opt)
%LOCAL_DIAGNOSE  Per-height numbers behind the ensemble, and the a+g onset.
%
%   Documents two things the render alone cannot show: how far apart the
%   replicate stop times are (the mechanism the retired MinReplicates rule
%   turned into survivor bias), and whether display smoothing moves the a+g
%   onset (it should not -- the movmean shrinks its window at the edge).

    fprintf('\n=== DIAGNOSE ===\n');
    models = unique([C.model], 'stable');
    for mi = 1:numel(models)
        sel = find([C.model] == models(mi));
        if isempty(sel), continue; end
        fprintf('\n  %s\n', models(mi));
        fprintf(['    %5s %4s  %-25s %-11s %-11s %s\n'], ...
                'h(mm)', 'n', 't_stop min/med/max (ms)', 'supp@end', ...
                'meas@end', 'old-rule end (ms)');
        for q = sel
            oldStr = 'never';   % the retired rule would not have truncated
            if isfinite(C(q).oldRuleEndMs)
                oldStr = sprintf('%.2f', C(q).oldRuleEndMs);
            end
            fprintf('    %5g %4d  %7.2f /%7.2f /%7.2f  %8d    %8d    %s\n', ...
                    C(q).height, C(q).nMax, ...
                    C(q).tStopMinMs, C(q).medStopMs, C(q).tStopMaxMs, ...
                    C(q).supportAtEnd, C(q).supportMeasAtEnd, oldStr);
        end
        fprintf(['      supp@end = replicates contributing at the endpoint ' ...
                 '(all, via rest\n      extension); meas@end = those still ' ...
                 'before their own stop. old-rule end\n      = where the ' ...
                 'retired MinReplicates rule would have truncated the curve.\n']);

        % One representative mid-height: the a+g onset, before vs after the
        % display smooth. If smoothing were shifting the onset, the first
        % non-NaN sample would move between these two rows.
        hs = [C(sel).height];
        [~, imid] = min(abs(hs - median(hs)));
        q = sel(imid);
        fprintf('\n    a+g onset at h = %g mm (first 5 post-impact grid samples)\n', ...
                C(q).height);
        gi = find(C(q).t_ms >= 0, 1, 'first');
        if isempty(gi)
            fprintf('      (no post-impact samples)\n');
        else
            ii = gi:min(gi+4, numel(C(q).t_ms));
            fprintf('      t (ms)   : %s\n', local_fmt(C(q).t_ms(ii)));
            fprintf('      raw      : %s\n', local_fmt(C(q).agRaw(ii)));
            fprintf('      smoothed : %s\n', local_fmt(C(q).ag(ii)));
            fprintf(['      (smoothing window %g ms; movmean shrinks at the ' ...
                     'edge, so the\n       first sample is not shifted)\n'], ...
                    opt.SmoothMs.ag);
        end
    end
    fprintf('\n');
end

function s = local_fmt(v)
    s = strjoin(arrayfun(@(x) sprintf('%9.2f', x), v(:).', ...
                         'UniformOutput', false), ' ');
end

% ═════════════════════════════════════════════════════════════════════════
%  LIMITS AND COLOUR
% ═════════════════════════════════════════════════════════════════════════
function LIMS = local_shared_limits(C, O, opt)
%LOCAL_SHARED_LIMITS  One set of x/y limits per row, shared across every model
%   and every figure, so the panels stay comparable.

    allT = vertcat(C.t_ms);
    if ~isempty(O), allT = [allT; vertcat(O.t_ms)]; end
    XL = [min(allT), max(allT)];

    zAll = vertcat(C.z);
    if ~isempty(O), zAll = [zAll; vertcat(O.z)]; end
    YL_z = local_pad_ylim(zAll);

    vAll = vertcat(C.v);
    if ~isempty(O), vAll = [vAll; vertcat(O.v)]; end
    YL_v = local_pad_ylim(vAll);
    YL_v(1) = 0;                 % velocity floor is zero: nothing below rest

    YL_ag = local_log_ylim(vertcat(C.ag), opt.AccelYMin);

    LIMS = struct('t', XL, 'z', YL_z, 'v', YL_v, 'ag_log', YL_ag);
end

function yl = local_pad_ylim(vals)
    vals = vals(isfinite(vals));
    if isempty(vals), yl = [-1 1]; return; end
    lo = min(vals); hi = max(vals); r = hi - lo;
    if r == 0, r = max(abs(hi), 1); end
    yl = [lo - 0.05*r, hi + 0.05*r];
end

function yl = local_log_ylim(vals, floorVal)
%LOCAL_LOG_YLIM  Log-axis bounds. The lower limit is the configured floor --
%   set by this rig's noise floor rather than by KD's 1e1 -- and the upper is
%   rounded up to cover the data.
    vals = vals(isfinite(vals) & vals > 0);
    hi = 3e4;
    if ~isempty(vals), hi = max(hi, 10^ceil(log10(max(vals)))); end
    yl = [floorVal, hi];
end

function cmap = local_colormap(nLevels)
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
    fig = figure('Color','w','Units','inches','Position',[1 1 3.4 7.0], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, 3, 1, 'Padding','compact', 'TileSpacing','compact');

    ax1 = nexttile(tl); local_draw_panel(ax1, 'z',  model, C, O, CLIM, CMAP, LIMS, opt);
    ax2 = nexttile(tl); local_draw_panel(ax2, 'v',  model, C, O, CLIM, CMAP, LIMS, opt);
    ax3 = nexttile(tl); local_draw_panel(ax3, 'ag', model, C, O, CLIM, CMAP, LIMS, opt);

    xlabel(ax3, 't - t_impact (ms)', 'FontSize', 8, 'Interpreter', 'none');
    set([ax1 ax2], 'XTickLabel', []);
end

function fig = local_draw_grid3x3(models, C, O, CLIM, CMAP, LIMS, opt)
    nCol = numel(models);
    fig = figure('Color','w','Units','inches','Position',[1 1 7.0 6.4], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl  = tiledlayout(fig, 3, nCol, 'Padding','compact', 'TileSpacing','compact');

    ROWS = {'z','v','ag'};
    for r = 1:3
        for c = 1:nCol
            ax = nexttile(tl, (r-1)*nCol + c);
            local_draw_panel(ax, ROWS{r}, models(c), C, O, CLIM, CMAP, LIMS, opt);
            if r < 3, set(ax, 'XTickLabel', []); end
            if c > 1, set(ax, 'YTickLabel', []); end
            if r == 3 && c == 1
                xlabel(ax, 't - t_impact (ms)', 'FontSize', 8, 'Interpreter', 'none');
            end
        end
    end
end

function local_draw_panel(ax, rowKey, model, C, O, CLIM, CMAP, LIMS, opt)
%LOCAL_DRAW_PANEL  One (row, model) panel -- shared by both layouts and both
%   styles, since C/O are already in a common shape.

    hold(ax, 'on');

    % Drawn FIRST so they sit behind the data.
    xline(ax, 0, '--', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.75, ...
          'HandleVisibility', 'off');
    if strcmp(rowKey, 'z')
        % Solid thin line at depth = 0: the bed surface, a physical datum.
        yline(ax, 0, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end

    mkMap = containers.Map({'Default','Tight','Wide'}, {'o','s','^'});
    % Deliberately different from depth_scaling.m's Tight=^/Wide=s convention:
    % this figure's marker mapping is specified independently in its own brief.

    inCol = find([C.model] == model);
    for k = inCol
        y = C(k).(rowKey);
        if isempty(y) || all(isnan(y)), continue; end
        if strcmp(rowKey, 'ag')
            y = local_prep_log(y, opt.AccelYMin);
            if all(isnan(y)), continue; end
        end
        col = local_v0_color(C(k).medV0, CLIM, CMAP);
        plot(ax, C(k).t_ms, y, '-', 'LineWidth', 1.5, 'Color', col);

        if strcmp(rowKey, 'z')
            % Stop marker at the unified endpoint: the height's median t_stop.
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

    if ismember(rowKey, {'z','v'}) && ~isempty(O)
        for k = find([O.model] == model)
            y = O(k).(rowKey);
            if isempty(y) || all(isnan(y)), continue; end
            col = local_v0_color(O(k).medV0, CLIM, CMAP);
            plot(ax, O(k).t_ms, y, '-', 'LineWidth', 0.75, 'Color', col);
        end
    end

    switch rowKey
        case 'ag'
            set(ax, 'YScale', 'log');
            ylim(ax, LIMS.ag_log);
            ylabel(ax, 'a + g (cm/s^2)', 'FontSize', 8, 'Interpreter', 'none');
        case 'v'
            ylim(ax, LIMS.v);
            ylabel(ax, 'v (cm/s)', 'FontSize', 8, 'Interpreter', 'none');
        otherwise
            ylim(ax, LIMS.z);
            ylabel(ax, 'depth (cm)', 'FontSize', 8, 'Interpreter', 'none');
    end
    xlim(ax, LIMS.t);

    box(ax, 'on');
    grid(ax, 'on');
    set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
            'FontSize', 8, 'LineWidth', 0.5, 'Layer', 'top');
end

function yPlot = local_prep_log(y, floorVal)
%LOCAL_PREP_LOG  Prepare a curve for the log axis.
%
%   Sub-floor samples are masked to NaN so they do not crawl along the bottom
%   of the panel -- EXCEPT the single sample immediately before the first
%   above-floor sample, which is CLIPPED to the floor instead of masked. That
%   one retained point is what lets MATLAB draw the rising connector from the
%   axis floor up into the first resolved value: without it the curve would
%   begin abruptly at the peak, and with a NaN there the connector would be
%   broken.
%
%   Pre-impact samples are near zero (free fall), so they sit below the floor
%   and this transition is exactly the impact onset.
%
%   CONSIDERED AND REJECTED: inserting a synthetic onset point at a fixed
%   value (e.g. 1e1) to force a visible rise. Every plotted sample here is
%   measured; the sub-resolution part of the rise is represented by the
%   clipped connector, and no fabricated datum appears in the figure.

    yPlot = y;
    sub   = ~isfinite(y) | y <= floorVal;
    first = find(~sub, 1, 'first');
    yPlot(sub) = NaN;
    if ~isempty(first) && first > 1
        prev = first - 1;
        if isfinite(y(prev))       % a measured sample that sits below the floor
            yPlot(prev) = floorVal;
        end
    end
end

% ═════════════════════════════════════════════════════════════════════════
function csvPath = local_write_stops_csv(C, opt)
%LOCAL_WRITE_STOPS_CSV  Per-height reportable scalars for the cross-geometry
%   comparison: the median stop time and median v0 behind each plotted curve.
    if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
    csvPath = fullfile(opt.OutDir, 'fig1_kinematics_stops.csv');
    fid = fopen(csvPath, 'w');
    if fid < 0
        error('fig_kinematics:csvFailed', 'Could not write %s', csvPath);
    end
    fprintf(fid, '# written by scripts/fig_kinematics.m on %s\n', ...
            datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '# median over replicates; t_stop from each trial''s kin.t_stop_s\n');
    fprintf(fid, 'model,height_mm,median_v0_cm_s,median_t_stop_ms,n\n');
    for k = 1:numel(C)
        fprintf(fid, '%s,%g,%.4f,%.4f,%d\n', ...
                C(k).model, C(k).height, C(k).medV0, C(k).medStopMs, C(k).nMax);
    end
    fclose(fid);
end

function paths = local_export(fig, outDir, stem)
    if ~isfolder(outDir), mkdir(outDir); end
    pdfPath = fullfile(outDir, [stem '.pdf']);
    pngPath = fullfile(outDir, [stem '.png']);
    exportgraphics(fig, pdfPath, 'ContentType', 'vector');
    exportgraphics(fig, pngPath, 'Resolution', 600);
    paths = [string(pdfPath); string(pngPath)];
end

function s = local_tern(c,a,b), if c, s=a; else, s=b; end, end
