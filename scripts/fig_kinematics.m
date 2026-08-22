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
%       'PreCapMs'          pre-impact context cap, ms (default 10). Keep it
%                           an exact multiple of GridMs (10/0.2 = 50) so t = 0
%                           lands ON a grid point: that is what makes the
%                           first ensemble a+g sample the impact sample itself
%                           rather than a point already partway up the rise.
%       'GridMs'            ensemble time step, ms (default 0.2)
%       'MinReplicates'     minimum trials for a height to be PLOTTED AT ALL
%                           (default 3). It no longer decides where a curve
%                           ends -- see ENSEMBLE below.
%       'Average'           'median' (default) | 'mean', across replicates
%       'SmoothMs'          struct('z',1.0,'v',2.0,'ag',1.5): display movmean
%                           window per row, ms, applied AFTER aggregation
%       'AccelScale'        a+g y-axis: 'linear' (default) | 'log'
%       'AccelYMin'         a+g LOG-axis floor, cm/s^2 (default 1e2); ignored
%                           on the linear axis
%       'ShowZeroDrop'      include zero-drop trials (default FALSE -- they
%                           are reserved for scalar measurements)
%       'IntrusionThreshCm' zero-drop depth-range threshold for an individual
%                           overlay (default 0.1 cm; only used when
%                           ShowZeroDrop is true)
%       'ColorBy'           'height' (default) | 'v0'. What the curve colour
%                           encodes. 'height' ranks the nominal drop height on
%                           the ladder shared by all three geometries, so the
%                           three figures step through identical colours by
%                           construction. 'v0' maps the median impact speed
%                           through one shared range. Neither touches the
%                           measured v0, which is the velocity panel's data.
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
%   a+g is deliberately NOT rest-extended across the grid. The pipeline does
%   not compute post-stop acceleration, and filling a value forward would be
%   fabricating a trace rather than displaying one.
%
%   EVERY a+g VALUE BETWEEN THE ANCHORS IS MEASURED. kd_kinematics masks a to
%   NaN outside [impact_index, stopFrame], so the computed part of the panel
%   spans exactly that interval; nothing inside it is reconstructed, and no
%   pre-impact acceleration is estimated from the velocity or anywhere else.
%   The two boundary anchors below are the only drawn values that do not come
%   from the trace, and both are exact rather than inferred.
%
%   BOUNDARY ANCHORS AT BOTH ENDS. The a+g curve is drawn between two states
%   that are known exactly without being recorded, because the definition
%   a + g = g - a fixes them:
%
%       onset, as t -> 0-    free fall, a = g   ->  a + g = 0
%       stop, at t_stop      at rest,   a = 0   ->  a + g = g
%
%   Neither is in the stored trace -- kd_kinematics computes a only over
%   [impact_index, stopFrame] -- and neither is a measurement, a
%   reconstruction or a placeholder. They are boundary values, and the curve
%   is anchored to them: a leading segment from (0, 0) to the first computed
%   sample, and a trailing segment from the last computed sample down to g.
%
%   What this buys: the steep rise is shown as the measured jump from its true
%   starting state rather than materialising at ~g out of nowhere, and the drop
%   across the terminal segment reads as the per-height acceleration
%   discontinuity -- how much net upward acceleration the grains were still
%   supplying when the rod stopped -- against a dotted reference line at g.
%
%   VELOCITY IS ANCHORED THE SAME WAY, to its own rest value. At the unified
%   stop the median v ends small-POSITIVE by construction -- about half the
%   replicates are already at rest at 0 there and the rest are still moving,
%   so the median sits just above zero -- and a final segment to (stop, 0)
%   closes that. Same helper, same endpoint. Unlike the a+g drop, this one
%   carries no physical content: it removes a construction artefact rather
%   than displaying a measured discontinuity.
%
%   The terminal segment ends at the SAME t as the depth and velocity curves
%   (one stop per height per geometry, shared by that figure's three panels)
%   and nothing is drawn past it. The three panels do NOT share a start, and
%   are not meant to. 'Diagnose' asserts the ENDPOINTS are identical.
%
%   AXIS CHOICE ('AccelScale', default 'linear'). The plotted a+g spans roughly
%   1e3-2.2e4 cm/s^2 -- about one decade -- so a linear axis preserves the
%   two-bump structure and the drop near the stop that a log axis compresses
%   into near-invisibility. Log axes earn their place across several decades:
%   KD 2007 spanned 1e1-1e4 at 20 us sampling, which this rig's ~1 ms
%   acceleration resolution does not reach. 'log' is kept as a one-argument
%   fallback should a reviewer ask for KD-style axes.
%
%   SMOOTHING is a movmean applied AFTER aggregation, window per row
%   ('SmoothMs'), NaN-aware. The pre-smoothing NaN mask is re-applied
%   afterwards, which is also what stops a line drawing a connector to the
%   axis at a curve's start or end.
%
%   The a+g row is smoothed EDGE-PRESERVING; depth and velocity are not, and
%   do not need to be. movmean's 'shrink' rule truncates the window at a
%   boundary but still averages about half a window there. Depth and velocity
%   begin PreCapMs before impact, where the curve is flat, so that costs
%   nothing. a+g begins AT impact, on the steepest millisecond of the trace,
%   where it cost everything: at the impact sample v is maximal, so a is ~0
%   and a+g is ~g, and averaging that measured onset together with the rise
%   above it reported the mid-rise value (0.3-1.6e4) as the starting point
%   instead of ~0.1e4. Edge preservation sets the window's half-width to the
%   distance from the edge of the finite support, so the first and last
%   samples are returned unchanged and the window stays symmetric -- nothing
%   is phase-shifted. 'Diagnose' asserts the first sample is untouched.
%
%   COLOUR ('ColorBy') encodes one of two things, through a reference built
%   ONCE from the FULL loaded set, before the 'Models' filter -- so it cannot
%   depend on which geometries a given call happens to draw. Regenerating one
%   figure on its own must not remap it.
%
%     'height' (default)  the rung of the nominal drop height on the ladder
%                         shared by all three geometries: 25 mm coldest,
%                         365 mm warmest. Because the index is a RANK on a
%                         fixed ladder, the three figures step through
%                         identical colours by construction, whatever v0 each
%                         geometry happened to reach at a given height.
%     'v0'                median impact speed through one shared range, the
%                         earlier behaviour, retained unchanged.
%
%   Nothing about this modifies a measured v0. It remains the plotted data in
%   the velocity panel, the medV0 on every curve, and the number every
%   analysis uses. Zero-drop curves are excluded from both references: their
%   v0 is 0 by quarantine rather than by measurement, which would anchor the
%   v0 scale at zero and put a phantom rung at the foot of the ladder.
%
%   The reference is printed in the run summary, and each figure's extreme
%   curves are printed as (height or v0) -> colormap index, so the shared-
%   reference claim is checkable from the log alone.
%
%   SUPPORT GATES AT BOTH ENDS. A row is drawn only where at least
%   MinReplicates replicates contribute to it, the same threshold that decides
%   whether a height is plotted at all.
%
%   Before impact this applies to ALL THREE rows. Trials carry different
%   amounts of pre-impact recording -- Stage A's window starts where the red
%   markers appear -- so the pointwise median steps as the contributing set
%   changes, which is the kink that was visible in the Default velocity curve
%   near -5 ms. It is an artefact of who is in the average, not of the motion.
%
%   After impact only a+g needs it, and the reason is below. Depth and
%   velocity are rest-extended, so their support is full there by
%   construction.
%
%   Neither gate can punch an interior hole: support is non-increasing away
%   from t = 0 in both directions -- going back, a trial drops out at the
%   start of its own recording; going forward, at its own stop -- so each gate
%   can only trim a contiguous head or tail.
%
%   THE a+g TAIL IS SUPPORT-GATED. a+g is never rest-extended, so past the
%   earliest replicate stop the pointwise median is taken over fewer and fewer
%   trials and starts to wander. The computed curve is therefore drawn only
%   where at least MinReplicates replicates still contribute; the terminal
%   segment runs from that last well-supported sample to (unified stop, g).
%   The gate always removes a contiguous tail, never a mid-curve hole, since
%   the count is non-increasing after impact. Nothing is clamped: a+g below g
%   mid-penetration is physical and is preserved exactly.
%
%   ZERO-DROP TRIALS are excluded by default: they are reserved for scalar
%   measurements and are not shown in these time histories. 'ShowZeroDrop'
%   restores them, including the individual-intrusion overlay rule.
%
%   Base MATLAB only.

% ── options ──────────────────────────────────────────────────────────────
opt.Root               = '';
opt.Models             = ["Default","Tight","Wide"];
opt.PreCapMs           = 10;
opt.GridMs             = 0.2;
opt.MinReplicates      = 3;
opt.Average            = 'median';
opt.SmoothMs           = struct('z',1.0, 'v',2.0, 'ag',1.5);
opt.AccelYMin          = 1e2;
opt.AccelScale         = 'linear';
opt.ShowZeroDrop       = false;
opt.IntrusionThreshCm  = 0.1;
opt.Layout             = 'per-model';
opt.Style              = 'mean';
opt.ColorBy            = 'height';
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
opt.AccelScale = lower(char(opt.AccelScale));
if ~ismember(opt.AccelScale, {'linear','log'})
    error('fig_kinematics:badAccelScale', ...
          'AccelScale must be ''linear'' or ''log'', got "%s".', opt.AccelScale);
end
opt.ColorBy = lower(char(opt.ColorBy));
if ~ismember(opt.ColorBy, {'height','v0'})
    error('fig_kinematics:badColorBy', ...
          'ColorBy must be ''height'' or ''v0'', got "%s".', opt.ColorBy);
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
fprintf('  AccelScale        : %s%s\n', opt.AccelScale, ...
        local_tern(strcmp(opt.AccelScale,'log'), ...
                   sprintf(' (floor %g cm/s^2)', opt.AccelYMin), ''));
fprintf('  ColorBy           : %s\n', opt.ColorBy);
fprintf('  ShowZeroDrop      : %s\n', local_tern(opt.ShowZeroDrop,'yes','no'));
fprintf('  OutDir            : %s\n', opt.OutDir);

% ── data ─────────────────────────────────────────────────────────────────
% Manual exclusions and the zero-drop SCALAR quarantine are applied by the
% loader. Windowing here always uses each trial's own stored kin.t_stop_s,
% never the table's (quarantined) t_stop_s column.
% KAll is the FULL loaded set, every model, kept in scope on purpose: the
% shared colour range is derived from it rather than from the models this
% invocation happens to draw. See local_global_color_ref.
KAll = load_kinematics_set(root);
K    = KAll(ismember(KAll.model, opt.Models), :);
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
         '  extension; movmean smoothing per row, edge-preserving on a+g so the\n' ...
         '  measured onset at impact is not averaged upward; a+g on a %s axis,\n' ...
         '  anchored at both ends to its boundary states (free fall a+g = 0 at\n' ...
         '  t = 0; rest a+g = g at the unified stop shared with depth and v);\n' ...
         '  quantitative fits use unmodified pipeline kinematics.\n'], ...
         opt.Average, opt.AccelScale);

if opt.Diagnose
    local_diagnose(C, opt);
end

% ── shared colour reference, ONE for every geometry ──────────────────────
% Derived from KAll -- the full loaded set, BEFORE the 'Models' filter -- so
% it is a property of the dataset, not of which figures this call is drawing.
% Computed once, here, and handed unchanged to every figure.
REF  = local_global_color_ref(KAll, opt);
CMAP = local_colormap(256);
if strcmp(opt.ColorBy, 'height')
    fprintf(['  colour by             : nominal drop height, %d rungs ' ...
             '%g - %g mm\n                          (shared ladder; ' ...
             'measured v0 is unchanged and is the v panel)\n'], ...
            numel(REF.ladder), min(REF.ladder), max(REF.ladder));
else
    fprintf(['  colour by             : median v0, %.1f - %.1f cm/s over ' ...
             '%d (model, height) group(s)\n                          ' ...
             'across ALL geometries, independent of the Models option\n'], ...
            REF.clim(1), REF.clim(2), REF.nGroups);
end

% ── shared axis limits, computed once across every model and row ─────────
LIMS = local_shared_limits(C, O, opt);
if strcmp(opt.AccelScale, 'log')
    fprintf('  shared a+g y-limits   : %.4g - %.4g cm/s^2 (log axis)\n', ...
            LIMS.ag_log(1), LIMS.ag_log(2));
else
    fprintf(['  shared a+g y-limits   : %.4g - %.4g cm/s^2 ' ...
             '(floor held at 0; top = 1.05 x global max)\n'], ...
            LIMS.ag_lin(1), LIMS.ag_lin(2));
end
fprintf('  shared t-limits       : %.4g - %.4g ms\n', LIMS.t(1), LIMS.t(2));

% ── figures ──────────────────────────────────────────────────────────────
figs    = struct();
written = strings(0,1);

% Printed per figure so the shared-range claim is checkable from the log
% alone: the same v0 must map to the same index everywhere, and every line
% must quote the same clim.
fprintf('\n--- colour check (v0 -> colormap index, %d levels) ---\n', size(CMAP,1));

if strcmp(opt.Layout, 'per-model')
    for c = 1:nCol
        if ~any([C.model] == opt.Models(c)), continue; end
        local_print_color_check(opt.Models(c), C, REF, CMAP, opt);
        f = local_draw_one_model(opt.Models(c), C, O, REF, CMAP, LIMS, opt);
        figs.(matlab.lang.makeValidName(char(opt.Models(c)))) = f;
        if opt.Save
            stem = sprintf('fig1_kinematics_%s', lower(char(opt.Models(c))));
            written = [written; local_export(f, opt.OutDir, stem)]; %#ok<AGROW>
        end
    end
else % grid3x3
    for c = 1:nCol
        local_print_color_check(opt.Models(c), C, REF, CMAP, opt);
    end
    fg = local_draw_grid3x3(opt.Models, C, O, REF, CMAP, LIMS, opt);
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
R.colorBy  = opt.ColorBy;
R.clim     = REF.clim;         % v0 range; populated in both modes
R.ladder   = REF.ladder;       % shared height ladder; populated in both modes
R.limits   = LIMS;
R.written  = written;
R.figures  = figs;
end

% ═════════════════════════════════════════════════════════════════════════
%  DATA: loading, pre-impact acceleration, gridding
% ═════════════════════════════════════════════════════════════════════════
function s = local_load_kin(kinPath)
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
    % a+g exists ONLY where the pipeline computes it: kd_kinematics masks a to
    % NaN outside [impact_index, stopFrame], and this figure shows exactly that
    % span. Pre-impact a+g is not reconstructed here.
    % g comes back from net_accel rather than being re-derived from calib, so
    % the rest state drawn at a + g = g cannot drift from the value the trace
    % itself was computed with.
    [s.ag, s.g_cm_s2] = net_accel(kin, calib);  % g - a; NaN outside [impact,stop]

    % The SAME stop scalar every analysis uses, and the measured rest depth.
    s.t_stop_ms  = kin.t_stop_s * 1000;
    s.d_final_cm = kin.d_final_cm;
end

function [T, shortRec] = local_trial_on_grid(s, grid)
%LOCAL_TRIAL_ON_GRID  One trial resampled onto the shared grid, rest-filled all
%   the way to the END OF THE GRID.
%
%   The rest state (v = 0, depth = this trial's d_final) starts at the trial's
%   REST ONSET and continues to the grid's end -- not merely to the trial's
%   last recorded sample. That distinction is the whole point: a trial whose
%   recording ends before the longest t_stop at its height used to contribute
%   NaN over the remainder, and where enough replicates did that the pointwise
%   median went NaN too, leaving a visible mid-curve gap (seen at t ~ 55-63 ms
%   in the lowest-v0 curves). Filling to the grid end makes every trial
%   contribute everywhere after its own onset, so the median is continuous from
%   onset to the ensemble endpoint by construction.
%
%   REST ONSET is min(t_stop, last recorded sample). Normally that is t_stop.
%   When a recording ends first, treating the trial as at rest from that point
%   is an APPROXIMATION -- the trial may still have been moving -- so those
%   trials are counted and reported by Diagnose rather than passing silently.
%
%   a+g is never rest-filled: kd_kinematics computes no post-stop acceleration
%   and this figure will not invent one.

    T = struct('z',  nan(size(grid)), ...
               'v',  nan(size(grid)), ...
               'ag', nan(size(grid)));

    finiteMeas = isfinite(s.t_ms) & (isfinite(s.z) | isfinite(s.v));
    if ~any(finiteMeas)
        shortRec = false;
        return
    end
    tLastMeas = max(s.t_ms(finiteMeas));

    restOnset = min(s.t_stop_ms, tLastMeas);
    shortRec  = tLastMeas < s.t_stop_ms;

    meas = grid <= restOnset;
    rest = ~meas;

    T.z(meas)  = local_interp(s.t_ms, s.z,  grid(meas));
    T.v(meas)  = local_interp(s.t_ms, s.v,  grid(meas));
    T.ag(meas) = local_interp(s.t_ms, s.ag, grid(meas));

    T.z(rest)  = s.d_final_cm;
    T.v(rest)  = 0;
    T.ag(rest) = NaN;
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
            loaded{j} = local_load_kin(K.kinPath(rows(j)));
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
        nShortRec = 0;
        for j = 1:nT
            [Tj, shortRec] = local_trial_on_grid(loaded{j}, grid);
            Z(:,j) = Tj.z;  V(:,j) = Tj.v;  A(:,j) = Tj.ag;
            nShortRec = nShortRec + shortRec;
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

        % ── a+g TAIL SUPPORT GATE ────────────────────────────────────────
        % na is how many replicates have a COMPUTED a+g at each grid point.
        % Unlike depth and v, a+g is never rest-extended -- the pipeline
        % computes no post-stop acceleration -- so past the earliest replicate
        % stop this count falls away, and the median is taken over fewer and
        % fewer trials until it is wandering over two or three. That is what
        % produced the dip to ~665 cm/s^2, below g, just before a terminal
        % drop: not physics, just a median with almost nothing left in it.
        %
        % So the computed a+g is drawn only where support >= MinReplicates,
        % the same threshold the figure already uses to decide whether a
        % height is worth plotting at all. The gate always cuts a contiguous
        % TAIL and never punches a hole mid-curve: after impact each trial
        % contributes over [0, its own stop], so the count is non-increasing
        % in t.
        %
        % Applied to the AGGREGATE, before smoothing, so the smoother's
        % trailing edge sits at the gate rather than reaching past it into
        % thin-support samples that will not be drawn.
        %
        % THIS IS NOT A FLOOR. Nothing is clamped to g anywhere: a+g below g
        % mid-penetration is real -- the grains are decelerating the rod less
        % than gravity accelerates it -- and those values are preserved
        % exactly. Only the thin-support tail is removed.
        agThin = na < opt.MinReplicates;
        aRaw(agThin) = NaN;

        % ── PRE-IMPACT SUPPORT GATE, all three rows ──────────────────────
        % The same rule, applied to the other end. Trials do not all carry the
        % same amount of pre-impact recording: Stage A's window starts where
        % the red markers appear, so at -10 ms a height may have three trials
        % contributing where at -2 ms it has ten. The pointwise median then
        % steps as the contributing SET changes, which is the kink visible in
        % the Default velocity curve around -5 ms. It is an artefact of who is
        % in the average, not a feature of the motion.
        %
        % So each row is drawn before impact only where at least
        % MinReplicates trials contribute -- its OWN support count, since the
        % three rows need not thin out together.
        %
        % Like the tail gate, this can only trim a contiguous HEAD and can
        % never punch an interior hole: going backwards from t = 0 a trial
        % drops out at the start of its own recording and never returns, so
        % the count is non-increasing away from impact in this direction too.
        %
        % Post-impact needs no such gate on depth and velocity: the rest
        % extension makes every replicate contribute at every point after its
        % onset, so their support is full there by construction. a+g is the
        % exception and is gated over the whole grid, above.
        pre = grid < 0;
        zRaw(pre & nz < opt.MinReplicates) = NaN;
        vRaw(pre & nv < opt.MinReplicates) = NaN;
        aRaw(pre & na < opt.MinReplicates) = NaN;

        zC = local_smooth_row(zRaw, opt.GridMs, opt.SmoothMs.z,  medStopMs, 'depth');
        vC = local_smooth_row(vRaw, opt.GridMs, opt.SmoothMs.v,  medStopMs, 'v');
        % edgePreserve for a+g only: its support starts AT impact, on the
        % steepest millisecond of the trace, where a shrunk movmean window
        % still averages the measured onset up into the rise.
        aC = local_smooth_row(aRaw, opt.GridMs, opt.SmoothMs.ag, medStopMs, 'a+g', [], true);

        % Unified endpoint: every panel ends at this height's median t_stop.
        past = grid > medStopMs;
        zC(past) = NaN;  vC(past) = NaN;  aC(past) = NaN;

        % THE shared endpoint in t, for all three panels of this height: the
        % last grid point at or before the median stop. Taken once, here, so
        % the a+g terminal segment cannot end anywhere else.
        tEndMs = max(grid(~past));

        % g used for the rest state is the one net_accel used per trial; they
        % are the same constant across a height in practice, and taking the
        % median rather than the first makes a stray calib harmless.
        gRest = median(cellfun(@(s) s.g_cm_s2, loaded));

        % Where the gated computed curve stops, and therefore where the
        % terminal segment starts. local_rest_terminal takes the last finite
        % sample, so gating aRaw above is all it takes for the segment to run
        % from the last well-supported value to (tEndMs, g) -- the unified
        % stop is unchanged, and so is the endpoint assertion.
        agGateMs = NaN;
        iGate = find(isfinite(aC), 1, 'last');
        if ~isempty(iGate), agGateMs = grid(iGate); end

        % Both rows are anchored to their rest value at the unified stop, by
        % the same helper: a+g to g, v to 0. See local_rest_terminal.
        [agTermT, agTermY] = local_rest_terminal(grid, aC, tEndMs, gRest);
        [vTermT,  vTermY]  = local_rest_terminal(grid, vC, tEndMs, 0);
        [agOnsT, agOnsY]   = local_ag_onset(grid, aC);

        % Where each row's drawn curve begins, after the pre-impact gate.
        % Reported by Diagnose beside the tail cut. A numeric triple rather
        % than a nested struct, for the reason local_interior_nan gives.
        preCutMs = [local_first_t(grid, zC), local_first_t(grid, vC), ...
                    local_first_t(grid, aC)];

        % Sub-g dip evidence, gathered here because it needs the per-replicate
        % matrix A, which Diagnose never sees. Numeric quadruple, not a nested
        % struct, for the reason local_interior_nan gives.
        dipInfo = local_dip_evidence(grid, aC, A, gRest);

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
            'nShortRec',    nShortRec, ...
            'interiorNaN',  local_interior_nan(zC, vC, aC), ...
            'agRaw',        aRaw, ...
            'tEndMs',       tEndMs, ...
            'agGateMs',     agGateMs, ...
            'preCutMs',     preCutMs, ...
            'dipInfo',      dipInfo, ...
            'gRest',        gRest, ...
            'agTermT',      agTermT, ...
            'agTermY',      agTermY, ...
            'vTermT',       vTermT, ...
            'vTermY',       vTermY, ...
            'agOnsetT',     agOnsT, ...
            'agOnsetY',     agOnsY, ...
            't_ms', grid, 'z', zC, 'v', vC, 'ag', aC);   %#ok<AGROW>

        % ── individual zero-drop overlays (only when explicitly enabled) ──
        if isZD && opt.ShowZeroDrop
            for j = 1:numel(loaded)
                sFull = loaded{j};
                depthRange = max(sFull.z, [], 'omitnan') - min(sFull.z, [], 'omitnan');
                if ~(depthRange > opt.IntrusionThreshCm), continue; end
                Tj = local_trial_on_grid(sFull, grid);   %#ok<ASGLU>
                zO = local_smooth_row(Tj.z, opt.GridMs, opt.SmoothMs.z, medStopMs, 'depth (overlay)');
                vO = local_smooth_row(Tj.v, opt.GridMs, opt.SmoothMs.v, medStopMs, 'v (overlay)');
                zO(past) = NaN;  vO(past) = NaN;
                O(end+1) = struct('model', K.model(rowsOk(j)), ...
                    'height', K.dropHeight_mm(rowsOk(j)), ...
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
        s = local_load_kin(K.kinPath(i));
        if isempty(s), continue; end

        idx = s.t_ms >= -opt.PreCapMs & s.t_ms <= s.t_stop_ms;
        if nnz(idx) < 2, continue; end
        t = s.t_ms(idx);

        zC = local_smooth_row(s.z(idx),  [], opt.SmoothMs.z,  s.t_stop_ms, 'depth', t);
        vC = local_smooth_row(s.v(idx),  [], opt.SmoothMs.v,  s.t_stop_ms, 'v',     t);
        aC = local_smooth_row(s.ag(idx), [], opt.SmoothMs.ag, s.t_stop_ms, 'a+g',   t, true);

        % NO a+g support gate here, and agGateMs is NaN to say so. A single
        % trial has support 1 at every point it covers, so applying the
        % ensemble gate (support >= MinReplicates, default 3) would blank the
        % entire a+g row in the one view whose whole purpose is showing what an
        % individual trial recorded.
        %
        % This view's endpoint is the trial's OWN last drawn sample, since it
        % ends at its own stop rather than at an ensemble one. Same terminal
        % rest segment, same rule.
        % No pre-impact gate either, for the same reason, so preCutMs is just
        % where each row's data happens to start. No sub-g dip statistics: the
        % evidence is a fraction OF REPLICATES, and there is one trial here.
        tEndMs = t(end);
        [agTermT, agTermY] = local_rest_terminal(t, aC, tEndMs, s.g_cm_s2);
        [vTermT,  vTermY]  = local_rest_terminal(t, vC, tEndMs, 0);
        [agOnsT,  agOnsY]  = local_ag_onset(t, aC);
        preCutMs = [local_first_t(t, zC), local_first_t(t, vC), ...
                    local_first_t(t, aC)];

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
            'nShortRec',    0, ...
            'interiorNaN',  local_interior_nan(zC, vC, aC), ...
            'agRaw',        s.ag(idx), ...
            'tEndMs',       tEndMs, ...
            'agGateMs',     NaN, ...
            'preCutMs',     preCutMs, ...
            'dipInfo',      [NaN NaN NaN 0], ...
            'gRest',        s.g_cm_s2, ...
            'agTermT',      agTermT, ...
            'agTermY',      agTermY, ...
            'vTermT',       vTermT, ...
            'vTermY',       vTermY, ...
            'agOnsetT',     agOnsT, ...
            'agOnsetY',     agOnsY, ...
            't_ms', t, 'z', zC, 'v', vC, 'ag', aC);
    end
end

function [tSeg, ySeg] = local_ag_onset(tAxis, ag)
%LOCAL_AG_ONSET  The free-fall boundary value, as one leading a+g segment.
%
%   During free fall the rod's acceleration IS gravity: a = g. Depth is
%   positive into the bed, so the net grain acceleration is
%
%       a + g = g - a = g - g = 0        identically, as t -> 0-
%
%   Zero is therefore not an assumption or a placeholder -- it is the exact
%   value of this quantity for every falling trial, at every height, fixed by
%   the definition in src/net_accel.m. It is a boundary state, not a
%   measurement, in the same sense as the rest state g the terminal segment
%   drops to: the pipeline computes no a outside [impact_index, stopFrame], so
%   neither endpoint appears in the stored trace, and both are known exactly
%   without being recorded.
%
%   Anchoring here is what makes the rise READ correctly. Without it the curve
%   simply materialises at ~g and climbs; with it the steep rise from 0 into
%   the first peak is the measured jump shown from its true starting state.
%
%   Implemented as a leading segment rather than by prepending to the ag
%   vector itself because t_ms is ONE axis shared by the depth, velocity and
%   a+g rows -- growing a+g alone would desynchronise it from the other two
%   and from agRaw, which the Diagnose onset assertion compares against.
%
%   The segment ends on the first computed post-impact sample, so no measured
%   value is displaced. When that sample sits at t = 0 (it does whenever
%   PreCapMs is a multiple of GridMs) the segment is vertical: the jump is
%   instantaneous at impact, which is exactly what the data says.
    tSeg = [];  ySeg = [];
    i = find(isfinite(ag), 1, 'first');
    if isempty(i), return; end
    if tAxis(i) < 0, return; end          % nothing sensible to anchor to
    tSeg = [0;      tAxis(i)];
    ySeg = [0;      ag(i)];
end

function [tSeg, ySeg] = local_rest_terminal(tAxis, y, tEndMs, yRest)
%LOCAL_REST_TERMINAL  The rest state, appended as one final segment.
%
%   Used by two rows, with the rest value each one takes at the stop:
%
%       a + g  ->  g     at rest the bed carries the rod's weight, a = 0
%       v      ->  0     at rest the rod is not moving
%
%   Both are MEASURED states, not extrapolations, and both are the one point
%   after the stop that is known rather than guessed. Neither row reaches its
%   rest value on its own: the pipeline computes no acceleration past
%   stopFrame, so a+g simply ends wherever it last had support; and the
%   velocity median ends SMALL-POSITIVE by construction, because at the
%   ensemble stop about half the replicates are already at rest at 0 and the
%   rest are still moving, so the median sits just above zero. Anchoring
%   removes that construction artefact rather than hiding a real residual
%   velocity.
%
%   The drop across the a+g segment carries physical content -- it is the
%   measured stopping discontinuity, how much net upward acceleration the
%   grains were still supplying at the instant the rod stopped. The drop
%   across the v segment does not; it is closing a median artefact.
%
%   tEndMs is the caller's UNIFIED endpoint -- the same t at which every panel
%   of that height ends. The segment therefore terminates exactly where the
%   others do, and nothing is drawn past it. When the last computed sample
%   already sits at tEndMs the segment is vertical, which is the sharp
%   terminal drop; when the row ran out earlier (the a+g support gate) it
%   slopes down to the same endpoint instead.
    tSeg = [];  ySeg = [];
    if ~isfinite(tEndMs) || ~isfinite(yRest), return; end
    i = find(isfinite(y), 1, 'last');
    if isempty(i), return; end
    if tAxis(i) > tEndMs, return; end        % nothing to append past the stop
    tSeg = [tAxis(i); tEndMs];
    ySeg = [y(i);     yRest];
end

function t = local_first_t(tAxis, y)
%LOCAL_FIRST_T  The time of a row's first finite sample; NaN if it has none.
%   The counterpart of local_last_t, for reporting where the pre-impact
%   support gate cut each row.
    t = NaN;
    i = find(isfinite(y), 1, 'first');
    if ~isempty(i), t = tAxis(i); end
end

function info = local_dip_evidence(grid, ag, A, g)
%LOCAL_DIP_EVIDENCE  Replicate-level evidence for a post-impact sub-g dip.
%
%   Returns [tMs, value, fractionAbove, n]:
%       tMs            time of the deepest post-impact point of the drawn a+g
%       value          a+g there
%       fractionAbove  fraction of the replicates CONTRIBUTING at that point
%                      whose own a > 0, i.e. whose own a+g < g
%       n              how many replicates contributed there
%   All NaN (n = 0) when the curve never dips below g after impact.
%
%   WHY THIS IS NOT A DEFECT TO BE FIXED. a + g = g - a, so a + g < g means
%   simply a > 0: the rod is still accelerating downward, because the grain
%   force on it is momentarily less than its own weight. That is ordinary
%   during penetration and must be preserved exactly -- nothing here clamps
%   or floors a+g.
%
%   What WOULD be a defect is the median dipping below g while most of its
%   replicates sit above it, which is the signature of an aggregate driven by
%   one or two outliers rather than by the population. So the fraction is
%   reported, and Diagnose flags it below 0.7. This is deliberately a
%   diagnostic and not a filter: it decides nothing about what is drawn.
%
%   Read from the PER-REPLICATE matrix, before smoothing or gating -- the
%   question is what the trials did, not what the display did to them.
    info = [NaN NaN NaN 0];
    idx = find(grid >= 0 & isfinite(ag));
    if isempty(idx) || ~isfinite(g), return; end

    [val, k] = min(ag(idx));
    i = idx(k);
    if ~(val < g), return; end               % no sub-g dip: nothing to report

    col = A(i, :);
    contrib = isfinite(col);
    n = sum(contrib);
    if n == 0, return; end

    info = [grid(i), val, sum(col(contrib) < g) / n, n];
end

function C = local_empty_curve()
    C = struct('model',{}, 'height',{}, 'medV0',{}, 'isZeroDrop',{}, ...
               'nMax',{}, 'medStopMs',{}, 'tStopMinMs',{}, 'tStopMaxMs',{}, ...
               'supportAtEnd',{}, 'supportMeasAtEnd',{}, 'oldRuleEndMs',{}, ...
               'nShortRec',{}, 'interiorNaN',{}, 'agRaw',{}, ...
               'tEndMs',{}, 'agGateMs',{}, 'preCutMs',{}, 'dipInfo',{}, ...
               'gRest',{}, 'agTermT',{}, 'agTermY',{}, ...
               'vTermT',{}, 'vTermY',{}, 'agOnsetT',{}, 'agOnsetY',{}, ...
               't_ms',{}, 'z',{}, 'v',{}, 'ag',{});
end

function O = local_empty_overlay()
    O = struct('model',{}, 'height',{}, 'medV0',{}, 't_ms',{}, 'z',{}, 'v',{});
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

function counts = local_interior_nan(zC, vC, aC)
%LOCAL_INTERIOR_NAN  NaNs strictly between a curve's first and last finite
%   sample, per row, returned as [z v ag].
%
%   Under the rest-extension every ensemble curve must be continuous from its
%   onset to the ensemble endpoint, so a nonzero count is a defect rather than
%   a property of the data -- Diagnose asserts on it.
%
%   a+g is counted the same way but spans only [impact, stop] by construction;
%   its "onset" is its first computed sample, not the start of the grid.
%
%   Returned as a numeric triple rather than a struct: struct() field values
%   that are themselves structs are a needless hazard inside the struct(...)
%   calls that build C.
    counts = [local_count_interior(zC), ...
              local_count_interior(vC), ...
              local_count_interior(aC)];
end

function n = local_count_interior(y)
    n = 0;
    f = find(isfinite(y), 1, 'first');
    l = find(isfinite(y), 1, 'last');
    if isempty(f) || l <= f, return; end
    n = sum(~isfinite(y(f:l)));
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

function y = local_smooth_row(yRaw, gridMs, windowMs, medianStopMs, rowLabel, tMs, edgePreserve)
%LOCAL_SMOOTH_ROW  movmean AFTER aggregation, then re-mask to the original
%   support so smoothing cannot extend a curve past where the data ends.
%
%   movmean's 'shrink' endpoint rule truncates the window at the boundary --
%   it does NOT leave the boundary sample alone. At the first sample a window
%   of w still averages about w/2 samples. For depth and velocity that is
%   harmless: their support starts PreCapMs before impact, where the curve is
%   flat and the window is full well before anything interesting happens.
%
%   For a + g it is not harmless, and 'edgePreserve' exists for it. That row's
%   support starts abruptly AT impact -- kd_kinematics masks a outside
%   [impact_index, stopFrame] -- and the first millisecond is the steepest part
%   of the whole trace. A shrunk-but-still-multi-sample window at that boundary
%   averages the measured onset together with the rise above it and reports the
%   mid-rise value as the starting point.
%
%   With edgePreserve the window's half-width is the distance to the nearest
%   end of the finite support, capped at the full half-width: exactly one
%   sample at each end, growing symmetrically to the full window inside. The
%   first and last samples are therefore returned unchanged, and because the
%   window stays symmetric nothing is phase-shifted.
    if nargin < 7 || isempty(edgePreserve), edgePreserve = false; end
    if isempty(gridMs)
        if nargin < 6 || numel(tMs) < 2, y = yRaw; return; end
        gridMs = median(diff(tMs), 'omitnan');
    end
    if ~(gridMs > 0) || ~isfinite(gridMs), y = yRaw; return; end

    w = max(1, round(windowMs / gridMs));
    y = movmean(yRaw, w, 'omitnan');
    y(isnan(yRaw)) = NaN;
    if edgePreserve && w > 1
        y = local_ramp_edges(yRaw, y, w);
    end

    if isfinite(medianStopMs) && medianStopMs > 0 && windowMs > 0.10*medianStopMs
        warning('fig_kinematics:overSmooth', ...
            ['%s smoothing window (%.2f ms) exceeds 10%% of the median stop ' ...
             'time (%.2f ms). The curve may be over-smoothed.'], ...
            rowLabel, windowMs, medianStopMs);
    end
end

function y = local_ramp_edges(yRaw, y, w)
%LOCAL_RAMP_EDGES  Re-average the samples near each end of a row's finite
%   support with a symmetric window that shrinks to a single sample at the edge.
%
%   Half-width = distance to the nearer end of the support, capped at floor(w/2)
%   -- so the edge sample is itself, the next one averages three, and so on
%   until the full movmean window applies. Written as one pass over the support
%   rather than two edge loops so that a support shorter than the window cannot
%   have one edge's ramp overwrite the other's.
%
%   NaN-aware in the same sense as movmean(...,'omitnan'): a window containing
%   only NaN yields NaN, which the caller's mask would restore anyway.
    n = numel(yRaw);
    f = find(isfinite(yRaw), 1, 'first');
    l = find(isfinite(yRaw), 1, 'last');
    if isempty(f) || l <= f, return; end
    half = floor(w/2);
    if half < 1, return; end

    for i = f:l
        k = min(i - f, l - i);
        if k >= half, continue; end      % full window already applies here
        lo  = max(f, i - k);
        hi  = min(min(l, n), i + k);
        seg = yRaw(lo:hi);
        seg = seg(isfinite(seg));
        if isempty(seg), y(i) = NaN; else, y(i) = mean(seg); end
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  DIAGNOSTICS
% ═════════════════════════════════════════════════════════════════════════
function local_diagnose(C, opt)
%LOCAL_DIAGNOSE  Per-height numbers behind the ensemble, and the a+g onset.
%
%   Documents what the render alone cannot show, and asserts on three things
%   that must hold:
%
%     * how far apart the replicate stop times are -- the mechanism the retired
%       MinReplicates rule turned into survivor bias
%     * CONTINUITY: no interior NaN in any curve, after the rest fill
%     * ENDPOINTS: depth, velocity and the a+g terminal segment all end at the
%       same t for each (model, height)
%     * ONSET: display smoothing leaves the first a+g sample exactly as the
%       aggregate produced it

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

        % Where the a+g support gate cut each height's computed curve, and how
        % much of the interval to the unified stop the terminal segment then
        % spans. A large gap means a long thin-support tail was removed -- the
        % wandering few-replicate median that dipped below g.
        fprintf('\n    support gates (support >= %d), per height:\n', ...
                opt.MinReplicates);
        fprintf(['      %5s  %-28s  %s\n'], 'h(mm)', ...
                'pre-impact cut  z / v / a+g', 'a+g tail cut -> unified stop');
        for q = sel
            pc = C(q).preCutMs;
            if isfinite(C(q).agGateMs)
                tailStr = sprintf('%7.2f -> %7.2f ms (gap %5.2f)', ...
                    C(q).agGateMs, C(q).tEndMs, C(q).tEndMs - C(q).agGateMs);
            else
                tailStr = 'no gate applied';
            end
            fprintf('      %5g  %7.2f /%7.2f /%7.2f      %s\n', ...
                    C(q).height, pc(1), pc(2), pc(3), tailStr);
        end
        fprintf(['      pre-impact cut = where each row starts once thin-support\n' ...
                 '      leading samples are dropped; the rows need not agree.\n']);

        % Sub-g dips: evidence, not a correction. a + g < g means a > 0, the
        % grain force momentarily below the foot's weight, which is ordinary
        % during penetration. What would NOT be ordinary is the median dipping
        % below g while most of its own replicates stay above it -- an
        % aggregate driven by one or two outliers. Hence the fraction.
        anyDip = false;
        for q = sel
            d = C(q).dipInfo;
            if ~isfinite(d(2)), continue; end
            if ~anyDip
                fprintf('\n    sub-g dips (a+g < g after impact; a > 0 is physical):\n');
                anyDip = true;
            end
            flag = '';
            if isfinite(d(3)) && d(3) < 0.70
                flag = '   <-- FLAG: below 0.70, dip may be outlier-driven';
            end
            fprintf(['      h = %4g mm : min a+g %8.1f at %7.2f ms; ' ...
                     '%.0f%% of %d contributing replicates have a > 0%s\n'], ...
                    C(q).height, d(2), d(1), 100*d(3), d(4), flag);
        end
        if ~anyDip
            fprintf('\n    sub-g dips: none (no curve falls below g after impact)\n');
        end

        % Continuity: after rest-filling to the grid end, every ensemble curve
        % must be gap-free from its onset to the endpoint. A nonzero count is
        % the mid-curve break this fill was introduced to remove.
        bad = false;
        for q = sel
            c = C(q).interiorNaN;
            if any(c > 0)
                bad = true;
                fprintf(['    CONTINUITY FAIL h = %g mm: interior NaNs ' ...
                         'z=%d v=%d a+g=%d\n'], C(q).height, c(1), c(2), c(3));
            end
        end
        if ~bad
            fprintf('    continuity: OK -- no interior NaNs in any curve (z, v, a+g)\n');
        end

        % Endpoint identity: depth, velocity and the a+g terminal segment must
        % all finish at the SAME t for a given (model, height). One stop per
        % height per geometry, shared across that figure's three panels -- if
        % these ever diverged, the a+g panel would be telling a different story
        % about when the rod stopped than the two above it.
        %
        % ENDPOINTS ONLY, deliberately. The three panels do NOT share a start:
        % depth and velocity begin at -PreCapMs, while a+g begins at its
        % free-fall anchor at t = 0 because the pipeline computes no
        % acceleration before impact. That asymmetry is the design, so this
        % must never be widened into a startpoint check.
        badEnd = false;
        for q = sel
            % Where each row is actually LAST DRAWN: for the two anchored rows
            % that is the end of their terminal segment, not the end of the
            % computed series, which the a+g support gate may cut earlier.
            tz = local_last_t(C(q).t_ms, C(q).z);
            if ~isempty(C(q).vTermT)
                tv = C(q).vTermT(end);
            else
                tv = local_last_t(C(q).t_ms, C(q).v);
            end
            if ~isempty(C(q).agTermT)
                ta = C(q).agTermT(end);
            else
                ta = local_last_t(C(q).t_ms, C(q).ag);
            end
            % Half a sample of whatever axis this curve is on -- the ensemble
            % grid under Style 'mean', the trial's own frames under 'trials'.
            % They should agree exactly; the tolerance only stops a rounding
            % difference being reported as a failure.
            tol = 0.5 * median(diff(C(q).t_ms), 'omitnan');
            if ~(tol > 0), tol = 0.5 * opt.GridMs; end
            if ~(isfinite(tz) && isfinite(tv) && isfinite(ta)) || ...
                    max(abs([tz tv ta] - C(q).tEndMs)) > tol
                badEnd = true;
                fprintf(['    ENDPOINT FAIL h = %g mm: z ends %.3f, v ends ' ...
                         '%.3f, a+g ends %.3f, unified %.3f ms\n'], ...
                        C(q).height, tz, tv, ta, C(q).tEndMs);
            end
        end
        if ~badEnd
            fprintf(['    endpoints: OK -- z, v and a+g all end at the same ' ...
                     'unified stop per height\n']);
        end

        nShort = sum([C(sel).nShortRec]);
        if nShort > 0
            fprintf(['    NOTE: %d trial(s) had recordings ending before their ' ...
                     'own t_stop and were\n          treated as at rest from ' ...
                     'the last recorded sample -- an approximation.\n'], nShort);
        end

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
            fprintf(['      (smoothing window %g ms, edge-preserving: the ' ...
                     'window''s half-width is\n       the distance to the ' ...
                     'edge of the support, so the first sample above is\n' ...
                     '       returned unchanged. movmean''s own ''shrink'' ' ...
                     'rule does NOT do this --\n       it truncates the ' ...
                     'window but still averages ~w/2 samples, which lifted\n' ...
                     '       the onset off g and into the rise.)\n'], ...
                    opt.SmoothMs.ag);
            d0 = abs(C(q).ag(gi) - C(q).agRaw(gi));
            if d0 > 1e-9
                fprintf(['      ONSET FAIL: smoothed first sample differs ' ...
                         'from the aggregate by %.3g\n'], d0);
            end
        end
    end
    fprintf('\n');
end

function s = local_fmt(v)
    s = strjoin(arrayfun(@(x) sprintf('%9.2f', x), v(:).', ...
                         'UniformOutput', false), ' ');
end

function t = local_last_t(tAxis, y)
%LOCAL_LAST_T  The time of a row's last finite sample; NaN if it has none.
    t = NaN;
    i = find(isfinite(y), 1, 'last');
    if ~isempty(i), t = tAxis(i); end
end

% ═════════════════════════════════════════════════════════════════════════
%  LIMITS AND COLOUR
% ═════════════════════════════════════════════════════════════════════════
function LIMS = local_shared_limits(C, O, opt)
%LOCAL_SHARED_LIMITS  One set of x/y limits per row, shared across every model
%   and every figure, so the three geometry figures are directly comparable.
%
%   Depth and velocity carry a 5% margin beyond the data extent, so no curve or
%   stop marker sits on a panel edge. The depth panel needs this most: its
%   markers sit exactly at each curve's maximum and were being clipped by the
%   frame.
%
%   a+g is different and deliberately so: its floor is a hard 0 with no
%   downward padding, because 0 is the free-fall boundary value the onset
%   segments anchor to rather than an arbitrary data minimum. Only the top gets
%   the 5% headroom.
%
%   The x-limit runs from -PreCapMs to the LATEST ensemble endpoint across all
%   geometries -- the last point actually drawn, not the end of the internal
%   grid, which extends past every curve's termination to the longest replicate
%   stop and would leave a wide empty margin on the right.

    MARGIN = 0.05;

    lastDrawn = -inf;
    for k = 1:numel(C)
        idx = find(isfinite(C(k).z) | isfinite(C(k).v) | isfinite(C(k).ag), 1, 'last');
        if ~isempty(idx), lastDrawn = max(lastDrawn, C(k).t_ms(idx)); end
        % The a+g terminal segment ends at the same unified stop as z and v, so
        % this should never extend the range -- included so that if it ever did,
        % the segment would be inside the frame rather than clipped off it.
        if ~isempty(C(k).agTermT)
            lastDrawn = max(lastDrawn, max(C(k).agTermT));
        end
    end
    for k = 1:numel(O)
        idx = find(isfinite(O(k).z) | isfinite(O(k).v), 1, 'last');
        if ~isempty(idx), lastDrawn = max(lastDrawn, O(k).t_ms(idx)); end
    end
    if ~isfinite(lastDrawn), lastDrawn = 0; end

    tLo = -opt.PreCapMs;
    XL  = [tLo, lastDrawn + MARGIN*(lastDrawn - tLo)];

    zAll = vertcat(C.z);
    if ~isempty(O), zAll = [zAll; vertcat(O.z)]; end
    YL_z = local_pad_ylim(zAll, MARGIN);

    vAll = [vertcat(C.v); vertcat(C.vTermY)];
    if ~isempty(O), vAll = [vAll; vertcat(O.v)]; end
    YL_v = local_pad_ylim(vAll, MARGIN);
    YL_v(1) = 0;                 % velocity floor is zero: nothing below rest

    % Every drawn a+g value: the curves, the terminal segments reaching down to
    % g, and the onset segments reaching down to 0. All three are part of the
    % same line, so all three set the range.
    agAll = [vertcat(C.ag); vertcat(C.agTermY); vertcat(C.agOnsetY)];

    % ONE shared a+g range across all three geometry figures. The floor is a
    % hard 0 -- not padded downward and never raised -- because 0 is the
    % free-fall boundary value the onset segments anchor to, and lifting the
    % floor would cut them off. The top is 1.05x the global maximum over every
    % plotted curve, so the tallest peak in any geometry clears the frame by
    % 5% and the same a+g reads at the same height in all three figures.
    agMax = max(agAll(isfinite(agAll)));
    if isempty(agMax) || ~(agMax > 0), agMax = 1; end
    YL_ag_lin = [0, 1.05 * agMax];
    YL_ag_log = local_log_ylim(agAll, opt.AccelYMin);

    % The rest baseline drawn behind the a+g data. Taken from the curves so it
    % is the same g those traces were computed with (net_accel's second
    % output), not a constant restated here.
    gRef = NaN;
    gAll = [C.gRest];
    gAll = gAll(isfinite(gAll));
    if ~isempty(gAll), gRef = median(gAll); end

    LIMS = struct('t', XL, 'z', YL_z, 'v', YL_v, ...
                  'ag_lin', YL_ag_lin, 'ag_log', YL_ag_log, 'gRef', gRef);
end

function yl = local_pad_ylim(vals, margin)
    vals = vals(isfinite(vals));
    if isempty(vals), yl = [-1 1]; return; end
    lo = min(vals); hi = max(vals); r = hi - lo;
    if r == 0, r = max(abs(hi), 1); end
    yl = [lo - margin*r, hi + margin*r];
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

function REF = local_global_color_ref(KAll, opt)
%LOCAL_GLOBAL_COLOR_REF  The one colour reference every figure maps through.
%
%   Returns a struct with everything both encodings need:
%       .clim    [lo hi] median-v0 range, for ColorBy 'v0'
%       .ladder  sorted unique drop heights, for ColorBy 'height'
%       .nGroups how many (model, height) groups fed the two
%
%   THE BUG THIS FIXES (the 'v0' half). The range used to be taken from the
%   curves that had actually been built, and those exist only for the models
%   in 'Models' (fig_kinematics filters K by it right after loading). So the
%   range was a property of the invocation, not of the dataset: drawing all
%   three geometries in one call gave one mapping, and regenerating a single
%   geometry on its own -- the natural way to redraw one figure -- gave a
%   different one. Same physical v0, different colour, and the three figures
%   stopped being comparable. Nothing was wrong inside the per-model loop; the
%   range handed to it was already wrong.
%
%   The fix, and the reason BOTH references are built here: they come from the
%   FULL loaded set, before that filter, so neither can depend on which models
%   are being drawn.
%
%   Reproduces each curve's medV0 and height without loading a single trace:
%   group by (model, height), take the median v0 over the group, keep groups
%   with at least MinReplicates trials. Zero-drop rows are excluded outright,
%   whether or not they are drawn -- their v0 is 0 by quarantine rather than by
%   measurement (load_kinematics_set: no fall, so no impact speed), and
%   including it would anchor the v0 scale at zero and put a phantom rung at
%   the bottom of the height ladder.
%
%   The group set is a SUPERSET of the curves that end up drawn: a group can
%   clear MinReplicates here and still lose its curve later if a trial proves
%   unreadable or has no finite t_stop. For 'v0' that can only widen the range
%   slightly; for 'height' it can only add a rung no curve occupies, which
%   shifts nothing, since a height's colour is its rank on the ladder and the
%   ladder is the same in every figure either way. Both are stable across
%   invocations, which is the property being bought.
    REF = struct('clim', [0 1], 'ladder', zeros(0,1), 'nGroups', 0);

    R = KAll(~KAll.isZeroDrop, :);
    if isempty(R), return; end

    g   = findgroups(R.model, R.dropHeight_mm);
    med = splitapply(@(x) median(x, 'omitnan'), R.v0_cm_s, g);
    hgt = splitapply(@(x) x(1),                 R.dropHeight_mm, g);
    cnt = splitapply(@numel,                    R.v0_cm_s, g);

    big = cnt >= opt.MinReplicates;
    medBig = med(big & isfinite(med));
    hgtBig = hgt(big & isfinite(hgt));

    % The ladder is the set of nominal heights, shared by every geometry: one
    % fixed height -> colour map, so the three figures step through the same
    % colours in the same order by construction. Heights the models do not
    % share still occupy their own rung -- the ladder is the union, not the
    % intersection, so a height present in only one geometry does not shift
    % the colours of the heights around it in the others.
    REF.ladder = unique(hgtBig(:));

    if isempty(medBig)
        % No group is big enough to be plotted; fall back to the raw spread so
        % a 'trials' run still gets a sane scale.
        medBig = R.v0_cm_s(isfinite(R.v0_cm_s));
        if isempty(REF.ladder)
            REF.ladder = unique(R.dropHeight_mm(isfinite(R.dropHeight_mm)));
        end
    end
    if isempty(medBig), return; end

    REF.nGroups = numel(medBig);
    REF.clim = [min(medBig), max(medBig)];
    if ~(diff(REF.clim) > 0)
        REF.clim = REF.clim(1) + [0, max(1, abs(REF.clim(1)))];
    end
end

function local_print_color_check(label, C, REF, CMAP, opt)
%LOCAL_PRINT_COLOR_CHECK  The colour mapping actually used, for this figure's
%   extreme curves, in whichever mode is active.
%
%   The check the shared-reference claim has to survive: run the three
%   geometries and compare.
%
%     'height'  the SAME height must print the SAME index in every figure.
%               Because the index is a rank on a fixed shared ladder, the
%               three figures step through identical colours by construction,
%               and any difference here means the ladder was rebuilt.
%     'v0'      a v0 appearing in two figures must print the same index, and
%               every line must quote the same clim.
    inCol = find([C.model] == label);
    if isempty(inCol), return; end
    n = size(CMAP, 1);

    if strcmp(opt.ColorBy, 'height')
        h = [C(inCol).height];
        h = h(isfinite(h));
        if isempty(h), return; end
        fprintf(['  %-8s  low  h %6.4g mm -> idx %3d   high h %6.4g mm -> ' ...
                 'idx %3d   (ladder %d rungs, %g - %g mm)\n'], label, ...
                min(h), local_height_index(min(h), REF.ladder, n), ...
                max(h), local_height_index(max(h), REF.ladder, n), ...
                numel(REF.ladder), min(REF.ladder), max(REF.ladder));
    else
        v = [C(inCol).medV0];
        v = v(isfinite(v));
        if isempty(v), return; end
        fprintf(['  %-8s  low  v0 %7.2f -> idx %3d   high v0 %7.2f -> idx ' ...
                 '%3d   (clim %.2f - %.2f)\n'], label, ...
                min(v), local_v0_index(min(v), REF.clim, n), ...
                max(v), local_v0_index(max(v), REF.clim, n), ...
                REF.clim(1), REF.clim(2));
    end
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

function idx = local_v0_index(v0, clim, nLevels)
%LOCAL_V0_INDEX  v0 -> colormap row. The ONE place the mapping is defined, so
%   the verification printout cannot drift from what the drawing actually does.
    f = (v0 - clim(1)) / (clim(2) - clim(1));
    if ~isfinite(f), f = 0; end
    f = min(max(f, 0), 1);
    idx = 1 + round(f * (nLevels - 1));
end

function idx = local_height_index(h, ladder, nLevels)
%LOCAL_HEIGHT_INDEX  Nominal drop height -> colormap row, by RANK on the
%   shared ladder rather than by value.
%
%   Rank, not linear interpolation of the height itself, is what makes the
%   three geometries step through identical colours: every figure draws the
%   same rungs in the same order, so rung r gets the same colour everywhere
%   regardless of which heights a particular geometry happens to have. It also
%   spreads the colormap evenly over the heights actually used instead of
%   bunching them wherever the ladder's spacing is tight.
%
%   Coldest at the bottom rung, warmest at the top: 25 mm -> index 1,
%   365 mm -> index nLevels, on the campaign ladder.
%
%   A height not on the ladder -- a zero-drop overlay at h = 0, which the
%   ladder deliberately excludes -- falls to the nearest rung, so it takes the
%   coldest colour rather than an error or a blank.
    idx = 1;
    if isempty(ladder) || ~isfinite(h), return; end
    n = numel(ladder);
    if n < 2, return; end
    [~, r] = min(abs(ladder(:) - h));        % nearest rung; exact when on it
    idx = 1 + round((r - 1) / (n - 1) * (nLevels - 1));
end

function col = local_curve_color(crv, REF, CMAP, opt)
%LOCAL_CURVE_COLOR  One curve's colour, in whichever encoding is active.
%   The single dispatch point, so the drawing and the verification printout
%   cannot disagree about what a curve's colour is.
    n = size(CMAP, 1);
    if strcmp(opt.ColorBy, 'height')
        col = CMAP(local_height_index(crv.height, REF.ladder, n), :);
    else
        col = CMAP(local_v0_index(crv.medV0, REF.clim, n), :);
    end
end

% ═════════════════════════════════════════════════════════════════════════
%  DRAWING
% ═════════════════════════════════════════════════════════════════════════
function fig = local_draw_one_model(model, C, O, REF, CMAP, LIMS, opt)
%LOCAL_DRAW_ONE_MODEL  Three stacked panels on a square canvas.
%   Every panel keeps full numeric tick labels on BOTH axes -- the x labels are
%   not blanked on the upper two. Repeating them costs a little space and
%   removes any doubt about which panel a reader is looking at.
    fig = figure('Color','w','Units','inches','Position',[1 1 7.0 7.0], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl = tiledlayout(fig, 3, 1, 'Padding','compact', 'TileSpacing','compact');

    ax1 = nexttile(tl); local_draw_panel(ax1, 'z',  model, C, O, REF, CMAP, LIMS, opt);
    ax2 = nexttile(tl); local_draw_panel(ax2, 'v',  model, C, O, REF, CMAP, LIMS, opt);
    ax3 = nexttile(tl); local_draw_panel(ax3, 'ag', model, C, O, REF, CMAP, LIMS, opt);

    % The geometry is named on the figure, not only in the filename: these are
    % three separate files and a reader should not have to check the path.
    title(ax1, model, 'Interpreter', 'none');

    xlabel(ax1, 't - t_impact (ms)', 'Interpreter', 'none');
    xlabel(ax2, 't - t_impact (ms)', 'Interpreter', 'none');
    xlabel(ax3, 't - t_impact (ms)', 'Interpreter', 'none');
end

function fig = local_draw_grid3x3(models, C, O, REF, CMAP, LIMS, opt)
    nCol = numel(models);
    fig = figure('Color','w','Units','inches','Position',[1 1 7.0 6.4], ...
                 'Visible', local_tern(opt.Show,'on','off'));
    tl  = tiledlayout(fig, 3, nCol, 'Padding','compact', 'TileSpacing','compact');

    ROWS = {'z','v','ag'};
    for r = 1:3
        for c = 1:nCol
            ax = nexttile(tl, (r-1)*nCol + c);
            local_draw_panel(ax, ROWS{r}, models(c), C, O, REF, CMAP, LIMS, opt);
            % Full numeric tick labels on both axes of every panel: nothing is
            % blanked, even where a row or column repeats them.
            if r == 3
                xlabel(ax, 't - t_impact (ms)', 'Interpreter', 'none');
            end
            if r == 1
                title(ax, models(c), 'Interpreter', 'none');
            end
        end
    end
end

function local_draw_panel(ax, rowKey, model, C, O, REF, CMAP, LIMS, opt)
%LOCAL_DRAW_PANEL  One (row, model) panel -- shared by both layouts and both
%   styles, since C/O are already in a common shape.

    % House style shared with depth_scaling.m (grid on, box on, MATLAB default
    % ticks and fonts) -- see src/apply_fig_style.m.
    apply_fig_style(ax);

    % Drawn FIRST so they sit behind the data.
    % Impact is the reference every panel is drawn against, so its line is
    % black rather than a light grey that reads as chart furniture.
    xline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 0.9, ...
          'HandleVisibility', 'off');
    if strcmp(rowKey, 'z')
        % Solid thin line at depth = 0: the bed surface, a physical datum.
        yline(ax, 0, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.5, ...
              'HandleVisibility', 'off');
    end
    if strcmp(rowKey, 'ag') && isfinite(LIMS.gRef)
        % The rest baseline, a + g = g. Dotted, light and unlabelled: it is
        % there so the terminal drops can be seen landing on a marked value,
        % not to be read off. Deliberately NOT the removed dashed g line of
        % the earlier draft, which sat on top of the data at full weight.
        yline(ax, LIMS.gRef, ':', 'Color', [0.60 0.60 0.60], 'LineWidth', 0.75, ...
              'HandleVisibility', 'off');
    end

    mkMap = containers.Map({'Default','Tight','Wide'}, {'o','s','^'});
    % Deliberately different from depth_scaling.m's Tight=^/Wide=s convention:
    % this figure's marker mapping is specified independently in its own brief.

    inCol = find([C.model] == model);
    for k = inCol
        y = C(k).(rowKey);
        if isempty(y) || all(isnan(y)), continue; end
        if strcmp(rowKey, 'ag') && strcmp(opt.AccelScale, 'log')
            % A log axis cannot render non-positive values; mask them so the
            % line breaks rather than silently vanishing.
            y(y <= 0) = NaN;
            if all(isnan(y)), continue; end
        end
        col = local_curve_color(C(k), REF, CMAP, opt);
        plot(ax, C(k).t_ms, y, '-', 'LineWidth', 2.0, 'Color', col);

        if strcmp(rowKey, 'ag') && ~isempty(C(k).agTermT)
            % The measured rest state: a = 0 at rest, so a + g = g. Same colour
            % and weight as the curve because it IS the curve's last segment --
            % it ends at the same unified stop as the depth and velocity
            % panels, and nothing is drawn past it. See local_rest_terminal.
            tSeg = C(k).agTermT;  ySeg = C(k).agTermY;
            if strcmp(opt.AccelScale, 'log'), ySeg(ySeg <= 0) = NaN; end
            plot(ax, tSeg, ySeg, '-', 'LineWidth', 2.0, 'Color', col, ...
                 'HandleVisibility', 'off');
        end

        if strcmp(rowKey, 'v') && ~isempty(C(k).vTermT)
            % The rest state for velocity: 0 at the same unified stop. Exactly
            % parallel to the a+g segment above, and for a narrower reason --
            % the median ends small-positive by construction, because half the
            % replicates are already at rest at 0 there and half are still
            % moving. This closes that construction artefact; it is not
            % asserting the rod decelerated any faster than measured.
            plot(ax, C(k).vTermT, C(k).vTermY, '-', 'LineWidth', 2.0, ...
                 'Color', col, 'HandleVisibility', 'off');
        end

        if strcmp(rowKey, 'ag') && ~isempty(C(k).agOnsetT)
            % The free-fall boundary state: a = g while falling, so a + g = 0
            % exactly. Same colour and weight -- it is the curve's first
            % segment. See local_ag_onset. Masked out on a log axis, which
            % cannot draw zero at all.
            tSeg = C(k).agOnsetT;  ySeg = C(k).agOnsetY;
            if strcmp(opt.AccelScale, 'log'), ySeg(ySeg <= 0) = NaN; end
            if any(isfinite(ySeg))
                plot(ax, tSeg, ySeg, '-', 'LineWidth', 2.0, 'Color', col, ...
                     'HandleVisibility', 'off');
            end
        end

        if strcmp(rowKey, 'z')
            % Stop marker at the unified endpoint: the height's median t_stop.
            last = find(~isnan(y), 1, 'last');
            if ~isempty(last)
                mk = 'd';
                if isKey(mkMap, char(model)), mk = mkMap(char(model)); end
                plot(ax, C(k).t_ms(last), y(last), mk, 'MarkerSize', 7, ...
                     'MarkerFaceColor', col, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
                     'LineWidth', 0.5);
            end
        end
    end

    if ismember(rowKey, {'z','v'}) && ~isempty(O)
        for k = find([O.model] == model)
            y = O(k).(rowKey);
            if isempty(y) || all(isnan(y)), continue; end
            col = local_curve_color(O(k), REF, CMAP, opt);
            plot(ax, O(k).t_ms, y, '-', 'LineWidth', 0.75, 'Color', col);
        end
    end

    isLogAg = strcmp(rowKey, 'ag') && strcmp(opt.AccelScale, 'log');
    switch rowKey
        case 'ag'
            if isLogAg
                set(ax, 'YScale', 'log');
                ylim(ax, LIMS.ag_log);
            else
                set(ax, 'YScale', 'linear');
                ylim(ax, LIMS.ag_lin);   % 0 to 1.05x global max; MATLAB adds
            end                          % the x10^4 multiplier on the labels
            ylabel(ax, 'a + g (cm/s^2)', 'Interpreter', 'none');
        case 'v'
            ylim(ax, LIMS.v);
            ylabel(ax, 'v (cm/s)', 'Interpreter', 'none');
        otherwise
            ylim(ax, LIMS.z);
            ylabel(ax, 'depth (cm)', 'Interpreter', 'none');
    end
    xlim(ax, LIMS.t);

    % ── ticks ────────────────────────────────────────────────────────────
    % Explicit major ticks on both axes, every one labelled (no XTickLabel is
    % set, so MATLAB labels them all). No minor ticks: the shared house style
    % in src/apply_fig_style.m leaves them at MATLAB's default off, and turning
    % them on here would desynchronise this figure from depth_scaling.
    % A grid3x3 panel is a third the width of a per-model one, so it takes a
    % third of the labels before the time axis starts to collide.
    maxXT = 14;
    if strcmp(opt.Layout, 'grid3x3'), maxXT = 8; end
    xt = local_time_ticks(LIMS.t, maxXT);
    if numel(xt) >= 2, set(ax, 'XTick', xt); end
    if ~isLogAg
        % Log ticks are decade powers and MATLAB already picks them well.
        yt = local_round_ticks(ylim(ax), 6);
        if numel(yt) >= 2, set(ax, 'YTick', yt); end
    end
end

function ticks = local_time_ticks(xl, maxTicks)
%LOCAL_TIME_TICKS  Major time ticks on a round millisecond step.
%   Prefers 5 ms and steps up through round values only when 5 would crowd the
%   axis. The pre-impact window is 10 ms by default, so a 5 ms step puts ticks
%   at -10 and -5 and no closer: the crowding to avoid is on the long
%   post-impact side, not before impact. Which step wins therefore depends on
%   how far the slowest geometry's curves run, not on the pre-window.
    ticks = local_step_ticks(xl, [5 10 20 25 50 100 200], maxTicks);
end

function ticks = local_round_ticks(yl, maxTicks)
%LOCAL_ROUND_TICKS  Major ticks on a round step of 1, 2, 2.5 or 5 x 10^n.
%   The standard nice-number rule: take the smallest round step that keeps the
%   tick count at or under maxTicks.
    ticks = [];
    r = yl(2) - yl(1);
    if ~(r > 0) || ~all(isfinite(yl)), return; end
    raw = r / max(2, maxTicks - 1);
    mag = 10^floor(log10(raw));
    step = 10*mag;
    for m = [1 2 2.5 5 10]
        if m*mag >= raw, step = m*mag; break; end
    end
    ticks = (ceil(yl(1)/step) : floor(yl(2)/step)) * step;
end

function ticks = local_step_ticks(lims, candidates, maxTicks)
%LOCAL_STEP_TICKS  First candidate step giving at most maxTicks ticks.
    ticks = [];
    if ~all(isfinite(lims)) || ~(lims(2) > lims(1)), return; end
    step = candidates(end);
    for s = candidates
        n = floor(lims(2)/s) - ceil(lims(1)/s) + 1;
        if n >= 2 && n <= maxTicks, step = s; break; end
    end
    ticks = (ceil(lims(1)/step) : floor(lims(2)/step)) * step;
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
