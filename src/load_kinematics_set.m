function K = load_kinematics_set(root, varargin)
% LOAD_KINEMATICS_SET  Load saved kinematics into one table for analysis.
%
%   K = load_kinematics_set(root)
%       'model'      ''|'Default'|'Tight'|'Wide'   geometry from tag suffix
%       'condition'  ''|'GB/full'|...              material/container
%       'exclude'    get_manual_exclusions()        trialTags to drop
%       'series'     false                         also return z/v/a+g arrays
%
%   Reads <container>/kinematics/*_kin_scalars.csv (and *_kin.mat when
%   'series' is true). Returns one row per trial with, at minimum:
%     trialTag, model, condition, dropHeight_mm, fps, impactFrame, stopFrame,
%     v0_cm_s, d_final_cm, a_stop_cm_s2, impactDistPx, bedX
%   plus the *_meas companions described under ZERO-DROP QUARANTINE below,
%   and dropHeight_asRecorded / relabeled described under HEIGHT RELABELING.
%
%   HEIGHT RELABELING. dropHeight_mm is corrected at read time from
%   get_relabel_map, for the handful of trials whose recorded height label is
%   known to be wrong. Nothing on disk is modified -- see get_relabel_map for
%   why. Two columns record what happened:
%
%     dropHeight_asRecorded  the label as it came off the CSV
%     relabeled              true where dropHeight_mm came from the map
%
%   dropHeight_mm is the corrected height and is what every height axis, bin
%   and fit should use. dropHeight_asRecorded is provenance, not a measurement.
%
%   ZERO-DROP QUARANTINE. A d0 trial is released FROM CONTACT: no fall, no
%   impact, no deceleration. kd_kinematics is a dynamic-impact pipeline and is
%   correct for its domain; h = 0 is outside it. Everything it derives from the
%   impact/stop event is therefore UNDEFINED on those rows, not merely noisy:
%
%     v0_cm_s        set to 0        (protocol: no fall, so no impact speed)
%     d_final_cm     set to NaN
%     t_stop_s       set to NaN
%     a_stop_cm_s2   set to NaN
%
%   The raw numbers are preserved for QA in v0_meas_cm_s, d_final_meas_cm,
%   t_stop_meas_s and a_stop_meas_cm_s2. A large v0_meas_cm_s flags a trial
%   that was dropped rather than released. NEVER fit against the *_meas
%   columns; they are diagnostics, not measurements.
%
%   BACKLOG: d0 needs a dedicated quasi-static extraction -- the settled
%   penetration under the foot's own weight, with no impact event. Until that
%   exists, zero-drop rows carry no usable depth.
%
%   dropHeight_mm is the physical drop height and is the column every height
%   axis should use. trialTag is an identifier only.
%
%   Replaces load_default_gb, which hardcoded GB/full + Default geometry and
%   one manual exclusion list. The inclusion rules that were baked into that
%   function are available here as options, so the analysis scripts no longer
%   each carry their own copy:
%
%   OPTIONS
%       'model'       '' (all) | 'Default' | 'Tight' | 'Wide'
%       'condition'   '' (all) | 'GB/full' | 'GB/shallow' | 'CHIN/dense' | ...
%       'exclude'     cellstr/string array of trialTags to drop outright.
%                     DEFAULTS to get_manual_exclusions(), the reviewed list,
%                     so every analysis through this loader drops the same
%                     trials without keeping its own copy. Passing this option
%                     REPLACES that list rather than adding to it -- pass
%                     [get_manual_exclusions(); "extra_tag"] to extend it, or
%                     strings(0,1) to disable manual exclusion entirely.
%       'series'      false; when true also loads kin.z / kin.v / kin.t_s and
%                     returns them as a cell column per trial
%       'dropNaNDepth' true;  drop h > 0 trials with non-finite d_final_cm
%       'dropGlitch'   true;  drop h > 0 trials failing the monotonic-|v| test
%       'depthCut'     [];    drop h > 0 trials deeper than this (cm)
%       'v0RatioMax'   1.15;  free-fall bound for the v0 consistency guard,
%                             as a multiple of sqrt(2*g*h). WARNS, never drops
%       'verbose'      true;  print the per-rule breakdown
%
%   Option names are matched case-insensitively, and an unrecognised name is an
%   ERROR rather than a silently ignored no-op -- a mistyped 'v0RatioMax' would
%   otherwise disable the consistency guard without any sign that it had.
%
%   GEOMETRY. Two tag conventions coexist: the validated drop trials carry no
%   model suffix (165mm_T03_full) while the newer zero-drop trials do
%   (0mm_T01_full_default). Tight and Wide always carry _tight / _wide. So an
%   unsuffixed tag is Default, not "unknown".
%
%   ZERO-DROP TRIALS (h = 0) are returned with isZeroDrop true and are exempt
%   from the automatic rules, which assume a fall and a deceleration that a
%   released-from-rest trial does not have. Manual exclusions still apply. The
%   NaN-depth rule in particular must stay exempt, or the quarantine above
%   would delete the very rows it is meant to preserve.
%
%   The returned table always carries the diagnostic columns keep and reason,
%   so an exclusion can be traced rather than silently applied. Rows failing a
%   rule are REMOVED from K; pass 'verbose' to see the counts.

opt.model        = '';
opt.condition    = '';
opt.exclude      = get_manual_exclusions();
opt.series       = false;
opt.dropNaNDepth = true;
opt.dropGlitch   = true;
opt.depthCut     = [];
opt.v0RatioMax   = 1.15;
opt.verbose      = true;

% Case-insensitive, and unknown names stop rather than being absorbed as new
% struct fields nothing ever reads. The v0 guard below is the reason this
% matters: an option that silently does nothing is worse than no option.
optNames = fieldnames(opt);
for i = 1:2:numel(varargin)
    j = find(strcmpi(optNames, varargin{i}), 1);
    if isempty(j)
        error('load_kinematics_set:unknownOption', ...
              'Unknown option "%s". Valid options: %s', ...
              string(varargin{i}), strjoin(optNames', ', '));
    end
    opt.(optNames{j}) = varargin{i+1};
end

exclude = string(opt.exclude(:));

fprintf('\n=== load_kinematics_set ===\n');

% ── gather candidates ────────────────────────────────────────────────────
F = dir(fullfile(root,'03_RESULTS','**','*_kin_scalars.csv'));
F = F(~[F.isdir]);
if isempty(F)
    error('load_kinematics_set:noFiles', ...
          'No *_kin_scalars.csv under %s', fullfile(root,'03_RESULTS'));
end

n = numel(F);
tag   = strings(n,1);  cond  = strings(n,1);  mdl = strings(n,1);
hLab  = nan(n,1);      d     = nan(n,1);      v0  = nan(n,1);
tstop = nan(n,1);      phi   = nan(n,1);      fps = nan(n,1);
impF  = nan(n,1);      stopF = nan(n,1);      aSt = nan(n,1);
idp   = nan(n,1);      bedX  = nan(n,1);      kp  = strings(n,1);

for i = 1:n
    T = readtable(fullfile(F(i).folder,F(i).name));
    if height(T) < 1, continue; end
    tag(i)   = string(erase(F(i).name,'_kin_scalars.csv'));
    cond(i)  = local_str(T,'condition');
    mdl(i)   = local_model_from_tag(tag(i));
    d(i)     = local_num(T,'d_final_cm');
    v0(i)    = abs(local_num(T,'v0_cm_s'));
    tstop(i) = local_num(T,'t_stop_s');
    phi(i)   = local_num(T,'phi');
    fps(i)   = local_num(T,'fps');
    impF(i)  = local_num(T,'impact_index');
    stopF(i) = local_num(T,'stopFrame');
    aSt(i)   = local_num(T,'a_stop_cm_s2');
    idp(i)   = local_num(T,'impactDistPx');
    bedX(i)  = local_num(T,'bedX');
    hLab(i)  = local_num(T,'dropHeight_mm');
    kp(i)    = string(strrep(fullfile(F(i).folder,F(i).name), ...
                             '_kin_scalars.csv','_kin.mat'));
end

ok    = tag ~= "";
tag=tag(ok); cond=cond(ok); mdl=mdl(ok); hLab=hLab(ok); d=d(ok); v0=v0(ok);
tstop=tstop(ok); phi=phi(ok); fps=fps(ok); impF=impF(ok); stopF=stopF(ok);
aSt=aSt(ok); idp=idp(ok); bedX=bedX(ok); kp=kp(ok);

K = table(tag, mdl, cond, hLab, fps, impF, stopF, ...
          v0, d, tstop, aSt, phi, idp, bedX, kp, ...
    'VariableNames', {'trialTag','model','condition','dropHeight_mm', ...
        'fps','impactFrame','stopFrame', ...
        'v0_cm_s','d_final_cm','t_stop_s','a_stop_cm_s2','phi', ...
        'impactDistPx','bedX','kinPath'});

% ── duplicate trialTags are a hard error ─────────────────────────────────
% One trial must contribute exactly one row. A tag appearing twice means two
% _kin_scalars.csv files describe the same trial -- typically a stale copy left
% under an old batch or model folder after a re-run or a rename. Silently
% keeping both double-weights that trial in every per-height mean and every
% fit, so this stops rather than warns: the caller cannot tell from the results
% that it happened, and the fix is a deliberate deletion, not a filter.
[uTags, ~, uIdx] = unique(K.trialTag);
counts = accumarray(uIdx, 1);
dupTags = uTags(counts > 1);
if ~isempty(dupTags)
    msg = sprintf(['%d trialTag(s) appear more than once. One trial must ' ...
                   'produce exactly one row.\n\n'], numel(dupTags));
    for i = 1:numel(dupTags)
        rows = find(K.trialTag == dupTags(i));
        msg  = [msg sprintf('  %s  (%d copies)\n', dupTags(i), numel(rows))]; %#ok<AGROW>
        for r = rows(:).'
            % kinPath was derived from the scalars CSV; invert that to name the
            % file the user actually has to delete.
            csvPath = strrep(K.kinPath(r), '_kin.mat', '_kin_scalars.csv');
            msg = [msg sprintf('      %s\n', csvPath)]; %#ok<AGROW>
        end
    end
    msg = [msg sprintf(['\nDelete the stale copies (keep exactly one per trial) ' ...
                        'and re-run.\nDo not filter them out here: which copy is ' ...
                        'current is a decision only you can make.\n'])];
    error('load_kinematics_set:duplicateTrialTag', '%s', msg);
end

% ── read-time height relabeling ──────────────────────────────────────────
% A few trials carry a recorded height that review established is wrong.
% get_relabel_map holds the corrections, one line per trial with its reason
% and date; this is the ONLY place they are applied, and nothing on disk is
% touched (the trialTag itself still says the old height -- deliberately, so
% provenance back to the capture session survives).
%
% Applied HERE, before isZeroDrop and before the quarantine, because the
% corrected height is the physical truth: if a relabel ever moved a trial to
% or from h = 0, every rule downstream must see the corrected value, not the
% recorded one.
R = get_relabel_map();
K.dropHeight_asRecorded = K.dropHeight_mm;
K.relabeled             = false(height(K),1);
if height(R) > 0
    [tf, loc] = ismember(K.trialTag, R.trialTag);
    K.dropHeight_mm(tf) = R.dropHeight_mm(loc(tf));
    K.relabeled         = tf;   % true = this height came from the map
    if any(tf)
        fprintf('  relabeled: %d trial(s) (see get_relabel_map)\n', sum(tf));
        if opt.verbose
            for i = find(tf(:).')
                fprintf('             %-28s %g -> %g mm\n', K.trialTag(i), ...
                    K.dropHeight_asRecorded(i), K.dropHeight_mm(i));
            end
        end
    end
    % A mapped tag that matches no file is a correction doing nothing at all.
    % Report it the way the manual-exclusion list reports its own misses.
    missR = R.trialTag(~ismember(R.trialTag, K.trialTag));
    if ~isempty(missR)
        fprintf('  relabel entries matching no file : %d\n', numel(missR));
        fprintf('             %s\n', missR);
    end
end

K.isZeroDrop = K.dropHeight_mm == 0;

% ── zero-drop quarantine ─────────────────────────────────────────────────
% A d0 trial is released FROM CONTACT: there is no fall, no impact, and no
% deceleration transient. kd_kinematics is a DYNAMIC-IMPACT pipeline and is
% correct for its domain -- h = 0 is simply outside it. Everything it derives
% from the impact/stop event is therefore undefined here, not merely noisy:
%
%   v0_cm_s       the velocity peak of a release transient, not an impact speed
%   d_final_cm    depth at a "stop" located by a zero-crossing that never
%                 belonged to a deceleration
%   t_stop_s      extrapolated from that same non-event
%   a_stop_cm_s2  the slope of a fit to it
%
% All four are set to their protocol value (v0 = 0) or to NaN, and the raw
% numbers are preserved in *_meas companion columns for QA. A large
% v0_meas_cm_s, for instance, flags a trial that was dropped rather than
% released.
%
% BACKLOG: d0 needs a dedicated quasi-static extraction -- the settled
% penetration under the foot's own weight, measured without an impact event.
% Until that exists these rows carry no usable depth, and nothing downstream
% should treat the quarantined values as measurements.
K.v0_meas_cm_s      = K.v0_cm_s;
K.d_final_meas_cm   = K.d_final_cm;
K.t_stop_meas_s     = K.t_stop_s;
K.a_stop_meas_cm_s2 = K.a_stop_cm_s2;

K.v0_cm_s(K.isZeroDrop)       = 0;      % protocol: no fall, so no impact speed
K.d_final_cm(K.isZeroDrop)    = NaN;    % undefined until quasi-static d0 exists
K.t_stop_s(K.isZeroDrop)      = NaN;
K.a_stop_cm_s2(K.isZeroDrop)  = NaN;

% ── phi: take the CURRENT value, not the one baked into the CSV ──────────
% _kin_scalars.csv carries meta.phi as it stood at TRACKING time, which for
% this dataset is the retired single-measurement figure (0.629 / 0.636 /
% 0.276 / 0.409). get_substrate_properties holds the finalised five-prep
% means (0.624 / 0.643 / 0.280 / 0.402), so re-read it here and overwrite.
K.phi_csv = K.phi;                       % kept for provenance
for i = 1:height(K)
    parts = split(K.condition(i), "/");
    if numel(parts) ~= 2, continue; end
    try
        sub = get_substrate_properties(char(parts(1)), char(parts(2)));
        if isstruct(sub) && isfield(sub,'phi') && isfinite(sub.phi)
            K.phi(i) = sub.phi;
        end
    catch
        % leave the CSV value in place if the lookup fails
    end
end
nPhi = sum(isfinite(K.phi) & isfinite(K.phi_csv) & abs(K.phi-K.phi_csv) > 1e-6);
if nPhi > 0
    fprintf('  phi refreshed from get_substrate_properties : %d row(s)\n', nPhi);
end

fprintf('kin_scalars found         : %d\n', height(K));

% ── model / condition filter ─────────────────────────────────────────────
if ~isempty(opt.model)
    K = K(strcmpi(K.model, string(opt.model)), :);
    fprintf('  model == %-10s      : %d\n', string(opt.model), height(K));
end
if ~isempty(opt.condition)
    K = K(K.condition == string(opt.condition), :);
    fprintf('  condition == %-10s  : %d\n', string(opt.condition), height(K));
end
if isempty(K)
    error('load_kinematics_set:noneSelected', ...
          'No trials left after the model/condition filter.');
end

% ── rules ────────────────────────────────────────────────────────────────
exManual = ismember(K.trialTag, exclude);

exNaN = false(height(K),1);
if opt.dropNaNDepth
    exNaN = ~K.isZeroDrop & ~isfinite(K.d_final_cm);
end

exGlitch = false(height(K),1);
if opt.dropGlitch
    for i = 1:height(K)
        if K.isZeroDrop(i) || exManual(i) || exNaN(i), continue; end
        exGlitch(i) = local_isglitch(K.kinPath(i));
    end
end

exCut = false(height(K),1);
if ~isempty(opt.depthCut)
    exCut = ~K.isZeroDrop & K.d_final_cm > opt.depthCut;
end

K.keep   = ~(exManual | exNaN | exGlitch | exCut);
K.reason = strings(height(K),1);
K.reason(exManual) = "manual";
K.reason(exNaN)    = "NaNDEPTH";
K.reason(exGlitch) = "GLITCH";
K.reason(exCut)    = "depth cut";

if opt.verbose
    fprintf('\n--- exclusions applied ---\n');
    fprintf('  manual   : %2d of %d listed  (%d listed tags matched no file)\n', ...
        sum(exManual), numel(exclude), sum(~ismember(exclude, K.trialTag)));
    miss = exclude(~ismember(exclude, K.trialTag));
    if ~isempty(miss), fprintf('             %s\n', miss); end
    if opt.dropNaNDepth, fprintf('  NaNDEPTH : %2d  (h > 0 only)\n', sum(exNaN));    end
    if opt.dropGlitch,   fprintf('  GLITCH   : %2d  (h > 0 only)\n', sum(exGlitch)); end
    if ~isempty(opt.depthCut)
        fprintf('  depth cut: %2d  (> %.2f cm)\n', sum(exCut), opt.depthCut);
    end
end

K = K(K.keep, :);

% ── v0 consistency guard ─────────────────────────────────────────────────
% Every kept h > 0 trial is checked against the free-fall speed its label
% implies. This is a GUARD, not a rule: it warns and counts, it never drops.
% Which trials to remove is a review decision that belongs in
% get_manual_exclusions with a reason attached, not in an automatic filter
% that would quietly reshape a bin the next time a threshold moved.
%
% Two failure modes, and they are not symmetric:
%
%   HIGH SIDE  v0 > v0RatioMax * sqrt(2*g*h) is a physical impossibility once
%              the tolerance is allowed for -- a body released from h cannot
%              arrive faster than free fall. A violation means the LABEL is
%              wrong (contamination between height bins) or the TIME BASE is
%              wrong (fps, or dropped frames as in campaign 1). Either way the
%              trial's v0 and everything derived from it are suspect, and this
%              must never pass silently: it is exactly the failure that the
%              relabeling work of 2026-08-21 was chasing.
%
%   LOW SIDE   v0 < 0.5 * sqrt(2*g*h) is not impossible -- release friction, a
%              snagged line, a late trigger. It is a detector for a spoiled
%              release (the known 85 mm Default trial at 39.9 cm/s), so the
%              threshold is deliberately loose: half of free fall is far below
%              any honest air-drag loss over these heights.
%
% Real drops land a few percent BELOW free fall, so the healthy ratio is just
% under 1. Two campaigns sit systematically above it because their release
% height reference was offset -- see the data notes in README.md; those
% warnings are documented behaviour, not bad trials.
G_CM_S2 = 980;    % as calib.g_cm_s2 (get_calibration.m) and net_accel's fallback

vFree = sqrt(2 * G_CM_S2 * K.dropHeight_mm / 10);   % cm/s, from the CORRECTED height
chkV0 = ~K.isZeroDrop & K.dropHeight_mm > 0 & isfinite(vFree) & isfinite(K.v0_cm_s);

K.v0_freefall_ratio = nan(height(K),1);             % v0 / sqrt(2gh); NaN where undefined
K.v0_freefall_ratio(chkV0) = K.v0_cm_s(chkV0) ./ vFree(chkV0);

hiV0 = chkV0 & K.v0_cm_s > opt.v0RatioMax * vFree;
loV0 = chkV0 & K.v0_cm_s < 0.5 * vFree;

for i = find(hiV0(:).')
    warning('load_kinematics_set:v0AboveFreeFall', ...
        ['%s: v0 = %.1f cm/s EXCEEDS the free-fall bound %.1f cm/s ' ...
         '(%.2f x sqrt(2gh) at h = %g mm, limit %.2f). A drop cannot be ' ...
         'faster than free fall: suspect the height label or the time base.'], ...
        K.trialTag(i), K.v0_cm_s(i), opt.v0RatioMax*vFree(i), ...
        K.v0_freefall_ratio(i), K.dropHeight_mm(i), opt.v0RatioMax);
end
for i = find(loV0(:).')
    warning('load_kinematics_set:v0BelowHalfFreeFall', ...
        ['%s: v0 = %.1f cm/s is below half the free-fall speed %.1f cm/s ' ...
         '(%.2f x sqrt(2gh) at h = %g mm). Suspect a snagged or spoiled ' ...
         'release.'], ...
        K.trialTag(i), K.v0_cm_s(i), vFree(i), ...
        K.v0_freefall_ratio(i), K.dropHeight_mm(i));
end

if opt.verbose
    fprintf('\n--- v0 consistency guard (warns, never drops) ---\n');
    fprintf('  checked            : %3d kept h > 0 trials\n', sum(chkV0));
    fprintf('  above free fall    : %3d  (v0 > %.2f x sqrt(2gh))\n', ...
        sum(hiV0), opt.v0RatioMax);
    fprintf('  below half free fall: %2d  (v0 < 0.50 x sqrt(2gh))\n', sum(loV0));
end

% ── optional time series ─────────────────────────────────────────────────
if opt.series
    z = cell(height(K),1); v = cell(height(K),1); t = cell(height(K),1);
    for i = 1:height(K)
        [z{i}, v{i}, t{i}] = local_series(K.kinPath(i));
    end
    K.z = z; K.v = v; K.t_s = t;
end

% ── report ───────────────────────────────────────────────────────────────
if opt.verbose
    fprintf('\n--- loaded set ---\n');
    DR = K(~K.isZeroDrop,:);  D0 = K(K.isZeroDrop,:);
    fprintf('  drop-height trials : %3d over %d heights\n', ...
        height(DR), numel(unique(DR.dropHeight_mm)));
    fprintf('  zero-drop trials   : %3d\n', height(D0));
    cs = unique(K.condition);
    for c = 1:numel(cs)
        fprintf('    %-16s %3d\n', cs(c), sum(K.condition==cs(c)));
    end
    if ~isempty(D0)
        % Reads the preserved column: d_final_cm is NaN on these rows by
        % quarantine, so the raw numbers live in d_final_meas_cm.
        fprintf('\n  QUARANTINED d0 (h = 0, released from contact):\n');
        fprintf('    %.3f +/- %.3f cm   (n = %d, range %.3f - %.3f)\n', ...
            mean(D0.d_final_meas_cm), std(D0.d_final_meas_cm), height(D0), ...
            min(D0.d_final_meas_cm), max(D0.d_final_meas_cm));
        fprintf(['    NOT A MEASUREMENT: zero-drop kinematics come from the\n' ...
                 '    dynamic-impact pipeline, which has no impact to work from\n' ...
                 '    at h = 0. Quarantined pending a dedicated quasi-static d0\n' ...
                 '    extraction (backlog). d_final_cm, t_stop_s and a_stop_cm_s2\n' ...
                 '    are NaN on these rows; the raw values are in the *_meas\n' ...
                 '    columns for QA only.\n']);
    end
    fprintf('\n');
end
end

% ─────────────────────────────────────────────────────────────────────────
function m = local_model_from_tag(tag)
% Tight and Wide always carry an explicit suffix; everything else is Default.
if endsWith(tag,"_tight"),   m = "Tight";
elseif endsWith(tag,"_wide"), m = "Wide";
else,                         m = "Default";
end
end

function tf = local_isglitch(kinPath)
%LOCAL_ISGLITCH  |v| should fall monotonically from v0 to zero between impact
%   and stop. An increase above 15% of v0 marks a tracking artefact.
tf = false;
if ~isfile(kinPath), return; end
try
    S = load(kinPath,'kin');
catch
    return
end
k = S.kin;
if ~all(isfield(k,{'v','impact_index','stopFrame','v0_cm_s'})), return; end
i1 = max(1, k.impact_index);
i2 = min(numel(k.v), k.stopFrame);
if i2 <= i1 + 1, return; end
seg = abs(k.v(i1:i2));
seg = seg(isfinite(seg));
if numel(seg) < 3, return; end
tf = any(diff(seg) > 0.15*abs(k.v0_cm_s));
end

function [z, v, t] = local_series(kinPath)
z = []; v = []; t = [];
if ~isfile(kinPath), return; end
try
    S = load(kinPath,'kin');
catch
    return
end
k = S.kin;
if isfield(k,'z_cm'),    z = k.z_cm;    elseif isfield(k,'z'), z = k.z; end
if isfield(k,'v_cm_s'),  v = k.v_cm_s;  elseif isfield(k,'v'), v = k.v; end
if isfield(k,'t_s'),     t = k.t_s;     end
end

function v = local_num(T,name)
v = NaN;
if ismember(name,T.Properties.VariableNames)
    x = T.(name)(1);
    if isnumeric(x)||islogical(x), v = double(x); else, v = str2double(string(x)); end
end
end

function s = local_str(T,name)
s = "";
if ismember(name,T.Properties.VariableNames)
    s = strtrim(string(T.(name)(1)));
end
end
