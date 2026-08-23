function M = export_master_dataset(root, varargin)
%EXPORT_MASTER_DATASET  One row per trial, every quantity the scaling study needs.
%
%   READ-ONLY with respect to the results tree: it reads *_kin.mat,
%   *_kin_scalars.csv and the foot-area CSVs, and writes only to OutDir.
%
%   WHY THIS EXISTS. The scaling analysis is to be built from the trial-level
%   data, not from per-height curves or frozen figures. This script gathers
%   every trial the pipeline has produced -- kept, excluded and quarantined --
%   into one table with the exclusion decision, its reason and the provenance
%   carried on each row, so that any later selection is explicit and
%   reversible. Nothing is filtered out here; the keep_reviewed column IS the
%   reviewed selection, and it is derived from load_kinematics_set with its
%   default rules so the two cannot drift apart.
%
%   USAGE
%       M = export_master_dataset('D:\ME_GRANULAB\JerboaImpact');
%       M = export_master_dataset(root, 'Series', false);
%
%   OPTIONS (name-value)
%       'OutDir'    output folder             (default <root>/03_RESULTS/_exports)
%       'Series'    also write the per-frame series (long CSV + MAT) (default true)
%       'RawRoot'   capture tree; when given, each trial's raw .avi is located
%                   with find_raw_video and its file timestamp is recorded as
%                   the capture time (default getenv('JERBOA_RAW_ROOT'); '' skips)
%       'MassG'     projectile mass, g        (default 65, as track_tracers_2)
%       'Stamp'     timestamp the file names  (default true)
%
%   OUTPUTS (in OutDir; <stamp> = yyyymmdd_HHMMSS)
%       master_trials_<stamp>.csv       one row per trial, all scalar columns
%       master_trials_<stamp>.mat       T (the table), S (per-trial series),
%                                       P (export provenance)
%       master_series_<stamp>.csv       long format: one row per (trial, frame)
%       master_dictionary_<stamp>.csv   column, units, source, definition
%
%   COLUMN GROUPS -- see the dictionary for every column.
%       identity     trialTag, model, condition, material, container, batch,
%                    trialNum, heightLabel
%       height       dropHeight_mm (corrected), dropHeight_asRecorded, relabeled,
%                    isZeroDrop, v0_ff_cm_s, v0_freefall_ratio
%       kinematics   v0_cm_s, d_final_cm, t_stop_s, a_stop_cm_s2, d_max_cm,
%                    rebound_cm, apg_peak_cm_s2, z_at_apg_peak_cm, n_frames_impact_to_stop
%       quarantine   *_meas companions (raw pipeline values on h = 0 rows)
%       selection    keep_reviewed, exclude_reason, excl_manual, excl_nanDepth,
%                    excl_glitch, v0_guard
%       substrate    phi, phi_csv, rho_particle_g_cm3, rho_bulk_g_cm3
%       geometry     A_bare_cm2, A_hull_cm2, hull_over_bare, A_bare_at_dfinal_cm2,
%                    A_hull_at_dfinal_cm2, a_local_at_dfinal_cm2, mass_g
%       provenance   fps, dt_s, impactFrame, stopFrame, windowStart/End,
%                    eventFrameFirst/Last, mmPerPx, impactDistPx, bedX, g_cm_s2,
%                    refMarkerID, rodAngle_*, gitCommit, processedOn, detectPass,
%                    kinMethod, kinPath, capture_datetime
%
%   Base MATLAB only.

% ── options ──────────────────────────────────────────────────────────────
opt.OutDir  = fullfile(root, '03_RESULTS', '_exports');
opt.Series  = true;
opt.RawRoot = getenv('JERBOA_RAW_ROOT');
opt.MassG   = 65;
opt.Stamp   = true;
optNames = fieldnames(opt);
for i = 1:2:numel(varargin)
    j = find(strcmpi(optNames, varargin{i}), 1);
    if isempty(j)
        error('export_master_dataset:unknownOption', ...
              'Unknown option "%s". Valid: %s', string(varargin{i}), strjoin(optNames', ', '));
    end
    opt.(optNames{j}) = varargin{i+1};
end
G_CM_S2     = 980;      % as load_kinematics_set and get_calibration
V0_RATIO_HI = 1.15;     % load_kinematics_set default v0RatioMax
V0_RATIO_LO = 0.50;     % load_kinematics_set low-side threshold

stamp = datestr(now, 'yyyymmdd_HHMMSS');
sfx   = '';
if opt.Stamp, sfx = ['_' stamp]; end

fprintf('\n=== export_master_dataset ===\n');
fprintf('  root    : %s\n', root);
fprintf('  OutDir  : %s\n', opt.OutDir);
fprintf('  RawRoot : %s\n', local_tern(isempty(opt.RawRoot), '(none; capture times skipped)', opt.RawRoot));

% ── 1. every trial, no rule applied ──────────────────────────────────────
% Relabeling, zero-drop quarantine and the phi refresh still happen inside
% the loader; only the row-dropping rules are switched off.
K = load_kinematics_set(root, 'exclude', strings(0,1), ...
                        'dropNaNDepth', false, 'dropGlitch', false, 'verbose', false);

% ── 2. the reviewed selection, from the loader's defaults ────────────────
% Second call with the default rules. Its warnings duplicate the first call's,
% so they are silenced here only.
ws = warning; restore = onCleanup(@() warning(ws)); %#ok<NASGU>
warning('off', 'load_kinematics_set:v0AboveFreeFall');
warning('off', 'load_kinematics_set:v0BelowHalfFreeFall');
Kdef = load_kinematics_set(root, 'verbose', false);
clear restore

n = height(K);
manualList = get_manual_exclusions();
excl_manual   = ismember(K.trialTag, manualList);
excl_nanDepth = ~K.isZeroDrop & ~isfinite(K.d_final_cm);
keep_reviewed = ismember(K.trialTag, Kdef.trialTag);
% The loader applies exactly three rules by default (manual, NaN depth,
% glitch). A row absent from the reviewed set for neither of the first two
% reasons failed the glitch test.
excl_glitch   = ~keep_reviewed & ~excl_manual & ~excl_nanDepth;

exclude_reason = strings(n,1);
exclude_reason(excl_glitch)   = "GLITCH";
exclude_reason(excl_nanDepth) = "NaNDEPTH";
exclude_reason(excl_manual)   = "manual";

v0_guard = strings(n,1);
v0_guard(isfinite(K.v0_freefall_ratio) & K.v0_freefall_ratio > V0_RATIO_HI) = "above_freefall";
v0_guard(isfinite(K.v0_freefall_ratio) & K.v0_freefall_ratio < V0_RATIO_LO) = "below_half_freefall";

v0_ff = sqrt(2 * G_CM_S2 * K.dropHeight_mm / 10);

fprintf('  trials found      : %d\n', n);
fprintf('  reviewed set      : %d kept  (%d manual, %d NaN depth, %d glitch)\n', ...
        sum(keep_reviewed), sum(excl_manual), sum(excl_nanDepth), sum(excl_glitch));

% ── 3. per-trial reads of _kin.mat and _kin_scalars.csv ──────────────────
P = local_empty_provenance(n);
S = struct('trialTag', cell(n,1), 'frame', [], 't_s', [], 'z_cm', [], ...
           'v_cm_s', [], 'a_cm_s2', [], 'a_plus_g_cm_s2', [], ...
           'accel_se_cm_s2', [], 'winHalfFrames', [], 'rodBedDist_px', []);
fprintf('  reading %d _kin.mat ...\n', n);
for i = 1:n
    kp = char(K.kinPath(i));
    if ~isfile(kp)
        warning('export_master_dataset:missingKin', '%s: no _kin.mat at %s', K.trialTag(i), kp);
        continue
    end
    L = load(kp);          % meta, kin, calib
    kin = L.kin; meta = L.meta; calib = L.calib;

    % series (all frames where depth is defined; already masked by kd_kinematics)
    z = kin.z(:); v = kin.v(:); a = kin.a(:); apg = kin.a_plus_g(:); t = kin.t_s(:);
    fr = (1:numel(z)).';
    ok = isfinite(z);
    S(i).trialTag        = K.trialTag(i);
    S(i).frame           = fr(ok);
    S(i).t_s             = t(ok);
    S(i).z_cm            = z(ok);
    S(i).v_cm_s          = v(ok);
    S(i).a_cm_s2         = a(ok);
    S(i).a_plus_g_cm_s2  = apg(ok);
    S(i).accel_se_cm_s2  = local_field(kin, 'accel_se', nan(size(z)));  S(i).accel_se_cm_s2 = S(i).accel_se_cm_s2(ok);
    S(i).winHalfFrames   = local_field(kin, 'winHalfFrames', nan(size(z))); S(i).winHalfFrames = S(i).winHalfFrames(ok);
    S(i).rodBedDist_px   = local_field(kin, 'rodBedDist_px', nan(size(z))); S(i).rodBedDist_px = S(i).rodBedDist_px(ok);

    % series-derived scalars, over [impact, end of retained record]
    i1 = max(1, kin.impact_index); i2 = min(numel(z), kin.stopFrame);
    seg = i1:numel(z);
    [dmax, kmax] = max(z(seg), [], 'omitnan');
    P.d_max_cm(i)      = dmax;
    P.frame_at_dmax(i) = seg(kmax);
    P.rebound_cm(i)    = dmax - kin.d_final_cm;
    segIS = i1:i2;
    [apk, kapk] = max(apg(segIS), [], 'omitnan');
    P.apg_peak_cm_s2(i)    = apk;
    P.z_at_apg_peak_cm(i)  = z(segIS(kapk));
    P.n_frames_impact_to_stop(i) = i2 - i1 + 1;
    P.dt_s(i)          = local_field(kin, 'dt', NaN);
    P.refMarkerID(i)   = local_field(kin, 'refMarkerID', NaN);
    P.kinMethod(i)     = string(local_field(kin, 'method', ''));

    % calibration
    P.mmPerPx(i)  = local_field(calib, 'mmPerPx', NaN);
    P.g_cm_s2(i)  = local_field(calib, 'g_cm_s2', NaN);

    % meta / provenance
    P.material(i)    = string(local_field(meta, 'material', ''));
    P.container(i)   = string(local_field(meta, 'container', ''));
    P.batch(i)       = string(local_field(meta, 'batchName', ''));
    P.trialNum(i)    = local_field(meta, 'trialNum', NaN);
    P.heightLabel(i) = string(local_field(meta, 'heightLabel', ''));
    P.nFrames(i)     = local_field(meta, 'nFrames', NaN);
    P.nTracked(i)    = local_field(meta, 'nTracked', NaN);
    P.firstValidFrame(i) = local_field(meta, 'firstValidFrame', NaN);
    P.windowStart(i) = local_field(meta, 'windowStart', NaN);
    P.windowEnd(i)   = local_field(meta, 'windowEnd', NaN);
    P.autoWindow(i)  = isequal(local_field(meta, 'autoWindow', false), true);
    ef = local_field(meta, 'eventFrameRange', [NaN NaN]);
    if numel(ef) == 2, P.eventFrameFirst(i) = ef(1); P.eventFrameLast(i) = ef(2); end
    prov = local_field(meta, 'provenance', struct());
    P.gitCommit(i)   = string(local_field(prov, 'gitCommit', ''));
    P.processedOn(i) = string(local_field(prov, 'processedOn', ''));
    P.detectPass(i)  = string(local_field(prov, 'detectPass', ''));
    P.frameSource(i) = string(local_field(prov, 'frameSource', ''));

    % rod-angle scalars: only in the scalars CSV
    csvp = strrep(kp, '_kin.mat', '_kin_scalars.csv');
    if isfile(csvp)
        T1 = readtable(csvp);
        P.rodAngle_peak_abs_deg(i)  = local_num(T1, 'rodAngle_peak_abs_deg');
        P.rodAngle_net_deg(i)       = local_num(T1, 'rodAngle_net_deg');
        P.rodAngle_range_deg(i)     = local_num(T1, 'rodAngle_range_deg');
        P.rodAngle_noise_sd_deg(i)  = local_num(T1, 'rodAngle_noise_sd_deg');
    end

    % capture time from the raw clip's file timestamp (session provenance)
    if ~isempty(opt.RawRoot)
        try
            vp = find_raw_video(opt.RawRoot, meta, char(K.model(i)));
            d  = dir(vp);
            P.capture_datetime(i) = string(datestr(d.datenum, 'yyyy-mm-dd HH:MM:SS'));
            P.rawVideoPath(i)     = string(vp);
        catch
            % leave empty: provenance only, never a reason to stop an export
        end
    end
end

% ── 4. substrate properties (densities; phi already refreshed by the loader)
rho_p = nan(n,1); rho_b = nan(n,1);
for i = 1:n
    parts = split(K.condition(i), "/");
    if numel(parts) ~= 2, continue; end
    sub = get_substrate_properties(char(parts(1)), char(parts(2)));
    if sub.ok, rho_p(i) = sub.rho_particle_g_cm3; rho_b(i) = sub.rho_bulk_g_cm3; end
end

% ── 5. foot areas: single source = the schematic script's own outputs ────
figDir   = fullfile(root, '03_RESULTS', '_figures');
areaCsv  = fullfile(figDir, 'foot_areas_computed.csv');
sweepCsv = fullfile(figDir, 'foot_area_vs_depth.csv');
if ~isfile(areaCsv) || ~isfile(sweepCsv)
    error('export_master_dataset:noAreas', ...
          ['Foot-area CSVs not found under %s.\nRun fig_foot_schematic(''Root'', root) ' ...
           'first: the export reads the computed areas rather than typing them in.'], figDir);
end
A  = readtable(areaCsv,  'CommentStyle', '#');
SW = readtable(sweepCsv, 'CommentStyle', '#');
A_bare = nan(n,1); A_hull = nan(n,1);
A_bare_d = nan(n,1); A_hull_d = nan(n,1); a_loc_d = nan(n,1);
for i = 1:n
    r = find(strcmpi(string(A.model), K.model(i)), 1);
    if isempty(r), continue; end
    A_bare(i) = A.A_bare_cm2(r); A_hull(i) = A.A_hull_cm2(r);
    m = char(K.model(i));
    if isfinite(K.d_final_cm(i))
        zc = SW.z_mm / 10;
        A_bare_d(i) = interp1(zc, SW.(['A_bare_'  m]), K.d_final_cm(i), 'linear', NaN);
        A_hull_d(i) = interp1(zc, SW.(['A_hull_'  m]), K.d_final_cm(i), 'linear', NaN);
        a_loc_d(i)  = interp1(zc, SW.(['a_local_' m]), K.d_final_cm(i), 'linear', NaN);
    end
end
cutY = A.CutY_mm(1);

% ── 6. assemble the table ────────────────────────────────────────────────
T = table();
T.trialTag              = K.trialTag;
T.model                 = K.model;
T.condition             = K.condition;
T.material              = P.material;
T.container             = P.container;
T.batch                 = P.batch;
T.trialNum              = P.trialNum;
T.heightLabel           = P.heightLabel;
T.dropHeight_mm         = K.dropHeight_mm;
T.dropHeight_asRecorded = K.dropHeight_asRecorded;
T.relabeled             = K.relabeled;
T.isZeroDrop            = K.isZeroDrop;
T.v0_ff_cm_s            = v0_ff;
T.v0_freefall_ratio     = K.v0_freefall_ratio;
T.v0_cm_s               = K.v0_cm_s;
T.d_final_cm            = K.d_final_cm;
T.t_stop_s              = K.t_stop_s;
T.a_stop_cm_s2          = K.a_stop_cm_s2;
T.d_max_cm              = P.d_max_cm;
T.rebound_cm            = P.rebound_cm;
T.apg_peak_cm_s2        = P.apg_peak_cm_s2;
T.z_at_apg_peak_cm      = P.z_at_apg_peak_cm;
T.n_frames_impact_to_stop = P.n_frames_impact_to_stop;
T.v0_meas_cm_s          = K.v0_meas_cm_s;
T.d_final_meas_cm       = K.d_final_meas_cm;
T.t_stop_meas_s         = K.t_stop_meas_s;
T.a_stop_meas_cm_s2     = K.a_stop_meas_cm_s2;
T.keep_reviewed         = keep_reviewed;
T.exclude_reason        = exclude_reason;
T.excl_manual           = excl_manual;
T.excl_nanDepth         = excl_nanDepth;
T.excl_glitch           = excl_glitch;
T.v0_guard              = v0_guard;
T.phi                   = K.phi;
T.phi_csv               = K.phi_csv;
T.rho_particle_g_cm3    = rho_p;
T.rho_bulk_g_cm3        = rho_b;
T.A_bare_cm2            = A_bare;
T.A_hull_cm2            = A_hull;
T.hull_over_bare        = A_hull ./ A_bare;
T.A_bare_at_dfinal_cm2  = A_bare_d;
T.A_hull_at_dfinal_cm2  = A_hull_d;
T.a_local_at_dfinal_cm2 = a_loc_d;
T.area_cutY_mm          = repmat(cutY, n, 1);
T.mass_g                = repmat(opt.MassG, n, 1);
T.fps                   = K.fps;
T.dt_s                  = P.dt_s;
T.impactFrame           = K.impactFrame;
T.stopFrame             = K.stopFrame;
T.frame_at_dmax         = P.frame_at_dmax;
T.nFrames               = P.nFrames;
T.nTracked              = P.nTracked;
T.firstValidFrame       = P.firstValidFrame;
T.windowStart           = P.windowStart;
T.windowEnd             = P.windowEnd;
T.autoWindow            = P.autoWindow;
T.eventFrameFirst       = P.eventFrameFirst;
T.eventFrameLast        = P.eventFrameLast;
T.mmPerPx               = P.mmPerPx;
T.impactDistPx          = K.impactDistPx;
T.bedX                  = K.bedX;
T.g_cm_s2               = P.g_cm_s2;
T.refMarkerID           = P.refMarkerID;
T.rodAngle_peak_abs_deg = P.rodAngle_peak_abs_deg;
T.rodAngle_net_deg      = P.rodAngle_net_deg;
T.rodAngle_range_deg    = P.rodAngle_range_deg;
T.rodAngle_noise_sd_deg = P.rodAngle_noise_sd_deg;
T.gitCommit             = P.gitCommit;
T.processedOn           = P.processedOn;
T.detectPass            = P.detectPass;
T.frameSource           = P.frameSource;
T.kinMethod             = P.kinMethod;
T.capture_datetime      = P.capture_datetime;
T.rawVideoPath          = P.rawVideoPath;
T.kinPath               = K.kinPath;

T = sortrows(T, {'model', 'condition', 'dropHeight_mm', 'trialNum'});

% ── 7. provenance of the export itself ───────────────────────────────────
Pexp = struct();
Pexp.root          = root;
Pexp.exportedOn    = datestr(now, 'yyyy-mm-dd HH:MM:SS');
Pexp.codeCommit    = pipeline_commit(fileparts(fileparts(mfilename('fullpath'))));
Pexp.options       = opt;
Pexp.manualExclusions = manualList;
Pexp.relabelMap    = get_relabel_map();
Pexp.areaSource    = struct('areas', areaCsv, 'sweep', sweepCsv, 'cutY_mm', cutY);
Pexp.v0GuardThresholds = [V0_RATIO_LO V0_RATIO_HI];
Pexp.g_cm_s2       = G_CM_S2;

% ── 8. write ─────────────────────────────────────────────────────────────
if ~isfolder(opt.OutDir), mkdir(opt.OutDir); end
written = strings(0,1);

csvPath = fullfile(opt.OutDir, ['master_trials' sfx '.csv']);
writetable(T, csvPath);
written(end+1) = string(csvPath);

dictPath = fullfile(opt.OutDir, ['master_dictionary' sfx '.csv']);
writetable(local_dictionary(), dictPath);
written(end+1) = string(dictPath);

matPath = fullfile(opt.OutDir, ['master_trials' sfx '.mat']);
P = Pexp; %#ok<NASGU>
if opt.Series
    save(matPath, 'T', 'S', 'P', '-v7.3');
else
    save(matPath, 'T', 'P', '-v7.3');
end
written(end+1) = string(matPath);

if opt.Series
    serPath = fullfile(opt.OutDir, ['master_series' sfx '.csv']);
    local_write_series(serPath, S);
    written(end+1) = string(serPath);
end

% ── 9. summary ───────────────────────────────────────────────────────────
fprintf('\n--- master dataset ---\n');
fprintf('  rows               : %d  (%d reviewed-keep, %d excluded, %d zero-drop)\n', ...
        height(T), sum(T.keep_reviewed & ~T.isZeroDrop), sum(~T.keep_reviewed), sum(T.isZeroDrop));
mdls = unique(T.model);
for k = 1:numel(mdls)
    r = T.model == mdls(k) & T.keep_reviewed & ~T.isZeroDrop;
    if ~any(r), continue; end
    fprintf('  %-8s kept %3d  heights %2d  v0 %5.1f-%5.1f cm/s  d_final %.2f-%.2f cm  A_bare %.3f  A_hull %.3f\n', ...
            mdls(k), sum(r), numel(unique(T.dropHeight_mm(r))), ...
            min(T.v0_cm_s(r)), max(T.v0_cm_s(r)), min(T.d_final_cm(r)), max(T.d_final_cm(r)), ...
            T.A_bare_cm2(find(r,1)), T.A_hull_cm2(find(r,1)));
end
conds = unique(T.condition);
fprintf('  conditions         : %s\n', strjoin(cellstr(conds), ', '));
fprintf('  d_final > foot cut : %d of %d kept trials deeper than z = %.2f cm (bar + post in the bed)\n', ...
        sum(T.keep_reviewed & ~T.isZeroDrop & T.d_final_cm > (cutY + 61.192)/10), ...
        sum(T.keep_reviewed & ~T.isZeroDrop), (cutY + 61.192)/10);
fprintf('  v0 guard           : %d above free fall, %d below half free fall (flags only)\n', ...
        sum(T.v0_guard == "above_freefall"), sum(T.v0_guard == "below_half_freefall"));
fprintf('  capture times      : %d of %d resolved\n', sum(T.capture_datetime ~= ""), height(T));
fprintf('\n--- written ---\n'); fprintf('  %s\n', written);
fprintf('\n');

M = struct('T', T, 'S', S, 'P', Pexp, 'written', written);
end

% ═════════════════════════════════════════════════════════════════════════
function P = local_empty_provenance(n)
    nanv = nan(n,1); strv = strings(n,1);
    P = struct('d_max_cm',nanv, 'frame_at_dmax',nanv, 'rebound_cm',nanv, ...
        'apg_peak_cm_s2',nanv, 'z_at_apg_peak_cm',nanv, 'n_frames_impact_to_stop',nanv, ...
        'dt_s',nanv, 'refMarkerID',nanv, 'kinMethod',strv, 'mmPerPx',nanv, 'g_cm_s2',nanv, ...
        'material',strv, 'container',strv, 'batch',strv, 'trialNum',nanv, 'heightLabel',strv, ...
        'nFrames',nanv, 'nTracked',nanv, 'firstValidFrame',nanv, 'windowStart',nanv, ...
        'windowEnd',nanv, 'autoWindow',false(n,1), 'eventFrameFirst',nanv, 'eventFrameLast',nanv, ...
        'gitCommit',strv, 'processedOn',strv, 'detectPass',strv, 'frameSource',strv, ...
        'rodAngle_peak_abs_deg',nanv, 'rodAngle_net_deg',nanv, 'rodAngle_range_deg',nanv, ...
        'rodAngle_noise_sd_deg',nanv, 'capture_datetime',strv, 'rawVideoPath',strv);
end

function v = local_field(s, name, default)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function v = local_num(T, name)
    v = NaN;
    if ismember(name, T.Properties.VariableNames)
        x = T.(name)(1);
        if isnumeric(x) || islogical(x), v = double(x); else, v = str2double(string(x)); end
    end
end

function local_write_series(path, S)
    fid = fopen(path, 'w');
    if fid < 0, error('export_master_dataset:csvFailed', 'Could not write %s', path); end
    fprintf(fid, 'trialTag,frame,t_s,z_cm,v_cm_s,a_cm_s2,a_plus_g_cm_s2,accel_se_cm_s2,winHalfFrames,rodBedDist_px\n');
    for i = 1:numel(S)
        if isempty(S(i).frame), continue; end
        tag = char(S(i).trialTag);
        blk = [S(i).frame, S(i).t_s, S(i).z_cm, S(i).v_cm_s, S(i).a_cm_s2, ...
               S(i).a_plus_g_cm_s2, S(i).accel_se_cm_s2, S(i).winHalfFrames, S(i).rodBedDist_px];
        for r = 1:size(blk,1)
            fprintf(fid, '%s,%d,%.6f,%.5f,%.4f,%.3f,%.3f,%.3f,%g,%.3f\n', tag, blk(r,1), blk(r,2:end));
        end
    end
    fclose(fid);
end

function D = local_dictionary()
    C = {
    % column                  units      source                          definition
    'trialTag'                ''         'file name'                     'Unique trial identifier; the capture label. Identifier only: never a physical axis.'
    'model'                   ''         'load_kinematics_set'           'Foot geometry: Tight / Default / Wide, from the tag suffix (unsuffixed = Default).'
    'condition'               ''         'kin_scalars.csv'               'material/container, e.g. GB/full.'
    'material'                ''         'meta'                          'GB (glass beads) or CHIN (chinchilla dust).'
    'container'               ''         'meta'                          'full / shallow / as_poured / dense.'
    'batch'                   ''         'meta.batchName'                'Capture batch folder the clip came from.'
    'trialNum'                ''         'meta'                          'Trial number within the height ladder.'
    'heightLabel'             ''         'meta'                          'Height label as captured, e.g. 25mm.'
    'dropHeight_mm'           'mm'       'load_kinematics_set'           'Release height, CORRECTED by get_relabel_map at read time. Replicate-group key; a control setting, not a measured input.'
    'dropHeight_asRecorded'   'mm'       'kin_scalars.csv'               'Height label as recorded on disk. Provenance only.'
    'relabeled'               'logical'  'get_relabel_map'               'true where dropHeight_mm came from the relabel map.'
    'isZeroDrop'              'logical'  'load_kinematics_set'           'dropHeight_mm == 0: released from contact, no impact event.'
    'v0_ff_cm_s'              'cm/s'     'derived'                       'sqrt(2 g h) from the corrected height. Cross-check only, never the fitted predictor (rail friction unmeasured).'
    'v0_freefall_ratio'       ''         'load_kinematics_set'           'v0_cm_s / v0_ff_cm_s. Healthy just under 1; Tight and Wide campaigns carry documented release-height offsets (README).'
    'v0_cm_s'                 'cm/s'     'kd_kinematics'                 'Measured impact speed: fitted velocity at the impact frame (velocity-first, KD 2007). THE driving variable for scaling. 0 on zero-drop rows by protocol.'
    'd_final_cm'              'cm'       'kd_kinematics'                 'Rod depth at the stop frame, zeroed at impact, positive into the bed. NaN on zero-drop rows (quarantine).'
    't_stop_s'                's'        'kd_kinematics'                 'Stop time from KD linear extrapolation of v(t) to zero, relative to impact. NaN on zero-drop rows.'
    'a_stop_cm_s2'            'cm/s^2'   'kd_kinematics'                 'Slope of the pre-stop v(t) fit: KD acceleration discontinuity at stop. NaN on zero-drop rows.'
    'd_max_cm'                'cm'       'derived from kin.z'            'Maximum depth over the retained record (impact to stop + postCap). Compare with d_final_cm: the difference is rebound.'
    'rebound_cm'              'cm'       'derived'                       'd_max_cm - d_final_cm. Rebound after the stop; expected small.'
    'apg_peak_cm_s2'          'cm/s^2'   'derived from kin.a_plus_g'     'Peak net grain deceleration g - a between impact and stop. Display-grade quantity (adaptive-window fit), not a fitted parameter.'
    'z_at_apg_peak_cm'        'cm'       'derived'                       'Depth at apg_peak. For the later two-bump / A(z) diagnostic only.'
    'n_frames_impact_to_stop' ''         'derived'                       'stopFrame - impactFrame + 1. Time resolution of the deceleration record.'
    'v0_meas_cm_s'            'cm/s'     'load_kinematics_set'           'Raw pipeline v0 before the zero-drop quarantine. QA only; NEVER fit.'
    'd_final_meas_cm'         'cm'       'load_kinematics_set'           'Raw pipeline d_final before quarantine. QA only; NEVER fit.'
    't_stop_meas_s'           's'        'load_kinematics_set'           'Raw pipeline t_stop before quarantine. QA only.'
    'a_stop_meas_cm_s2'       'cm/s^2'   'load_kinematics_set'           'Raw pipeline a_stop before quarantine. QA only.'
    'keep_reviewed'           'logical'  'load_kinematics_set defaults'  'true = in the reviewed analysis set (manual list, NaN-depth and glitch rules applied). THE selection column for analysis.'
    'exclude_reason'          ''         'derived'                       'manual / NaNDEPTH / GLITCH / empty. Why keep_reviewed is false.'
    'excl_manual'             'logical'  'get_manual_exclusions'         'On the reviewed manual-exclusion list (reason and date in that file).'
    'excl_nanDepth'           'logical'  'rule'                          'h > 0 trial with non-finite d_final_cm.'
    'excl_glitch'             'logical'  'rule'                          'h > 0 trial failing the monotonic-|v| test (tracking artefact).'
    'v0_guard'                ''         'load_kinematics_set'           'above_freefall (ratio > 1.15) / below_half_freefall (ratio < 0.5) / empty. Warns only; never drops.'
    'phi'                     ''         'get_substrate_properties'      'Packing fraction, five-preparation mean (0.624 GB/full). Random +/- only; see memory notes on systematics.'
    'phi_csv'                 ''         'kin_scalars.csv'               'phi as it stood at tracking time (retired single-measurement value). Provenance only.'
    'rho_particle_g_cm3'      'g/cm^3'   'get_substrate_properties'      'True particle density used as the phi reference (2.50 glass, 2.35 chinchilla).'
    'rho_bulk_g_cm3'          'g/cm^3'   'get_substrate_properties'      'Bulk density = phi * rho_particle.'
    'A_bare_cm2'              'cm^2'     'fig_foot_schematic'            'Projected foot area at the foot/bar junction (y = area_cutY_mm), holes excluded. Identical across models.'
    'A_hull_cm2'              'cm^2'     'fig_foot_schematic'            'Convex-hull area of the same footprint. Grows with splay: 2.607 / 3.495 / 4.052.'
    'hull_over_bare'          ''         'derived'                       'A_hull_cm2 / A_bare_cm2.'
    'A_bare_at_dfinal_cm2'    'cm^2'     'foot_area_vs_depth.csv'        'Cumulative projected area of EVERYTHING below the surface at z = d_final (bar and post included). Interpolated; NaN beyond 4 cm.'
    'A_hull_at_dfinal_cm2'    'cm^2'     'foot_area_vs_depth.csv'        'Convex hull of the same at z = d_final.'
    'a_local_at_dfinal_cm2'   'cm^2'     'foot_area_vs_depth.csv'        'Cross-section of the model at the surface plane when the toe tip is at d_final.'
    'area_cutY_mm'            'mm'       'fig_foot_schematic'            'STL cut plane the constants were computed at (-50 = foot/bar junction = z 1.119 cm).'
    'mass_g'                  'g'        'option MassG'                  'Total projectile mass (65 g, as track_tracers_2 default).'
    'fps'                     '1/s'      'kin_scalars.csv'               'True frame rate used for the time base (post fps repair).'
    'dt_s'                    's'        'kin.dt'                        '1/fps as used by kd_kinematics.'
    'impactFrame'             ''         'kd_kinematics'                 'Impact frame index (velocity peak near the geometric bed crossing).'
    'stopFrame'               ''         'kd_kinematics'                 'Stop frame index (nearest frame to t_stop).'
    'frame_at_dmax'           ''         'derived'                       'Frame index of d_max_cm.'
    'nFrames'                 ''         'meta'                          'Frames in the tracked window.'
    'nTracked'                ''         'meta'                          'Frames with tracked markers.'
    'firstValidFrame'         ''         'meta'                          'First frame with exactly N markers detected (tracking frame 1).'
    'windowStart'             ''         'meta'                          'First raw frame of the processed window.'
    'windowEnd'               ''         'meta'                          'Last raw frame of the processed window.'
    'autoWindow'              'logical'  'meta'                          'Window chosen automatically.'
    'eventFrameFirst'         ''         'meta.eventFrameRange'          'First raw frame of the saved event-frame subset.'
    'eventFrameLast'          ''         'meta.eventFrameRange'          'Last raw frame of the saved event-frame subset.'
    'mmPerPx'                 'mm/px'    'calib'                         'Imaging calibration (0.1079 verified).'
    'impactDistPx'            'px'       'calib'                         'Reference-marker distance to the bed line at contact (-360 unified).'
    'bedX'                    'px'       'calib'                         'Bed line reference.'
    'g_cm_s2'                 'cm/s^2'   'calib'                         'g used by the pipeline (980).'
    'refMarkerID'             ''         'kd_kinematics'                 'Marker used as the impact reference (top marker on the rod).'
    'rodAngle_peak_abs_deg'   'deg'      'rod_angle'                     'Peak |rod tilt| over the impact-stop interval. Obliquity QA.'
    'rodAngle_net_deg'        'deg'      'rod_angle'                     'Net tilt change impact to stop.'
    'rodAngle_range_deg'      'deg'      'rod_angle'                     'Tilt range over the interval.'
    'rodAngle_noise_sd_deg'   'deg'      'rod_angle'                     'Tilt noise SD (pre-impact).'
    'gitCommit'               ''         'meta.provenance'               'Pipeline commit that produced the tracks.'
    'processedOn'             ''         'meta.provenance'               'Tracking date/time.'
    'detectPass'              ''         'meta.provenance'               'pass1 / pass2-backup.'
    'frameSource'             ''         'meta.provenance'               'stream / frames.'
    'kinMethod'               ''         'kin.method'                    'kd_kinematics method tag.'
    'capture_datetime'        ''         'raw .avi timestamp'            'File time of the raw clip (RawRoot given). Session provenance for the Default 325/345/365 split-bin question.'
    'rawVideoPath'            ''         'find_raw_video'                'Raw clip located for capture_datetime.'
    'kinPath'                 ''         'load_kinematics_set'           'Path of the _kin.mat this row was read from.'
    };
    D = cell2table(C, 'VariableNames', {'column', 'units', 'source', 'definition'});
end

function s = local_tern(c, a, b), if c, s = a; else, s = b; end, end
