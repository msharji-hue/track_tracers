function status = process_one_trial(cfg)
% PROCESS_ONE_TRIAL  Tracking-only pipeline for ONE video.
%
%   export frames -> detect circles -> find first frame with exactly N markers
%   -> redefine that as tracking frame 1 -> track_markers -> save + QA.
%
%   Deliberately contains NO kinematics: no Savitzky-Golay smoothing, no
%   velocity/acceleration, no (a+g), no force-law fitting, no rod rotation, and
%   no impact/stop event detection. Those live downstream (track_tracers_2.m ->
%   kd_kinematics.m, rod_angle.m) and operate on the saved tracks. This keeps a
%   smoothing or model problem from ever costing a 368-video batch.
%
%   Status semantics:
%     OK      : detections + tracks + metadata + QA all saved
%     PARTIAL : detections + tracks + metadata saved, but QA failed
%     FAILED  : no usable tracks (e.g. no frame with exactly N markers)
%
%   cfg fields: videoPath, framesDir, detDir, resultsDir, meta, params, calib,
%   nExpectedMarkers(8), interactive, makeQA, reuseFrames, reuseDetections,
%   filterType('sharpen').

    if ~isfield(cfg,'nExpectedMarkers'), cfg.nExpectedMarkers = 8;        end
    if ~isfield(cfg,'filterType'),       cfg.filterType       = 'sharpen';end
    if ~isfield(cfg,'makeQA'),           cfg.makeQA           = true;     end
    if ~isfield(cfg,'reuseFrames'),      cfg.reuseFrames      = false;    end
    if ~isfield(cfg,'reuseDetections'),  cfg.reuseDetections  = false;    end

    m = cfg.meta;

    % ── substrate bed properties (per material+condition; metadata only) ──
    sub = get_substrate_properties(m.material, m.container);
    m.rho_particle_g_cm3 = sub.rho_particle_g_cm3;
    m.rho_bulk_g_cm3     = sub.rho_bulk_g_cm3;
    m.phi                = sub.phi;
    if ~sub.ok
        fprintf('  NOTE: no substrate properties for (%s, %s) — phi/rho stored as NaN\n', ...
                m.material, m.container);
    end

    status = struct('ok',false, 'partial',false, 'tracksSaved',false, ...
                    'reason','', 'firstValidFrame',NaN, 'fps',NaN, ...
                    'nFrames',NaN, 'nTracked',NaN, ...
                    'framesDir',cfg.framesDir, 'detDir',cfg.detDir, ...
                    'resultsDir',cfg.resultsDir, 'meta',m);

    % ── 1) OPEN VIDEO ─────────────────────────────────────────────────────
    v   = open_video(cfg.videoPath);
    fps = v.FrameRate;
    status.fps = fps;
    m.fps_true = fps;

    % ── 2) EXPORT FRAMES ──────────────────────────────────────────────────
    %   The exported set is a WINDOW of the video, not necessarily the whole of
    %   it (see opts.autoWindow in process_trial). Every index downstream --
    %   detections, firstValidFrame, tracking frame 1 -- is relative to this
    %   window, so its absolute start is recorded on the meta as windowStart.
    %   Absolute video frame = windowStart + firstValidFrame + k - 2 for
    %   tracking frame k.
    haveFrames = isfolder(cfg.framesDir) && ...
                 ~isempty(dir(fullfile(cfg.framesDir, '*.png')));
    if cfg.reuseFrames && haveFrames
        fprintf('Reusing existing frames: %s\n', cfg.framesDir);
        % Take the window from the filenames ON DISK, not from cfg: reused
        % frames may have been exported under a different window, and trusting
        % cfg here would silently shift every absolute frame index.
        [startFrame, endFrame] = local_window_from_frames(cfg.framesDir);
        if isnan(startFrame)
            startFrame = 1;
            endFrame   = floor(v.Duration * v.FrameRate);
        end
        fprintf('  reused window: frames %d-%d\n', startFrame, endFrame);
    else
        fprintf('Exporting frames...\n');
        if ~exist(cfg.framesDir, 'dir'), mkdir(cfg.framesDir); end
        nAvail = floor(v.Duration * v.FrameRate);
        if cfg.interactive
            startFrame = pick_start_frame(v);
            endFrame   = pick_end_frame(v, startFrame);
        elseif isfield(cfg,'windowStart') && ~isempty(cfg.windowStart) && ...
               isfield(cfg,'windowEnd')   && ~isempty(cfg.windowEnd)
            startFrame = max(1, round(cfg.windowStart));
            endFrame   = min(nAvail, round(cfg.windowEnd));
        else
            startFrame = 1;
            endFrame   = nAvail;
        end
        export_frames(v, startFrame, endFrame, cfg.framesDir, cfg.filterType);
        fps_true = fps; %#ok<NASGU>
        save(fullfile(cfg.framesDir, 'video_meta.mat'), 'fps_true');
    end
    m.windowStart = startFrame;
    m.windowEnd   = endFrame;
    m.autoWindow  = isfield(cfg,'autoWindow') && isequal(cfg.autoWindow, true);
    status.windowStart = startFrame;
    status.windowEnd   = endFrame;

    % ── 3) DETECT CIRCLES ─────────────────────────────────────────────────
    detMat = fullfile(cfg.detDir, [m.trialTag '_detections.mat']);
    if cfg.reuseDetections && exist(detMat, 'file')
        fprintf('Reusing existing detections: %s\n', detMat);
        S           = load(detMat);
        detectOut   = S.det.detect;
        centersCell = detectOut.centersCell;
        nDetected   = detectOut.nDetected;
    else
        fprintf('Detecting circles...\n');
        params = cfg.params;
        params.heightLabel = m.heightLabel;
        if ~cfg.interactive, params.showPreviewEveryN = 0; end
        detectOut = detect_circles_per_frame(cfg.framesDir, params);

        if ~exist(cfg.detDir, 'dir'), mkdir(cfg.detDir); end
        saveInfo = struct( ...
            'material',      m.material,      'batchName',   m.batchName, ...
            'heightLabel',   m.heightLabel,   'trialNum',    m.trialNum, ...
            'container',     m.container,     'dropHeight_mm', m.dropHeight_mm, ...
            'trialTag',      m.trialTag,      'framesDir',   cfg.framesDir, ...
            'fps_export',    fps,             'detDir',      cfg.detDir);
        save_detections(detectOut, saveInfo);

        centersCell = detectOut.centersCell;
        nDetected   = detectOut.nDetected;
    end
    status.nFrames = numel(centersCell);

    % ── 4) FIRST FRAME WITH EXACTLY N MARKERS ─────────────────────────────
    firstValidFrame = find(nDetected(:) == cfg.nExpectedMarkers, 1, 'first');
    if isempty(firstValidFrame)
        status.reason = sprintf('no frame with exactly %d detected markers', ...
                                cfg.nExpectedMarkers);
        write_trial_log(cfg.resultsDir, m, status);
        return;                                  % FAILED
    end
    status.firstValidFrame = firstValidFrame;
    m.firstValidFrame      = firstValidFrame;
    fprintf('Earliest valid %d-marker frame found at frame %d\n', ...
            cfg.nExpectedMarkers, firstValidFrame);

    % ── 5) TRACK (tracking frame 1 := firstValidFrame) ───────────────────
    fprintf('Tracking marker positions...\n');
    detForTrack          = centersCell(firstValidFrame:end);
    firstFrameCenters    = detForTrack{1};
    % Tolerance must stay below one inter-marker spacing (~34.6 px) so a stale
    % last-known position cannot match a neighbouring marker's detection.
    if isfield(cfg.calib,'trackTolerancePx')
        tolPx = cfg.calib.trackTolerancePx;
    else
        tolPx = 25;   % fallback for older calibration structs
    end
    [trackedX, trackedY] = track_markers(detForTrack, firstFrameCenters, tolPx);
    status.nTracked      = size(trackedX, 2);
    m.nFrames            = status.nFrames;
    m.nTracked           = status.nTracked;

    % ── 6) SAVE TRACKS + METADATA ────────────────────────────────────────
    fprintf('Saving results...\n');
    tracks = struct('trackedX',trackedX, 'trackedY',trackedY, ...
                    'firstValidFrame',firstValidFrame, 'fps',fps, ...
                    'nFrames',status.nFrames, 'nTracked',status.nTracked);
    save_tracks(m, tracks, cfg.calib, cfg.resultsDir);
    status.tracksSaved = true;

    % ── 7) QA (never downgrades a good trial to FAILED) ──────────────────
    if cfg.makeQA
        try
            make_track_qa(m, tracks, nDetected, cfg.resultsDir);
            annotDir = fullfile(cfg.resultsDir, 'qa', [m.trialTag '_annotated']);
            make_annotated_frames(cfg.framesDir, detectOut, trackedX, trackedY, ...
                                  firstValidFrame, cfg.calib, annotDir, 'maxFrames', 300);
        catch MEq
            status.partial = true;
            status.reason  = ['QA failed: ' MEq.message];
            fprintf('  WARNING: %s (detections + tracks saved)\n', status.reason);
            write_trial_log(cfg.resultsDir, m, status);
            return;                              % PARTIAL
        end
    end

    status.ok = true;
    write_trial_log(cfg.resultsDir, m, status);

    if cfg.interactive
        assignin('base','trackedX',trackedX);
        assignin('base','trackedY',trackedY);
        assignin('base','tracks',  tracks);
    end
    fprintf('Done.\n');
end

% ─────────────────────────────────────────────────────────────────────────
function write_trial_log(resultsDir, m, status)
    logDir = fullfile(resultsDir, 'logs');
    if ~exist(logDir, 'dir'), mkdir(logDir); end
    fid = fopen(fullfile(logDir, [m.trialTag '_log.txt']), 'w');

    fprintf(fid, 'trial       : %s\n', m.trialTag);
    fprintf(fid, 'material    : %s\n', m.material);
    fprintf(fid, 'batch       : %s\n', m.batchName);
    fprintf(fid, 'dropHeight  : %g mm\n', m.dropHeight_mm);
    fprintf(fid, 'trialNum    : %d\n', m.trialNum);
    fprintf(fid, 'condition   : %s\n', m.container);
    if isfield(m,'phi')
        fprintf(fid, 'substrate   : rho_p=%.2f  rho_bulk=%.3f  phi=%.3f (g/cm3)\n', ...
                m.rho_particle_g_cm3, m.rho_bulk_g_cm3, m.phi);
    end
    if status.ok,          st = 'OK';
    elseif status.partial, st = 'PARTIAL (tracks saved, QA failed)';
    else,                  st = 'FAILED';
    end
    fprintf(fid, 'status      : %s\n', st);
    fprintf(fid, 'fps         : %s\n', fmt(status.fps, '%.4f'));
    fprintf(fid, 'nFrames     : %s\n', num2str(status.nFrames));
    fprintf(fid, 'firstValid  : %s\n', num2str(status.firstValidFrame));
    fprintf(fid, 'nTracked    : %s\n', num2str(status.nTracked));
    fprintf(fid, 'tracksSaved : %d\n', status.tracksSaved);
    if ~isempty(status.reason)
        fprintf(fid, 'reason      : %s\n', status.reason);
    end
    fprintf(fid, '\nNOTE: kinematics (v, a, a+g), force-law fits, event detection\n');
    fprintf(fid, '      (impact/stop), and rod bending are computed DOWNSTREAM\n');
    fprintf(fid, '      from the saved tracks; this pipeline is tracking-only.\n');
    fclose(fid);
end

function s = fmt(x, f)
    if isnan(x), s = 'NaN'; else, s = sprintf(f, x); end
end

function [lo, hi] = local_window_from_frames(framesDir)
%LOCAL_WINDOW_FROM_FRAMES  Absolute frame range of the PNGs already on disk.
%   export_frames names every file by its ABSOLUTE video frame index
%   (frame_%05d.png), so the exported window can be recovered from the
%   filenames alone. Returns NaN if nothing parses.
    lo = NaN; hi = NaN;
    d = dir(fullfile(framesDir, 'frame_*.png'));
    if isempty(d), return; end
    n = nan(numel(d),1);
    for i = 1:numel(d)
        t = regexp(d(i).name, 'frame_(\d+)\.png$', 'tokens', 'once');
        if ~isempty(t), n(i) = str2double(t{1}); end
    end
    n = n(isfinite(n));
    if isempty(n), return; end
    lo = min(n); hi = max(n);
end
