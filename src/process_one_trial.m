function status = process_one_trial(cfg)
% PROCESS_ONE_TRIAL  Full per-trial pipeline for ONE video.
%
%   export frames -> detect circles -> find first frame with exactly N markers
%   -> redefine that as tracking frame 1 -> track_markers -> SAVE RAW TRACKS
%   -> kinematics -> clip at t_stop + buffer -> save full bundle / QA / log.
%
%   Robustness:
%     * No 8-marker frame        -> clean FAILED (tracks not saved), caller continues.
%     * Tracking succeeds but     -> RAW positions are saved immediately, so a
%       kinematics/t_stop fails      kinematics failure returns PARTIAL with the
%                                    tracked positions preserved (not a total loss).
%
%   cfg fields: videoPath, framesDir, detDir, resultsDir, meta, params, calib,
%   nExpectedMarkers(8), interactive, makeFigures, reuseFrames, reuseDetections,
%   filterType('sharpen').
%
%   status fields: ok, partial, tracksSaved, reason, firstValidFrame, fps,
%   v0_cm_s, d_final_cm, t_stop_s, framesDir, detDir, resultsDir, meta.

    if ~isfield(cfg,'nExpectedMarkers'), cfg.nExpectedMarkers = 8;        end
    if ~isfield(cfg,'filterType'),       cfg.filterType       = 'sharpen';end
    if ~isfield(cfg,'makeFigures'),      cfg.makeFigures      = false;    end
    if ~isfield(cfg,'reuseFrames'),      cfg.reuseFrames      = false;    end
    if ~isfield(cfg,'reuseDetections'),  cfg.reuseDetections  = false;    end

    m = cfg.meta;

    % ── substrate bed properties (per material+condition; additive metadata) ─
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
                    'v0_cm_s',NaN, 'd_final_cm',NaN, 't_stop_s',NaN, ...
                    'framesDir',cfg.framesDir, 'detDir',cfg.detDir, ...
                    'resultsDir',cfg.resultsDir, 'meta',m);

    % ── 1) OPEN VIDEO ─────────────────────────────────────────────────────
    v   = open_video(cfg.videoPath);
    fps = v.FrameRate;
    status.fps = fps;
    m.fps_true = fps;

    % ── 2) EXPORT FRAMES (frame 1 .. last) ───────────────────────────────
    haveFrames = isfolder(cfg.framesDir) && ...
                 ~isempty(dir(fullfile(cfg.framesDir, '*.png')));
    if cfg.reuseFrames && haveFrames
        fprintf('Reusing existing frames: %s\n', cfg.framesDir);
    else
        fprintf('Exporting frames...\n');
        if ~exist(cfg.framesDir, 'dir'), mkdir(cfg.framesDir); end
        if cfg.interactive
            startFrame = pick_start_frame(v);
            endFrame   = pick_end_frame(v, startFrame);
        else
            startFrame = 1;
            endFrame   = floor(v.Duration * v.FrameRate);
        end
        export_frames(v, startFrame, endFrame, cfg.framesDir, cfg.filterType);
        fps_true = fps; %#ok<NASGU>
        save(fullfile(cfg.framesDir, 'video_meta.mat'), 'fps_true');
    end

    % ── 3) DETECT CIRCLES ─────────────────────────────────────────────────
    detMat = fullfile(cfg.detDir, [m.trialTag '_detections.mat']);
    if cfg.reuseDetections && exist(detMat, 'file')
        fprintf('Reusing existing detections: %s\n', detMat);
        S           = load(detMat);
        centersCell = S.det.detect.centersCell;
        nDetected   = S.det.detect.nDetected;
    else
        fprintf('Detecting circles...\n');
        params = cfg.params;
        params.heightLabel = m.heightLabel;
        if ~cfg.interactive, params.showPreviewEveryN = 0; end
        detectOut = detect_circles_per_frame(cfg.framesDir, params);

        if ~exist(cfg.detDir, 'dir'), mkdir(cfg.detDir); end
        saveInfo = struct( ...
            'material',      m.material, ...
            'batchName',     m.batchName, ...
            'heightLabel',   m.heightLabel, ...
            'trialNum',      m.trialNum, ...
            'container',     m.container, ...
            'dropHeight_mm', m.dropHeight_mm, ...
            'trialTag',      m.trialTag, ...
            'framesDir',     cfg.framesDir, ...
            'fps_export',    fps, ...
            'detDir',        cfg.detDir);
        save_detections(detectOut, saveInfo);

        centersCell = detectOut.centersCell;
        nDetected   = detectOut.nDetected;
    end

    % ── 4) FIRST FRAME WITH EXACTLY N MARKERS ─────────────────────────────
    firstValidFrame = find(nDetected(:) == cfg.nExpectedMarkers, 1, 'first');
    if isempty(firstValidFrame)
        status.reason = sprintf('no frame with exactly %d detected markers', ...
                                cfg.nExpectedMarkers);
        write_trial_log(cfg.resultsDir, m, status, []);
        return;
    end
    status.firstValidFrame = firstValidFrame;
    m.firstValidFrame      = firstValidFrame;
    fprintf('Earliest valid %d-marker frame found at frame %d\n', ...
            cfg.nExpectedMarkers, firstValidFrame);

    % ── 5) TRACK (frame 1 := firstValidFrame) + SAVE RAW IMMEDIATELY ──────
    fprintf('Tracking marker positions...\n');
    detForTrack          = centersCell(firstValidFrame:end);
    firstFrameCenters    = detForTrack{1};
    [trackedX, trackedY] = track_markers(detForTrack, firstFrameCenters, 100);

    trackDir = fullfile(cfg.resultsDir, 'tracks');
    if ~exist(trackDir, 'dir'), mkdir(trackDir); end
    rawMeta = m; %#ok<NASGU>
    save(fullfile(trackDir, [m.trialTag '_tracks_raw.mat']), ...
        'trackedX', 'trackedY', 'firstValidFrame', 'fps', 'rawMeta', '-v7.3');
    status.tracksSaved = true;

    % ── 6) KINEMATICS + STOP (+ buffer) ─ guarded so failure => PARTIAL ──
    dt            = 1 / fps;
    postCapFrames = round(fps * (cfg.calib.postCapMs / 1000));
    c             = cfg.calib;
    try
        [t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, impact_index, ...
         toeMarkerID, toePx, stopFrame, sgOrder, sgWindow] = ...
            toe_kinematics(trackedX, trackedY, c.lineA, c.lineB, c.lineC, ...
                           dt, c.mmPerPx, c.impactDistPx, ...
                           'postCapFrames', postCapFrames);
    catch ME
        status.partial = true;
        status.reason  = ['kinematics failed: ' ME.message];
        fprintf('  WARNING: %s (raw tracks saved)\n', status.reason);
        write_trial_log(cfg.resultsDir, m, status, []);
        return;
    end

    a_plus_g   = -a_smooth - c.g_cm_s2;
    v0_cm_s    = v_smooth(impact_index);
    d_final_cm = z_smooth(stopFrame);
    t_stop_s   = t_s(stopFrame);
    status.v0_cm_s    = v0_cm_s;
    status.d_final_cm = d_final_cm;
    status.t_stop_s   = t_stop_s;

    savedThroughFrame = min(stopFrame + postCapFrames, size(trackedX, 2));

    % ── 7) BUNDLE + SAVE ─────────────────────────────────────────────────
    fprintf('Saving results...\n');
    tracks = struct( ...
        'trackedX',          trackedX,          'trackedY',     trackedY,     ...
        't_s',               t_s,               'depthRod_cm',  depthRod_cm,  ...
        'z_smooth',          z_smooth,          'v_smooth',     v_smooth,     ...
        'a_smooth',          a_smooth,          'a_plus_g',     a_plus_g,     ...
        'impact_index',      impact_index,      'stopFrame',    stopFrame,    ...
        'savedThroughFrame', savedThroughFrame, 'postCapFrames',postCapFrames,...
        'toeMarkerID',       toeMarkerID,       'toePx',        toePx,        ...
        'sgOrder',           sgOrder,           'sgWindow',     sgWindow,     ...
        'firstValidFrame',   firstValidFrame);
    scalars = struct('v0_cm_s',v0_cm_s, 'd_final_cm',d_final_cm, 't_stop_s',t_stop_s);
    scalars.rho_particle_g_cm3 = m.rho_particle_g_cm3;
    scalars.rho_bulk_g_cm3     = m.rho_bulk_g_cm3;
    scalars.phi                = m.phi;

    % ── rod-bending diagnostic (additive QA metric; never fails the trial) ─
    status.bend = [];
    try
        bend = rod_bending(trackedX, trackedY, impact_index, stopFrame, ...
                           t_s, cfg.calib.mmPerPx);
        tracks.bending             = bend;
        scalars.bend_peak_rms_mm   = bend.peak_rms_mm;
        scalars.bend_peak_max_mm   = bend.peak_max_mm;
        scalars.bend_signed_pk_mm  = bend.signed_peak_mm;
        scalars.bend_at_stop_mm    = bend.bend_at_stop_mm;
        scalars.bend_curv_pk_1pmm  = bend.curv_peak_1pmm;
        scalars.bend_tilt_pk_deg   = bend.tilt_peak_deg;
        scalars.bend_angle_pk_deg  = bend.bend_angle_peak_deg;
        scalars.bend_seg_pk_deg    = bend.seg_angle_peak_deg;
        scalars.bend_t_peak_ms     = bend.t_peak_ms;
        scalars.bend_baseline_mm   = bend.baseline_rms_mm;
        scalars.bend_flag          = double(bend.bendFlag);
        scalars.tilt_flag          = double(bend.tiltFlag);
        status.bend = bend;
        if bend.bendFlag
            fprintf('  NOTE: rod bending flagged (peak RMS %.3f mm vs baseline %.3f mm)\n', ...
                    bend.peak_rms_mm, bend.baseline_rms_mm);
        end
    catch MEb
        fprintf('  WARNING: rod-bending diagnostic failed: %s\n', MEb.message);
        tracks.bending = struct('error', MEb.message);
    end

    save_tracks(m, tracks, scalars, cfg.calib, cfg.resultsDir);

    status.ok = true;
    write_trial_log(cfg.resultsDir, m, status, ...
                    struct('impact_index',impact_index,'stopFrame',stopFrame));

    if cfg.makeFigures
        make_track_qa(m, tracks, nDetected, cfg.resultsDir);
    end
    if cfg.interactive
        assignin('base','trackedX',trackedX);
        assignin('base','trackedY',trackedY);
        assignin('base','tracks',  tracks);
        assignin('base','scalars', scalars);
    end
    fprintf('Done.\n');
end

% ─────────────────────────────────────────────────────────────────────────
function write_trial_log(resultsDir, m, status, extra)
    logDir = fullfile(resultsDir, 'logs');
    if ~exist(logDir, 'dir'), mkdir(logDir); end
    fid = fopen(fullfile(logDir, [m.trialTag '_log.txt']), 'w');
    fprintf(fid, 'trial       : %s\n', m.trialTag);
    fprintf(fid, 'material    : %s\n', m.material);
    fprintf(fid, 'batch       : %s\n', m.batchName);
    fprintf(fid, 'dropHeight  : %g mm\n', m.dropHeight_mm);
    fprintf(fid, 'trialNum    : %d\n', m.trialNum);
    fprintf(fid, 'container   : %s\n', m.container);
    if isfield(m,'phi')
        fprintf(fid, 'substrate   : rho_p=%.2f  rho_bulk=%.3f  phi=%.3f (g/cm3)\n', ...
                m.rho_particle_g_cm3, m.rho_bulk_g_cm3, m.phi);
    end
    if status.ok,        st = 'OK';
    elseif status.partial, st = 'PARTIAL (tracks saved, kinematics failed)';
    else,                st = 'FAILED';
    end
    fprintf(fid, 'status      : %s\n', st);
    fprintf(fid, 'fps         : %s\n', fmt(status.fps, '%.4f'));
    fprintf(fid, 'firstValid  : %s\n', num2str(status.firstValidFrame));
    fprintf(fid, 'tracksSaved : %d\n', status.tracksSaved);
    if ~isempty(status.reason), fprintf(fid, 'reason      : %s\n', status.reason); end
    if status.ok && ~isempty(extra)
        fprintf(fid, 'impact_idx  : %d\n', extra.impact_index);
        fprintf(fid, 'stopFrame   : %d\n', extra.stopFrame);
        fprintf(fid, 'v0_cm_s     : %.4f\n', status.v0_cm_s);
        fprintf(fid, 'd_final_cm  : %.4f\n', status.d_final_cm);
        fprintf(fid, 't_stop_s    : %.6f\n', status.t_stop_s);
    end
    if isfield(status,'bend') && ~isempty(status.bend)
        b = status.bend;
        fprintf(fid, 'bend_peakRMS: %.4f mm  (max %.4f, signed %.4f)\n', ...
                b.peak_rms_mm, b.peak_max_mm, b.signed_peak_mm);
        fprintf(fid, 'bend_angle  : %.3f deg   seg_max: %.3f deg\n', ...
                b.bend_angle_peak_deg, b.seg_angle_peak_deg);
        fprintf(fid, 'rod_tilt    : %.3f deg\n', b.tilt_peak_deg);
        fprintf(fid, 'bend_base   : %.4f mm   t_peak: %.2f ms   bendFLAG: %d   tiltFLAG: %d\n', ...
                b.baseline_rms_mm, b.t_peak_ms, b.bendFlag, b.tiltFlag);
    end
    fclose(fid);
end

function s = fmt(x, f)
    if isnan(x), s = 'NaN'; else, s = sprintf(f, x); end
end
