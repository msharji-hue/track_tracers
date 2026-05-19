%% scripts/process_trial.m
% Unified trial processing pipeline.
%   
% Combines:
%   1) Crop & Export  X— pick start/end frames, export PNGs, auto-save fps
%   2) Select & Detect — detect circular markers from exported PNGs
%
% Outputs are kept separate:
%   Crop step   → user-named frames folder  (PNG frames + video_meta.mat)
%   Detect step → user-named results folder (detections .mat + .csv)
%
% fps is auto-saved to video_meta.mat during crop and loaded
% automatically during detect — no manual entry needed.

clear; close all; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

% ══════════════════════════════════════════════════════════════════════════
%  STEP SELECTION
% ══════════════════════════════════════════════════════════════════════════
stepChoice = menu('Select processing step', ...
    '1 — Crop & Export  (video → PNG frames)', ...
    '2 — Select & Detect  (PNG frames → detections)', ...
    '3 — Full pipeline  (crop then detect)');

if stepChoice == 0
    fprintf('Cancelled.\n'); return;
end

runCrop   = ismember(stepChoice, [1 3]);
runDetect = ismember(stepChoice, [2 3]);

% ══════════════════════════════════════════════════════════════════════════
%  STEP 1 — CROP & EXPORT
% ══════════════════════════════════════════════════════════════════════════
if runCrop
    fprintf('\n═══ CROP & EXPORT ═══\n');

    % ── Select video ──────────────────────────────────────────────────────
    inputDir  = uigetdir(pwd, 'Select folder containing your video');
    [file, ~] = uigetfile({'*.mp4;*.avi;*.mov;*.MP4', 'Video Files'}, ...
        'Select a video', inputDir);
    if isequal(file, 0), fprintf('Cancelled.\n'); return; end
    videoPath = fullfile(inputDir, file);

    % ── Select output ─────────────────────────────────────────────────────
    parentDir  = uigetdir(pwd, 'Select where to save frames folder');
    folderName = inputdlg('Name your frames folder:', 'Frames Folder', 1, {'frames_output'});
    if isempty(folderName), fprintf('Cancelled.\n'); return; end
    framesDir  = fullfile(parentDir, folderName{1});

    % ── Filter ────────────────────────────────────────────────────────────
    filterChoice = menu('Frame filter', 'sharpen', 'none', 'gaussian', 'grayscale');
    filterTypes  = {'sharpen', 'none', 'gaussian', 'grayscale'};
    filterType   = filterTypes{max(1, filterChoice)};

    % ── Open video + pick frames ───────────────────────────────────────────
    v          = open_video(videoPath);
    fps_true   = v.FrameRate;
    startFrame = pick_start_frame(v);
    endFrame   = pick_end_frame(v, startFrame);

    % ── Export frames ─────────────────────────────────────────────────────
    export_frames(v, startFrame, endFrame, framesDir, filterType);

    % ── Auto-save fps metadata ────────────────────────────────────────────
    save(fullfile(framesDir, 'video_meta.mat'), 'fps_true');
    fprintf('Saved video_meta.mat  →  fps = %.4f\n', fps_true);
    fprintf('Frames saved to: %s\n', framesDir);
end

% ══════════════════════════════════════════════════════════════════════════
%  STEP 2 — SELECT & DETECT
% ══════════════════════════════════════════════════════════════════════════
if runDetect
    fprintf('\n═══ SELECT & DETECT ═══\n');

    % ── Select frames folder ──────────────────────────────────────────────
    if runCrop
        fprintf('Using frames folder from crop step: %s\n', framesDir);
    else
        framesDir = uigetdir(pwd, 'Select folder containing PNG frames');
        if isequal(framesDir, 0), fprintf('Cancelled.\n'); return; end
    end

    % ── Load fps ──────────────────────────────────────────────────────────
    metaFile = fullfile(framesDir, 'video_meta.mat');
    if exist(metaFile, 'file')
        meta     = load(metaFile);
        fps_true = meta.fps_true;
        fprintf('Loaded FPS from metadata: %.4f\n', fps_true);
    else
        fps_true = 2250;
        warning('No video_meta.mat found in frames folder — using default fps = %.0f', fps_true);
    end

    % ── Select results folder ─────────────────────────────────────────────
    resultsDir = uigetdir(pwd, 'Select folder to save detection results');
    if isequal(resultsDir, 0), fprintf('Cancelled.\n'); return; end
    folderName = inputdlg('Name your results subfolder:', 'Output', 1, {'LDH_X_X'});
    if isempty(folderName), fprintf('Cancelled.\n'); return; end
    detDir = fullfile(resultsDir, folderName{1});
    if ~exist(detDir, 'dir'), mkdir(detDir); end

    % ── Trial info ────────────────────────────────────────────────────────
    fields = inputdlg( ...
        {'Material (GB or CHIN)', 'Batch', 'Height label', 'Trial number'}, ...
        'Trial Info', 1, {'GB', 'Batch 1', 'H10', '1'});
    if isempty(fields), fprintf('Cancelled.\n'); return; end

    trialInfo = struct( ...
        'material',    fields{1},             ...
        'batchName',   fields{2},             ...
        'heightLabel', fields{3},             ...
        'trialNum',    str2double(fields{4}), ...
        'fps_true',    fps_true);

    % ── Detection parameters ──────────────────────────────────────────────
    params = struct();
    params.radiusRange       = [9 33];
    params.sensitivity       = 0.85;
    params.edgeThresh        = 0.10;
    params.polarity          = 'bright';
    params.alphaG            = 0.50;
    params.betaB             = 0.50;
    params.doCLAHE           = true;
    params.medianK           = 3;
    params.showPreviewEveryN = 100;
    params.showCenters       = true;
    params.heightLabel       = trialInfo.heightLabel;

    % ── Run detection ─────────────────────────────────────────────────────
    detectOut = detect_circles_per_frame(framesDir, params);

    saveInfo = struct( ...
        'material',    trialInfo.material,    ...
        'batchName',   trialInfo.batchName,   ...
        'heightLabel', trialInfo.heightLabel, ...
        'trialNum',    trialInfo.trialNum,    ...
        'framesDir',   framesDir,             ...
        'fps_export',  trialInfo.fps_true,    ...
        'detDir',      detDir);

    det = save_detections(detectOut, saveInfo);
    fprintf('Detections saved to: %s\n', detDir);
end
