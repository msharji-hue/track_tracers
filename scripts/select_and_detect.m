%% select_and_detect.m
clear; close all; clc;

%% 1) USER INPUTS

% Trial identifiers
material    = 'GB';          % "GB" or "CHIN"
batchName   = 'Batch1';      % edit per batch number
heightLabel = 'H10';
trialNum    = 2;

% Clip type and file name (exported ProRes)
videoFile   = sprintf('%s_T%02d_%s.mov', heightLabel, trialNum);
fps_true     = 2715;         % true camera acquisition fps

% Detection settings
radiusRange = [9 35];        % pixels, tune later
sensitivity = 0.84;          % 0-1
edgeThresh  = 0.10;          % 0-1
polarity    = 'bright';      % dots should appear bright in redness image

% Redness image weights
alphaG = 0.50;               % A = R - alphaG*G - betaB*B
betaB  = 0.50;

% Minimal preprocessing
doCLAHE = true;
medianK = 3;                 % 1 disables; 3 is usually safe

% Visualization for debugging
showPreviewEveryN = 100;     % 0 disables; try 100 to preview every 100 frames
showCenters       = true;

%% 2) PROJECT PATHS

% Data root (videos, results, logs, notes)
projectRoot = ['/Users/muhannadalsharji/Library/CloudStorage/GoogleDrive-msharji@umich.edu/' ...
               'My Drive/Research /Videos/SP26/JerboaImpact_VideoExports'];

% Code root (Git repo, separate from data)
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

% Videos
rawDir = fullfile(projectRoot, '01_RAW', material, batchName);
exportsDir = fullfile(projectRoot, '02_EXPORTS', material, batchName, heightLabel);

% Results
resultsRoot = fullfile(projectRoot, '04_RESULTS', material, batchName);
logDir      = fullfile(projectRoot, '05_LOGS', material, batchName);

detDir = fullfile(resultsRoot, 'detections', heightLabel);
qaDir  = fullfile(resultsRoot, 'qa', heightLabel);

% Create folders if missing
if ~exist(detDir, 'dir'), mkdir(detDir); end
if ~exist(qaDir,  'dir'), mkdir(qaDir);  end
if ~exist(logDir, 'dir'), mkdir(logDir); end

fprintf('Input dir: %s\n', inDir);
fprintf('Results detections dir: %s\n', detDir);

%% 3) LOAD VIDEO + METADATA

videoPath = fullfile(inDir, videoFile);
if ~exist(videoPath, 'file')
    error('Video not found: %s', videoPath);
end

v = VideoReader(videoPath);

% Export metadata
fps_export = v.FrameRate;
vidDur_s   = v.Duration;
nFrames    = floor(fps_export * vidDur_s);

fprintf('Video: %s\n', videoFile);
fprintf('Export fps = %.6f, Duration = %.3f s, Approx frames = %d\n', ...
    fps_export, vidDur_s, nFrames);

% Frame indexing convention
frame0 = (0:nFrames-1)';   % 0-based frame counter for exported clip
iMat   = frame0 + 1;       % MATLAB indexing

%% 4) Call functions to detect and save frames

% --- Detection parameter struct ---
params = struct();                  params.alphaG = alphaG;
params.betaB = betaB;               params.doCLAHE = doCLAHE;
params.medianK = medianK;           params.radiusRange = radiusRange;
params.polarity = polarity;         params.sensitivity = sensitivity;
params.edgeThresh = edgeThresh;     params.showPreviewEveryN = showPreviewEveryN;
params.showCenters = showCenters;   params.heightLabel = heightLabel;

% --- Run detection ---
detectOut = detect_circles_per_frame(v, nFrames, params);

% --- Save info struct ---
saveInfo = struct();                saveInfo.material   = material;
saveInfo.material   = material;     saveInfo.batchName  = batchName;
saveInfo.heightLabel = heightLabel; saveInfo.trialNum   = trialNum;
saveInfo.videoFile  = videoFile;    saveInfo.videoPath  = videoPath;
saveInfo.fps_export = fps_export;   saveInfo.detDir     = detDir;
det = save_detections(detectOut, saveInfo);

%% 5) SIMPLE LOGGING (log only)

userNote = input('Optional note for this trial (press Enter to skip): ', 's');

logFile = fullfile(logDir, ['runlog_' datestr(now,'yyyymmdd') '.txt']);
ts = datestr(now,'yyyy-mm-dd HH:MM:SS');

msg = sprintf(['%s | %s %s %s T%02d | fps_export=%.3f | ' ...
               'nFrames=%d | note=%s'], ...
               ts, material, batchName, heightLabel, trialNum, ...
               fps_export, nFramesReadOK, userNote);

fid = fopen(logFile, 'a');
fprintf('Logged to: %s\n', logFile);