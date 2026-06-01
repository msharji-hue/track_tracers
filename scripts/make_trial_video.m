%% scripts/make_trial_video.m
% Annotated video: rotated portrait video (left) + z(t) and v(t) plots (right).
%
% Slowdown = fps_true / outputFps * slowFactor
%   e.g. 2009/30*3 = 201x

clear; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

% ── Select inputs ─────────────────────────────────────────────────────────
framesDir          = uigetdir(pwd, 'Select frames folder');
[detFile, detPath] = uigetfile('*.mat', 'Select detections .mat');
[kinFile, kinPath] = uigetfile('*.mat', 'Select kinematics .mat');

detS = load(fullfile(detPath, detFile));
kinS = load(fullfile(kinPath, kinFile));
det  = detS.det;
kin  = kinS.kinematics;

% ── Bed line ──────────────────────────────────────────────────────────────
frameFiles = sort({dir(fullfile(framesDir,'*.png')).name});
firstFrame = imread(fullfile(framesDir, frameFiles{kin.impact_index}));
figure; imshow(firstFrame);
title('Use data cursor to find two points on the bed line, then close');

ans1 = inputdlg('Enter bed line point 1 as [x, y]:','Bed Line',1,{'[100, 120]'});
ans2 = inputdlg('Enter bed line point 2 as [x, y]:','Bed Line',1,{'[500, 120]'});
close;

bedPoint1 = str2num(ans1{1});  %#ok<ST2NM>
bedPoint2 = str2num(ans2{1});  %#ok<ST2NM>
lineA = bedPoint1(2) - bedPoint2(2);
lineB = bedPoint2(1) - bedPoint1(1);
lineC = bedPoint1(1)*bedPoint2(2) - bedPoint2(1)*bedPoint1(2);

% ── Timing ────────────────────────────────────────────────────────────────
outputFps  = 30;
slowFactor = 2;      % each frame written this many times
preFrames  = round(det.meta.fps_export * 0.008);   % 8ms before impact
postFrames = round(det.meta.fps_export * 0.010);   % 10ms after stop

fprintf('fps_true=%.0f | outputFps=%d | slowFactor=%d | Slowdown=%.0fx\n',...
    det.meta.fps_export, outputFps, slowFactor,...
    det.meta.fps_export/outputFps*slowFactor);

% ── Parameters ────────────────────────────────────────────────────────────
params.fps_true    = det.meta.fps_export;
params.outputFps   = outputFps;
params.slowFactor  = slowFactor;
params.mmPerPx     = 0.1079;           % ← update to your calibration
params.trailLen    = 8;
params.trailAlpha  = 0.6;
params.trailSize   = 6;                % trail dot size (pixels)
params.markerColor = [0.20 0.80 1.00]; % cyan
params.bedColor    = [1.00 0.90 0.20]; % yellow
params.circLW      = 3;                % tracer circle line width
params.bedPoint1   = bedPoint1;
params.bedPoint2   = bedPoint2;
params.lineA       = lineA;
params.lineB       = lineB;
params.lineC       = lineC;
params.startFrame  = max(1, kin.impact_index - preFrames);
params.stopFrame   = min(numel(frameFiles), kin.stopFrame + postFrames);

% ── Preview ───────────────────────────────────────────────────────────────
fprintf('\nRendering preview...\n');
make_annotated_video(framesDir, det, kin, params, false);

% ── Save ──────────────────────────────────────────────────────────────────
answer = questdlg('Save video?','Confirm','Yes','No','Yes');
if strcmp(answer,'Yes')
    outDir            = uigetdir(pwd,'Select output folder');
    params.outputName = fullfile(outDir,'trial_annotated.mp4');
    make_annotated_video(framesDir, det, kin, params, true);
else
    fprintf('Not saved.\n');
end