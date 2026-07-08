%% scripts/make_display_video.m
% Quick qualitative display clip — NO plots, NO kinematics.
%
% Flow:  pick video -> pick first/last frame -> export sharpened frames
%        -> detect circles -> click the bed line -> slowed tracer-circle video.
%
% Slowdown = fps_true / outputFps * slowFactor   (e.g. 2009/30*2 = 134x)

clear; clc;
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';   % <- edit if needed
addpath(fullfile(codeDir, 'src'));

% ── Tunables ───────────────────────────────────────────────────────────────
outputFps  = 30;     % output video frame rate
slowFactor = 2;      % write each frame this many times (extra slowdown)

% ── 1) Pick the video ──────────────────────────────────────────────────────
[vfile, vpath] = uigetfile({'*.avi;*.AVI;*.mp4;*.MP4;*.mov;*.MOV','Video Files'}, ...
                           'Select a video');
if isequal(vfile,0), fprintf('Cancelled.\n'); return; end
videoPath = fullfile(vpath, vfile);
[~, vname] = fileparts(vfile);
v = open_video(videoPath);

% ── 2) Pick first / last frame (browse + confirm) ──────────────────────────
startFrame = pick_start_frame(v);
endFrame   = pick_end_frame(v, startFrame);

% ── 3) Export sharpened frames to a working folder ─────────────────────────
framesDir = fullfile(vpath, 'display_frames', vname);
if exist(framesDir,'dir'), rmdir(framesDir,'s'); end   % start clean
export_frames(v, startFrame, endFrame, framesDir, 'sharpen');

% ── 4) Detect circles (same defaults as the main pipeline) ─────────────────
params = struct();
params.radiusRange       = [9 33];
params.sensitivity       = 0.85;
params.edgeThresh        = 0.10;
params.polarity          = 'bright';
params.alphaG            = 0.50;
params.betaB             = 0.50;
params.doCLAHE           = true;
params.medianK           = 3;
params.showPreviewEveryN = 0;       % no popups during detection
params.showCenters       = true;
params.heightLabel       = vname;
fprintf('Detecting circles...\n');
detectOut = detect_circles_per_frame(framesDir, params);

% ── 5) Enter two bed-line points as pixel coordinates ──────────────────────
%     A reference frame is shown so you can read coords with the data cursor.
ff     = sort({dir(fullfile(framesDir,'*.png')).name});
refImg = imread(fullfile(framesDir, ff{end}));   % last frame: bed clearly visible
bedFig = figure('Name','Bed line reference'); imshow(refImg);
title('Use the data cursor to read two points on the bed line, then type them in');
datacursormode(bedFig, 'on');
a1 = inputdlg('Bed line point 1 as [x, y]:', 'Bed Line', 1, {'[100, 120]'});
a2 = inputdlg('Bed line point 2 as [x, y]:', 'Bed Line', 1, {'[500, 120]'});
if isempty(a1) || isempty(a2), close(bedFig); fprintf('Cancelled.\n'); return; end
close(bedFig);
bedPoint1 = str2num(a1{1});  %#ok<ST2NM>
bedPoint2 = str2num(a2{1});  %#ok<ST2NM>
lineA = bedPoint1(2) - bedPoint2(2);
lineB = bedPoint2(1) - bedPoint1(1);
lineC = bedPoint1(1)*bedPoint2(2) - bedPoint2(1)*bedPoint1(2);

% ── 6) Render params ───────────────────────────────────────────────────────
rp = struct();
rp.outputFps   = outputFps;
rp.slowFactor  = slowFactor;
rp.fps_true    = v.FrameRate;
rp.markerColor = [0.35 0.95 0.25];   % lime green tracer circles (set [1 1 1] for white)
rp.bedColor    = [1.00 0.90 0.20];   % yellow dashed bed line
rp.circLW      = 2;
rp.trailLen    = 8;                  % set 0 to disable the fading trail
rp.trailAlpha  = 0.5;
rp.trailSize   = 6;
rp.rotateCCW   = true;               % rot90(img,1), matching make_trial_video.m
rp.bedPoint1   = bedPoint1;          % original (unrotated) pixel coords
rp.bedPoint2   = bedPoint2;

fprintf('\nfps_true=%.0f | outputFps=%d | slowFactor=%d | Slowdown=%.0fx\n', ...
    v.FrameRate, outputFps, slowFactor, v.FrameRate/outputFps*slowFactor);

% ── 7) Preview, then optionally save ───────────────────────────────────────
fprintf('Rendering preview...\n');
make_tracer_video(framesDir, detectOut, rp, false);

if strcmp(questdlg('Save video?','Confirm','Yes','No','Yes'), 'Yes')
    outDir = uigetdir(vpath, 'Select output folder');
    if ~isequal(outDir,0)
        rp.outputName = fullfile(outDir, [vname '_tracer.mp4']);
        make_tracer_video(framesDir, detectOut, rp, true);
    end
else
    fprintf('Not saved.\n');
end