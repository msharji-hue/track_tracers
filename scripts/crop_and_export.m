% crop_and_export.m
% Picks start and end frames visually, saves all frames in between as PNGs.

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

% ── USER INPUTS ───────────────────────────────────────────────
inputDir   = uigetdir(pwd, 'Select folder containing your video');
[file, ~]  = uigetfile({'*.mp4;*.avi;*.mov', 'Video Files'}, 'Select a video', inputDir);
videoPath  = fullfile(inputDir, file);

parentDir  = uigetdir(pwd, 'Select where to save output folder');
folderName = inputdlg('Name your output folder:', 'Output Folder', 1, {'frames_output'});
outputDir  = fullfile(parentDir, folderName{1});
filterType = 'sharpen';   % 'none' | 'grayscale' | 'gaussian' | 'sharpen'
% ─────────────────────────────────────────────────────────────

v          = open_video(videoPath);
startFrame = pick_start_frame(v);
endFrame   = pick_end_frame(v, startFrame);
export_frames(v, startFrame, endFrame, outputDir, filterType);
 