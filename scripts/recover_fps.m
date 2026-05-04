%% scripts/recover_fps.m
% Backfills video_meta.mat into existing frames folders.
% Run once per trial before select_and_detect if crop_and_export
% was run before fps auto-save was added.
%
% Workflow:
%   1) Select the original video file
%   2) Select the corresponding frames folder
%   3) Repeat for each trial — cancel dialog when done
%
% After running this, select_and_detect will automatically
% read fps from video_meta.mat with no manual entry needed.

clear; clc;

codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

fprintf('=== recover_fps ===\n');
fprintf('Pairs each video to its frames folder and saves fps.\n');
fprintf('Cancel the video dialog when done.\n\n');

count = 0;

while true
    % ── Select original video ─────────────────────────────────────────────
    [file, path] = uigetfile({'*.mp4;*.avi;*.mov;*.MP4', 'Video Files'}, ...
        sprintf('Select original video (trial %d) — cancel when done', count+1));
    if isequal(file, 0)
        fprintf('\nCancelled — %d trial(s) processed.\n', count);
        break;
    end

    % ── Read FPS from video header ────────────────────────────────────────
    v        = VideoReader(fullfile(path, file));
    fps_true = v.FrameRate;
    fprintf('[%d]  %s  →  FPS = %.4f\n', count+1, file, fps_true);

    % ── Select corresponding frames folder ────────────────────────────────
    framesDir = uigetdir(path, ...
        sprintf('Select frames folder for:  %s', file));
    if isequal(framesDir, 0)
        fprintf('    Skipped (no folder selected).\n\n');
        continue;
    end

    % ── Save metadata ─────────────────────────────────────────────────────
    save(fullfile(framesDir, 'video_meta.mat'), 'fps_true');
    fprintf('    Saved video_meta.mat  →  %s\n\n', framesDir);
    count = count + 1;
end

fprintf('Done. Run select_and_detect normally — fps will load automatically.\n');
