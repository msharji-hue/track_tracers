%% select_and_detect.m
clear; close all; clc;
%% 1) USER INPUTS  
videoFile   = 'Test2_H18_edited.mov';
frameStep   = 3;             % show every Nth frame when picking frames
radiusRange = [10 40];       % pixels [10 40]
sensitivity = 0.84;          % 0–1 (higher -> more detections)%0.82 
edgeThresh  = 0.1;           % lower -> more sensitive edges
polarity    = 'bright';      % 'bright' or 'dark'
channel     = 'r';           % 'r','g','b' or anything else -> grayscale

speed_pct   = 2;             % Adobe Premiere speed/duration (%) for THIS video
doCLAHE     = true;          % local contrast boost (minimize distortion during impact)
medianK     = 3;             % 1 disables; small odd ints like 3 or 5
showCenters = true;          % draw blue '+' at centers
ylwWidth    = 1.2;           % circle outline width

%% 1.1) EXPERIMENT PATHS (EDIT PER RUN)

% base directory for this material / packing / height (up to H18)
baseDir = ['/Users/muhannadalsharji/Library/CloudStorage/' ...
    'GoogleDrive-msharji@umich.edu/My Drive/Research /Videos/Edited/New_Drop_Tests/Glass beads/H18'];

% trial number and name
trialNum  = 2;                           % change per run
trialName = sprintf('Trial_%02d', trialNum);

% full path for this specific trial (folder should contain the video)
trialDir = fullfile(baseDir, trialName);

% label for this height (used in detections filename)
heightLabel = 'H18';                      % keep in sync with folder / video

% make sure trial directory exists
if ~exist(trialDir, 'dir')
    mkdir(trialDir);
end

fprintf('Using trial directory: %s\n', trialDir);

%% 2) LOAD VIDEO + BASIC INFO
v = VideoReader(fullfile(trialDir, videoFile));
nFrames = floor(v.FrameRate * v.Duration);
fps     = v.FrameRate;        % camera frame rate (Hz) for t_sec
fprintf('Approx. %d frames in %s\n', nFrames, videoFile);

%% 3) PICK REFERNCE FRAME 

startIdx = [];

for k = 1:frameStep:nFrames
    f = read(v, k);
    figure(1); clf; imshow(f);
    title(sprintf('Candidate start frame %d of %d', k, nFrames));
    drawnow;

    answ = input('Use this as FIRST frame? (y=yes, n=no, q=quit): ','s');
    if strcmpi(answ,'y')
        startIdx = k;
        break
    elseif strcmpi(answ,'q')
        break
    end
end

if isempty(startIdx)
    error('No start frame selected.');
end

% from here on: use every frame from startIdx to end
selectedIdx = startIdx:nFrames;   % step = 1 (all frames)
nSel       = numel(selectedIdx);

%% 4) DETECT CIRCLES PER FRAME (with minimal pre-processing)

results = struct('frame',cell(1,nSel),'centers',[],'radii',[]);

for i = 1:nSel
    frame_idx = selectedIdx(i);
    f = read(v, frame_idx);

    % choose channel
    switch lower(channel)
        case 'r', A = im2double(f(:,:,1));
        case 'g', A = im2double(f(:,:,2));
        case 'b', A = im2double(f(:,:,3));
        otherwise, A = rgb2gray(f);
    end

    % pre-processing
    if doCLAHE, A = adapthisteq(A); end
    if medianK>1, A = medfilt2(A,[medianK medianK]); end

    % circle detection
    [centers, radii] = imfindcircles(A, radiusRange, ...
        'ObjectPolarity',polarity, ...
        'Sensitivity',sensitivity, ...
        'EdgeThreshold',edgeThresh);

    % store
    results(i).frame   = frame_idx;
    results(i).centers = centers;
    results(i).radii   = radii;

    % progress print
    if mod(i,100) == 0 || i == nSel
        fprintf('Processed %d / %d frames\n', i, nSel);
    end
end

%% 5) SAVE DETECTIONS (CSV: one row per circle)

% one detections file per trial, named by height
outFile = fullfile(trialDir, sprintf('detections_%s.csv', heightLabel));

fid = fopen(outFile,'w');
fprintf(fid, 'nFrame,frame_idx,t_ms,t_sec,x,y,r,fps,speed_pct\n');   % header

for i = 1:nSel
    nFrame    = i-1;                 % 0,1,2,...
    frame_idx = results(i).frame;    % actual video frame index
    t_sec     = frame_idx / fps;     % video timeline seconds
    t_ms      = t_sec * (speed_pct/100) * 1000;   % PHYSICAL milliseconds

    C = results(i).centers;          % Nx2 [x y]
    R = results(i).radii;            % Nx1 [r]

    if isempty(C)
        % write a single row to mark the frame (no detections)
        fprintf(fid, '%d,%d,%.3f,%.6f,NaN,NaN,NaN,%.3f,%.2f\n', ...
            nFrame, frame_idx, t_ms, t_sec, fps, speed_pct);
    else
        for k = 1:size(C,1)
            fprintf(fid, '%d,%d,%.3f,%.6f,%.6f,%.6f,%.6f,%.3f,%.2f\n', ...
                nFrame, frame_idx, t_ms, t_sec, C(k,1), C(k,2), R(k), fps, speed_pct);
        end
    end
end
fclose(fid);
disp(['Saved detections to ', outFile]);

