%% SECTION 1: USER INPUTS + ROD AXIS + y0
clear; close all; clc;

% mark tracers throughout video 

% user inputs
videoFile = 'Test2_H18_edited.mov';   % must live inside trialDir

%% 1.1) EXPERIMENT PATHS (EDIT PER RUN)

% base directory for this material / packing / height (up to H18)
baseDir = ['/Users/muhannadalsharji/Library/CloudStorage/' ...
    'GoogleDrive-msharji@umich.edu/My Drive/Research /Videos/Edited/New_Drop_Tests/Glass beads/H18'];

% trial number and name
trialNum  = 2;                           % change per run
trialName = sprintf('Trial_%02d', trialNum);

% full path for this specific trial (folder should contain the video +
% detections from select_and_detect)
trialDir = fullfile(baseDir, trialName);

% label for this height (used in detections filename)
heightLabel = 'H18';                      % keep in sync with folder / video

% make sure trial directory exists
if ~exist(trialDir, 'dir')
    mkdir(trialDir);
end

fprintf('Using trial directory: %s\n', trialDir);

% rod + marker geometry (mm)
rod_length_mm          = 28.15;
tip_to_first_center_mm = 0.40;
marker_diameter_mm     = 2.00;
marker_pitch_mm        = 4.00;

% load detections and get reference frame (nFrame = 0)
detFile = fullfile(trialDir, sprintf('detections_%s.csv', heightLabel));
det     = readtable(detFile);

% drop rows that represent frames with no circles (x = NaN)
det = det(~isnan(det.x), :);
ref_rows      = det.nFrame == 0;
ref_frame_idx = det.frame_idx(ref_rows);
ref_frame_idx = ref_frame_idx(1);  % assume consistent

% show reference frame
v       = VideoReader(fullfile(trialDir, videoFile));
ref_img = read(v, ref_frame_idx);

figure(1); clf; imshow(ref_img); hold on;
title(sprintf('Reference frame %d', ref_frame_idx));

% pick rod axis: TOP then TIP
disp('Click rod TOP, then rod TIP');
[xr, yr] = ginput(2);
rod_top = [xr(1) yr(1)];
rod_tip = [xr(2) yr(2)];

plot(rod_top(1), rod_top(2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
plot(rod_tip(1), rod_tip(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
plot([rod_top(1) rod_tip(1)], [rod_top(2) rod_tip(2)], 'b-', 'LineWidth', 2);

% pick y0 line (bed): two clicks along it
disp('Click two points along y0 line (bed surface)');
[xy, yy] = ginput(2);
plot([xy(1) xy(2)], [yy(1) yy(2)], 'y-', 'LineWidth', 2);

y0 = mean(yy);   % single y0 level in pixels

disp('Section 1 done.');

%% SECTION 2: MARKER IDS FROM REFERENCE FRAME

% user input: pixels to mm (same in x and y)
mm_per_px = 0.0869;  % double check this

% make sure detections are loaded
if ~exist('det','var')
    detFile = fullfile(trialDir, sprintf('detections_%s.csv', heightLabel));
    det     = readtable(detFile);
end

% use reference frame nFrame = 0
ref_rows = det.nFrame == 0;
x0  = det.x(ref_rows);
y0c = det.y(ref_rows);

% remove any NaNs
valid = ~isnan(x0) & ~isnan(y0c);
x0  = x0(valid);
y0c = y0c(valid);

if isempty(x0)
    error('No circles found in reference frame (nFrame = 0).');
end

% sort markers along rod axis (from TIP to TOP)
u = (rod_top - rod_tip);
u = u / norm(u);  % unit vector along rod

proj = ( [x0 y0c] - rod_tip ) * u.';   % projection along rod
[proj_sorted, idx] = sort(proj);       % small = near tip

x0  = x0(idx);
y0c = y0c(idx);

nMarkers = numel(x0);
markerID = (1:nMarkers)';   % 1 = closest to tip

% pixel and mm coordinates in reference frame
markers.ID     = markerID;
markers.x_px0  = x0;
markers.y_px0  = y0c;
markers.x_mm0  = x0 * mm_per_px;
markers.y_mm0  = y0c * mm_per_px;

% assign colors for plotting (same color per ID later)
markers.colors = lines(nMarkers);

% quick plot check
figure(1); hold on;
for k = 1:nMarkers
    c = markers.colors(k,:);
    plot(x0(k), y0c(k), 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
         'Color', c, 'MarkerFaceColor', 'none');
    text(x0(k)+5, y0c(k), sprintf('%d', markerID(k)), ...
         'Color', c, 'FontSize', 9);
end

disp(['Section 2 done: assigned ', num2str(nMarkers), ' marker IDs.']);

%% SECTION 3: TRACK MARKERS + DEPTH VS TIME

% unit vector along rod (tip -> top)
u = (rod_top - rod_tip);
u = u / norm(u);

nMarkers = numel(markers.ID);
frames   = unique(det.nFrame);
nF       = numel(frames);

% map analysis frames (nFrame) to actual video frame indices
frame_idx_vec = nan(1, nF);
for fi = 1:nF
    rows = det.nFrame == frames(fi);
    frame_idx_vec(fi) = det.frame_idx(find(rows, 1, 'first'));
end

x_px     = nan(nMarkers, nF);
y_px     = nan(nMarkers, nF);
time_ms  = nan(1, nF);
depth_mm = nan(1, nF);

for fi = 1:nF
    nFrame = frames(fi);
    rows   = det.nFrame == nFrame;

    xx = det.x(rows);
    yy = det.y(rows);
    tt = det.t_ms(rows);

    if ~isempty(tt)
        time_ms(fi) = tt(1);   % same time for all circles in this frame
    end

    % drop NaNs
    valid = ~isnan(xx) & ~isnan(yy);
    xx = xx(valid);
    yy = yy(valid);

    Nd = numel(xx);
    if Nd == 0
        continue
    end

    % project detections along rod axis and sort (tip -> top)
    proj_det = ([xx yy] - rod_tip) * u.';
    [~, idxDet] = sort(proj_det);   % small = near tip
    xx = xx(idxDet);
    yy = yy(idxDet);

    % number of detections AFTER sorting
    Nd = numel(xx);

    % assign detections to TOP markers only
    if Nd > nMarkers
        % keep only the top-most nMarkers detections (closest to rod_top)
        xx = xx(end-nMarkers+1:end);
        yy = yy(end-nMarkers+1:end);
        Nd = nMarkers;
    end

    top_idx = (nMarkers-Nd+1):nMarkers;

    x_px(top_idx, fi) = xx;
    y_px(top_idx, fi) = yy;
end

for fi = 1:nF
    visible = ~isnan(x_px(:,fi)) & ~isnan(y_px(:,fi));
    if ~any(visible)
        depth_mm(fi) = NaN;
        continue
    end

    pos  = [x_px(visible,fi) y_px(visible,fi)];
    pos0 = [markers.x_px0(visible) markers.y_px0(visible)];
    dvec = pos - pos0;
    ds   = mean(dvec * u.');

    tip_t = rod_tip + ds * u;
    top_t = rod_top + ds * u;

    y_tip = tip_t(2);
    y_top = top_t(2);

    if y_tip <= y0 && y_top <= y0
        L_above = rod_length_mm;       % rod fully above bed
    elseif y_tip >= y0 && y_top >= y0
        L_above = 0;                   % rod fully below bed
    else
        s_int   = (y0 - y_tip) / (y_top - y_tip);
        L_above = (1 - s_int) * rod_length_mm;
    end

    depth_mm(fi) = rod_length_mm - L_above;
end

x_mm = x_px * mm_per_px;
y_mm = y_px * mm_per_px;

disp('Section 3 done: marker tracks and depth vs time computed.');


%% SECTION 4: VELOCITY VS TIME (mm/ms)

v_mm_per_ms = nan(size(depth_mm));

for fi = 2:nF
    d1 = depth_mm(fi-1);
    d2 = depth_mm(fi);
    t1 = time_ms(fi-1);
    t2 = time_ms(fi);

    if isnan(d1) || isnan(d2) || isnan(t1) || isnan(t2)
        continue
    end

    dt = t2 - t1;
    if dt <= 0
        continue
    end

    v_mm_per_ms(fi) = (d2 - d1) / dt;   % mm/ms
end

disp('Section 4 done: velocity computed and saved.');


%% SECTION 5: ACCELERATION

% acceleration a(t) from v(t), central difference (mm/ms^2)
a_mm_per_ms2 = nan(size(v_mm_per_ms));

for i = 2:(nF-1)
    v1 = v_mm_per_ms(i-1);
    v2 = v_mm_per_ms(i+1);
    t1 = time_ms(i-1);
    t2 = time_ms(i+1);

    if any(isnan([v1 v2 t1 t2]))
        continue
    end

    dt = t2 - t1;
    if dt <= 0
        continue
    end

    a_mm_per_ms2(i) = (v2 - v1) / dt;   % mm/ms^2
end

% add gravity: g = 9.81 m/s^2 = 0.00981 mm/ms^2
g_mm_per_ms2 = 9.81e-3;
a_plus_g_mm_per_ms2 = a_mm_per_ms2 + g_mm_per_ms2;

disp('Section 5 done: acceleration computed and saved.');

%% SECTION 6: RENDER TRACER OVERLAY VIDEO

% save tracer overlay video directly in this trial folder
outVideoFile = fullfile(trialDir, sprintf('%s_tracers.mp4', trialName));
disp(['Saving tracer video to: ', outVideoFile]);

% set up video reader and writer
v  = VideoReader(fullfile(trialDir, videoFile));
vw = VideoWriter(outVideoFile, 'MPEG-4');
vw.FrameRate = v.FrameRate;   % match original edited video fps
open(vw);

% figure for rendering
hFig = figure(67);
set(hFig, 'Color','w');   % you can add 'Visible','off' here if you like

for fi = 1:nF
    frame_idx = frame_idx_vec(fi);
    frame     = read(v, frame_idx);

    figure(hFig); clf;
    imshow(frame); hold on;

    % draw bed line y0 (optional)
    plot([1 size(frame,2)], [y0 y0], 'y--', 'LineWidth', 1.2);

    % overlay tracers: hollow colored circle, center mark, and marker ID
    for k = 1:nMarkers
        xk = x_px(k,fi);
        yk = y_px(k,fi);

        if isnan(xk) || isnan(yk)
            continue
        end

        c  = markers.colors(k,:);
        id = markers.ID(k);

        % hollow colored circle
        plot(xk, yk, 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
             'Color', c, 'MarkerFaceColor', 'none');

        % marked center (black '+')
        plot(xk, yk, '+', 'MarkerSize', 6, 'LineWidth', 1.5, ...
             'Color', 'k');

        % marker ID label (slightly offset)
        text(xk + 5, yk - 5, sprintf('%d', id), ...
             'Color', c, 'FontSize', 9, 'FontWeight','bold');
    end

    F = getframe(hFig);      % capture current figure as a frame
    writeVideo(vw, F);       % append it to the video
end

close(vw);

disp(['Section 6 done: tracer overlay video saved to ', outVideoFile]);

%% SECTION 7: SAVE RESULTS

% save tracking/kinematics for this trial in the same folder
trackFile = fullfile(trialDir, 'track_results.mat');

save(trackFile, ...
    'nMarkers','markers', ...
    'x_px','y_px','x_mm','y_mm', ...
    'time_ms','depth_mm', ...
    'v_mm_per_ms', ...
    'a_mm_per_ms2','a_plus_g_mm_per_ms2', ...
    'rod_top','rod_tip','y0','mm_per_px','videoFile', ...
    'frame_idx_vec');

disp(['Section 7 done: saved results to ', trackFile]);




