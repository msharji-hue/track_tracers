%% SECTION 1: Track Tracers 
clear; close all; clc;

%% 1) SETUP
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

% --- Folder + trial info via UI ---
[file, detDir] = uigetfile('*.mat', 'Select detections .mat file');
detFile    = fullfile(detDir, file);

material   = inputdlg({'Material','Batch','Height label','Trial number'}, ...
                       'Trial Info', 1, {'GB','Batch1','H06','1'});
trialInfo  = struct('material',    material{1}, ...
                    'batchName',   material{2}, ...
                    'heightLabel', material{3}, ...
                    'trialNum',    str2double(material{4}));

trialID    = sprintf('T%02d', trialInfo.trialNum);
baseName   = sprintf('%s_%s_long', trialInfo.heightLabel, trialID);
outRoot    = uigetdir(pwd, 'Select results root folder');

% --- Create subfolders ---
subFolders = {'detections','tracks','kinematics','figures','qa','logs'};
for k = 1:numel(subFolders)
    folder = fullfile(outRoot, subFolders{k}, trialInfo.heightLabel);
    if ~exist(folder, 'dir'), mkdir(folder); end
end

% --- Load detections ---
S          = load(detFile);
detections = S.det.detect.centersCell;
if isempty(detections{1}), error('No detections in first frame.'); end
x0 = detections{1}(:,1);
y0 = detections{1}(:,2);

%% 2) BED LINE FROM TWO POINTS
bedPoint1 = [17, 0];   bedPoint2 = [17, 47];

x1 = bedPoint1(1);  y1 = bedPoint1(2);
x2 = bedPoint2(1);  y2 = bedPoint2(2);

% Implicit line form: A*x + B*y + C = 0
lineA = y1 - y2;   % = -47
lineB = x2 - x1;   % =   0  (vertical line)
lineC = x1*y2 - x2*y1;  % = 17*47 - 17*0 = 799

% Depth = perpendicular distance from rod center to bed line
% For a vertical line this simplifies to: depth = cx - x_line
% Positive = rod is to the RIGHT of the line (into the bed)

%% 3) ASSIGN MARKER IDS FROM FIRST FRAME
firstFrameCenters = detections{1};
if isempty(firstFrameCenters), error('No circles detected in first frame.'); end

numberOfMarkers  = size(firstFrameCenters, 1);
initialMarkers   = struct('id',     (1:numberOfMarkers)', ...
                          'x',      x0, ...
                          'y',      y0, ...
                          'colors', lines(numberOfMarkers));

%% 4) TRACK + KINEMATICS
dt      = 1 / trialInfo.fps_true;
mmPerPx = 0.1079;
g_cm_s2 = 980;

[trackedX, trackedY] = track_markers(detections, firstFrameCenters);
[t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, impact_index, toeMarkerID, toePx] = toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, dt, mmPerPx, -373.87);
a_plus_g = g_cm_s2 - a_smooth;
a_smooth_plot = -a_smooth;            % flip sign for plotting
a_smooth_plot(1:impact_index-1) = nan;
v0_cm_s  = v_smooth(1);
fprintf('Impact frame: %d  |  v0 = %.2f cm/s\n', impact_index, v0_cm_s);

stopFrame = find(v_smooth <= 0, 1, 'first');
if isempty(stopFrame), stopFrame = numel(t_s); end
%% PLOT — sanity check
figure; set(gcf, 'Position', [100 100 1200 900]);

subplot(2,2,1)
plot(t_s(1:stopFrame), depthRod_cm(1:stopFrame), 'k.', 'MarkerSize', 4); hold on
plot(t_s(1:stopFrame), z_smooth(1:stopFrame), 'r', 'LineWidth', 1.5)
xline(0, '--k', 'impact'); ylabel('depth (cm)'); legend('raw','poly fit');
grid on; axis tight

subplot(2,2,2)
plot(t_s(1:stopFrame), v_smooth(1:stopFrame), 'b', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k'); ylabel('v (cm/s)');
grid on; axis tight

subplot(2,2,3)
plot(t_s(1:stopFrame), a_plus_g(1:stopFrame), 'm', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k'); ylabel('a+g (cm/s²)'); xlabel('t (s)');
grid on; axis tight

subplot(2,2,4)
plot(t_s(1:stopFrame), a_smooth_plot(1:stopFrame), 'g', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k'); ylabel('a (cm/s²)'); xlabel('t (s)');
grid on; axis tight

sgtitle(sprintf('%s %s %s T%02d', trialInfo.material, trialInfo.batchName, trialInfo.heightLabel, trialInfo.trialNum));