%% track_tracers_2.m
clear; close all; clc;

%% 1) SETUP
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

[file, detDir] = uigetfile('*.mat', 'Select detections .mat file');
detFile        = fullfile(detDir, file);

material       = inputdlg({'Material','Batch','Height label','Trial number'}, ...
                           'Trial Info', 1, {'GB','Batch1','H06','1'});
trialInfo      = struct('material',    material{1}, ...
                        'batchName',   material{2}, ...
                        'heightLabel', material{3}, ...
                        'trialNum',    str2double(material{4}), ...
                        'fps_true',    2350);

trialID        = sprintf('T%02d', trialInfo.trialNum);
baseName       = sprintf('%s_%s_long', trialInfo.heightLabel, trialID);
outRoot        = uigetdir(pwd, 'Select results root folder');

subFolders = {'detections','tracks','kinematics','figures','qa','logs'};
for k = 1:numel(subFolders)
    mkdir(fullfile(outRoot, subFolders{k}, trialInfo.heightLabel));
end

S          = load(detFile);
detections = S.det.detect.centersCell;
if isempty(detections{1}), error('No detections in first frame.'); end
x0 = detections{1}(:,1);
y0 = detections{1}(:,2);

%% 2) BED LINE FROM TWO POINTS
bedPoint1 = [17, 0];   bedPoint2 = [17, 47];
x1 = bedPoint1(1);  y1 = bedPoint1(2);
x2 = bedPoint2(1);  y2 = bedPoint2(2);
lineA = y1 - y2;
lineB = x2 - x1;
lineC = x1*y2 - x2*y1;

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
v0_cm_s  = v_smooth(impact_index);
fprintf('Impact frame: %d  |  v0 = %.2f cm/s\n', impact_index, v0_cm_s);

stopFrame = find(v_smooth(impact_index:end) <= 0, 1, 'first') + impact_index - 1;
fprintf('Stop frame: %d  |  t_stop = %.4f s\n', stopFrame, t_s(stopFrame));

%% 5) PLOT — sanity check
figure;
subplot(3,1,1)
plot(t_s, depthRod_cm, 'k.', 'MarkerSize', 4); hold on
plot(t_s, z_smooth, 'r', 'LineWidth', 1.5)
xline(0, '--k', 'impact'); ylabel('depth (cm)'); legend('raw','poly fit'); grid on

subplot(3,1,2)
plot(t_s, v_smooth, 'b', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k'); ylabel('v (cm/s)'); grid on

subplot(3,1,3)
plot(t_s, a_plus_g, 'm', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k'); ylabel('a+g (cm/s²)'); xlabel('t (s)'); grid on
