%% track_tracers_2.m
clear; close all; clc;

%% 1) SETUP
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

[file, detDir] = uigetfile('*.mat', 'Select detections .mat file');
detFile        = fullfile(detDir, file);

material = inputdlg({'Material','Batch','Height label','Trial number', ...
                     'Freefall height (cm)','SG order','SG window'}, ...
                     'Trial Info', 1, ...
                     {'GB','Batch1','H06','1','15.24','4','91'});
trialInfo = struct('material',    material{1}, ...
                   'batchName',   material{2}, ...
                   'heightLabel', material{3}, ...
                   'trialNum',    str2double(material{4}), ...
                   'h_cm',        str2double(material{5}), ...
                   'sgOrder',     str2double(material{6}), ...
                   'sgWindow',    str2double(material{7}), ...
                   'fps_true',    2350);

outRoot   = uigetdir(pwd, 'Select results root folder');
subFolders = {'detections','tracks','kinematics','figures','qa','logs'};
for k = 1:numel(subFolders)
    mkdir(fullfile(outRoot, subFolders{k}, trialInfo.heightLabel));
end

S          = load(detFile);
detections = S.det.detect.centersCell;
if isempty(detections{1}), error('No detections in first frame.'); end
x0 = detections{1}(:,1);
y0 = detections{1}(:,2);

%% 2) BED LINE
bedPoint1 = [17, 0];  bedPoint2 = [17, 47];
lineA = bedPoint1(2) - bedPoint2(2);
lineB = bedPoint2(1) - bedPoint1(1);
lineC = bedPoint1(1)*bedPoint2(2) - bedPoint2(1)*bedPoint1(2);

%% 3) MARKER IDS
firstFrameCenters = detections{1};
if isempty(firstFrameCenters), error('No circles in first frame.'); end
numberOfMarkers = size(firstFrameCenters, 1);
initialMarkers  = struct('id',     (1:numberOfMarkers)', ...
                         'x',      x0, ...
                         'y',      y0, ...
                         'colors', lines(numberOfMarkers));

%% 4) TRACK + KINEMATICS
dt      = 1 / trialInfo.fps_true;
mmPerPx = 0.1079;
g_cm_s2 = 980;

[trackedX, trackedY] = track_markers(detections, firstFrameCenters, 100);
[t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, impact_index, toeMarkerID, toePx] = toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, dt, mmPerPx, -373.87);

v0_cm_s   = v_smooth(impact_index);
stopFrame = find(v_smooth(impact_index:end) <= 0, 1, 'first') + impact_index - 1;
d_final_cm = depthRod_cm(stopFrame);
t_stop_s   = t_s(stopFrame);
fprintf('v0 = %.2f cm/s  |  d_final = %.3f cm  |  t_stop = %.4f s\n', v0_cm_s, d_final_cm, t_stop_s);

%% 5) PLOT
a_plus_g = -a_smooth - g_cm_s2;
colors   = lines(size(trackedX, 1));

% --- Figure 1: Kinematics (z, v, a+g) ---
figure('Name', 'Kinematics');
subplot(3,1,1)
plot(t_s(1:stopFrame), depthRod_cm(1:stopFrame), 'k.', 'MarkerSize', 4); hold on
plot(t_s, z_smooth, 'r', 'LineWidth', 1.5)
xline(0, '--k', 'impact'); ylabel('depth (cm)'); legend('raw','poly fit'); grid on
title(sprintf('%s T%02d — Kinematics', trialInfo.heightLabel, trialInfo.trialNum));

subplot(3,1,2)
plot(t_s, v_smooth, 'b', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k'); ylabel('v (cm/s)'); grid on

subplot(3,1,3)
plot(t_s, a_plus_g, 'm', 'LineWidth', 1.5)
xline(0, '--k'); yline(0, '--k');
ylabel('a+g (cm/s²)'); xlabel('t (s)'); grid on

%% --- Figure 2: Marker tracks + detections ---
figure('Name', 'Tracking QA');
subplot(2,1,1)
for m = 1:size(trackedX,1)
    plot(trackedX(m,:), '.-', 'Color', colors(m,:)); hold on
end
xline(impact_index, 'r', 'impact'); xlabel('frame'); ylabel('x (px)');
title('all marker tracks'); grid on; legend(arrayfun(@(i) sprintf('M%d',i), 1:size(trackedX,1), 'UniformOutput', false), 'Location', 'eastoutside');

subplot(2,1,2)
nDetected = cellfun(@(c) size(c,1), detections);
bar(nDetected, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none');
xline(impact_index, 'r', 'impact'); xlabel('frame'); ylabel('n detections'); grid on
title('detections per frame');

% --- Figure 3: Phase space ---
figure('Name', 'Phase Space');

% use normalized time 0-1 for better color variation
t_norm = (t_s - t_s(impact_index)) / (t_s(stopFrame) - t_s(impact_index));
t_norm = max(0, min(1, t_norm));   % clip to 0-1

subplot(3,1,1)
scatter(z_smooth, v_smooth, 20, t_norm, 'filled');
colorbar; colormap(jet);
xlabel('depth (cm)'); ylabel('v (cm/s)');
title('phase space: v vs depth (color = normalized time)'); grid on

subplot(3,1,2)
scatter(z_smooth, a_plus_g, 20, t_norm, 'filled');
colorbar; colormap(jet);
xlabel('depth (cm)'); ylabel('a+g (cm/s²)');
title('a+g vs depth (color = normalized time)'); grid on

subplot(3,1,3)
scatter(v_smooth.^2, a_plus_g, 20, z_smooth, 'filled');
colorbar; colormap(parula);
xlabel('v² (cm²/s²)'); ylabel('a+g (cm/s²)');
title('Katsuragi Fig 3a: a+g vs v² (color = depth)'); grid on

%% 6) SAVE
a_plus_g   = -a_smooth - g_cm_s2;
kinematics = struct('t_s',          t_s, ...
                    'depthRod_cm',  depthRod_cm, ...
                    'z_smooth',     z_smooth, ...
                    'v_smooth',     v_smooth, ...
                    'a_smooth',     a_smooth, ...
                    'a_plus_g',     a_plus_g, ...
                    'trackedX',     trackedX, ...
                    'trackedY',     trackedY, ...
                    'impact_index', impact_index, ...
                    'stopFrame',    stopFrame, ...
                    'toeMarkerID',  toeMarkerID, ...
                    'toePx',        toePx);

H_cm    = trialInfo.h_cm + d_final_cm;
scalars = struct('v0_cm_s',    v0_cm_s, ...
                 'd_final_cm', d_final_cm, ...
                 't_stop_s',   t_stop_s, ...
                 'h_cm',       trialInfo.h_cm, ...
                 'H_cm',       H_cm);

fprintf('v0 = %.2f cm/s  |  d = %.3f cm  |  H = %.3f cm  |  t_stop = %.4f s\n', ...
    v0_cm_s, d_final_cm, H_cm, t_stop_s);

save_kinematics(trialInfo, kinematics, scalars, outRoot);

figDir = fullfile(outRoot, 'figures', trialInfo.heightLabel);
for figNum = 1:3
    fh      = figure(figNum);
    figBase = sprintf('%s_T%02d_fig%d', trialInfo.heightLabel, trialInfo.trialNum, figNum);
    exportgraphics(fh, fullfile(figDir, [figBase '.png']), 'Resolution', 300);
    savefig(fh, fullfile(figDir, [figBase '.fig']));
    fprintf('Figure %d saved: %s\n', figNum, figBase);
end