%% track_tracers_2.m
clear; close all; clc;

%% 1) SETUP
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

[file, detDir] = uigetfile('*.mat', 'Select detections .mat file');
detFile        = fullfile(detDir, file);
S_check        = load(detFile);
framesDir      = S_check.det.meta.framesDir;

metaFile = fullfile(framesDir, 'video_meta.mat');
if exist(metaFile, 'file')
    meta     = load(metaFile);
    fps_true = meta.fps_true;
    fprintf('Loaded FPS: %.4f\n', fps_true);
else
    fps_true = 2250;
    warning('No video_meta.mat — using default fps=%.0f', fps_true);
end

material  = inputdlg({'Material','Batch','Height label','Trial number', ...
                      'Freefall height (cm)'}, 'Trial Info', 1, ...
                      {'GB','Batch 4','LDH_0','1','1'});
trialInfo = struct('material',    material{1}, ...
                   'batchName',   material{2}, ...
                   'heightLabel', material{3}, ...
                   'trialNum',    str2double(material{4}), ...
                   'h_cm',        str2double(material{5}), ...
                   'fps_true',    fps_true);

outRoot = uigetdir(pwd, 'Select results root folder');
for k = {'detections','tracks','kinematics','figures','qa','logs'}
    mkdir(fullfile(outRoot, k{1}, trialInfo.heightLabel));
end

%% 2) PHYSICAL COORDINATE SYSTEM
mmPerPx      = 0.1079;
impactDistPx = -400;
bedPoint1    = [23, 0];
bedPoint2    = [23, 32];
lineA        = bedPoint1(2) - bedPoint2(2);
lineB        = bedPoint2(1) - bedPoint1(1);
lineC        = bedPoint1(1)*bedPoint2(2) - bedPoint2(1)*bedPoint1(2);

%% 3) TRACK
S          = load(detFile);
detections = S.det.detect.centersCell;
if isempty(detections{1}), error('No detections in first frame.'); end

firstFrameCenters    = detections{1};
[trackedX, trackedY] = track_markers(detections, firstFrameCenters, 100);
nMarkers             = size(trackedX, 1);
cmap_markers         = lines(nMarkers);

%% 4) KINEMATICS
dt            = 1 / fps_true;
g_cm_s2       = 980;
postCapFrames = round(fps_true * 0.010);   % 10ms after t_stop for video

[t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, ...
 impact_index, refMarkerID, rodBedDist_px, stopFrame, sgOrder, sgWindow] = ...
    toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, ...
                   dt, mmPerPx, impactDistPx, ...
                   'postCapFrames', postCapFrames);

a_plus_g   = - a_smooth - g_cm_s2;
v0_cm_s    = v_smooth(impact_index);
d_final_cm = z_smooth(stopFrame);
t_stop_s   = t_s(stopFrame);

fprintf('v0       = %.2f cm/s\n',  v0_cm_s);
fprintf('z_peak   = %.3f cm\n',    d_final_cm);
fprintf('t_stop   = %.4f s\n',     t_stop_s);

%% 5) FIGURES
% ── Shared x limits ───────────────────────────────────────────────────────
tPre   = max(t_s(1), -0.005);
tPost  = t_s(stopFrame);
tRange = [tPre, tPost];

% Fig 1 — Kinematics
figure('Name','Kinematics');

subplot(3,1,1)
plot(t_s(1:stopFrame), depthRod_cm(1:stopFrame), 'k.', 'MarkerSize',4); hold on
plot(t_s(1:stopFrame), z_smooth(1:stopFrame), 'r', 'LineWidth',1.5)
xline(0,'--','Color',[0.5 0.5 0.5],'Label','impact','LabelVerticalAlignment','bottom');
ylabel('depth (cm)'); xlim(tRange); legend('raw','SG','Location','northwest'); grid on
title(sprintf('%s T%02d — Kinematics', trialInfo.heightLabel, trialInfo.trialNum));

subplot(3,1,2)
plot(t_s(1:stopFrame), v_smooth(1:stopFrame), 'b', 'LineWidth',1.5)
xline(0,'--','Color',[0.5 0.5 0.5]);
yline(0,'--','Color',[0.65 0.65 0.65],'LineWidth',0.8);
ylabel('v (cm/s)'); xlim(tRange); grid on

subplot(3,1,3)
plot(t_s(1:stopFrame), a_plus_g(1:stopFrame), 'm', 'LineWidth',1.5)
xline(0,'--','Color',[0.5 0.5 0.5]);
yline(0,'--','Color',[0.65 0.65 0.65],'LineWidth',0.8);
ylabel('a+g (cm/s²)'); xlabel('t (s)'); xlim(tRange); grid on
%% Fig 2 — Tracking QA
figure('Name','Tracking QA');
subplot(2,1,1)
for m = 1:nMarkers
    plot(trackedX(m,:), '.-', 'Color',cmap_markers(m,:)); hold on
end
xline(impact_index,'r','impact'); xlabel('frame'); ylabel('x (px)');
title('marker tracks'); grid on
legend(arrayfun(@(i) sprintf('M%d',i),1:nMarkers,'UniformOutput',false), ...
    'Location','eastoutside');
subplot(2,1,2)
nDetected = cellfun(@(c) size(c,1), detections);
bar(nDetected,'FaceColor',[0.3 0.6 0.9],'EdgeColor','none');
xline(impact_index,'r','impact'); xlabel('frame'); ylabel('n detections'); grid on
title('detections per frame');

% Fig 3 — Phase space
figure('Name','Phase Space');
t_norm = max(0,min(1,(t_s-t_s(impact_index))./(t_s(stopFrame)-t_s(impact_index))));
subplot(3,1,1); scatter(z_smooth,v_smooth,20,t_norm,'filled');
colorbar; colormap(jet); xlabel('depth (cm)'); ylabel('v (cm/s)'); grid on
title('v vs depth (color = norm. time)');
subplot(3,1,2); scatter(z_smooth,a_plus_g,20,t_norm,'filled');
colorbar; xlabel('depth (cm)'); ylabel('a+g (cm/s²)'); grid on
title('a+g vs depth');
subplot(3,1,3); scatter(v_smooth.^2,a_plus_g,20,z_smooth,'filled');
colorbar; colormap(parula); xlabel('v² (cm²/s²)'); ylabel('a+g (cm/s²)'); grid on
title('Katsuragi Fig 3a: a+g vs v²');

% Fig 4 — Frame QA
figure('Name','Frame QA');

frameFiles = sort({dir(fullfile(framesDir,'*.png')).name});
if isempty(frameFiles)
    warning('No PNGs in stored framesDir — please relocate.');
    framesDir  = uigetdir(pwd, 'Select frames folder for this trial');
    frameFiles = sort({dir(fullfile(framesDir,'*.png')).name});
end

theta_c     = linspace(0, 2*pi, 60);
frameIdxs   = [1, impact_index, stopFrame];
frameTitles = {sprintf('First frame (%d)', 1), ...
               sprintf('Impact  (%d)',  impact_index), ...
               sprintf('Stop    (%d)',  stopFrame)};

for k = 1:3
    fi = frameIdxs(k);
    subplot(2,3,k);
    imshow(imread(fullfile(framesDir, frameFiles{fi}))); hold on
    plot([bedPoint1(1) bedPoint2(1)],[bedPoint1(2) bedPoint2(2)],'--y','LineWidth',2);
    cF = detections{fi};
    rF = S.det.detect.radiiCell{fi};
    for m = 1:nMarkers
        tx = trackedX(m, fi);
        ty = trackedY(m, fi);
        if ~isfinite(tx) || ~isfinite(ty), continue; end
        col_m = cmap_markers(m,:);
        if ~isempty(cF)
            [~, bestIdx] = min(sqrt((cF(:,1)-tx).^2 + (cF(:,2)-ty).^2));
            r_m = rF(bestIdx);
        else
            r_m = 15;
        end
        plot(tx + r_m.*cos(theta_c), ty + r_m.*sin(theta_c), '-','Color',col_m,'LineWidth',2);
        plot(tx, ty, '+','Color',col_m,'MarkerSize',10,'LineWidth',1.5);
    end
    title(frameTitles{k},'FontSize',10);
end

% Marker x vs time (bottom row, full width)
subplot(2,3,[4 5 6]); hold on
for m = 1:nMarkers
    plot(t_s, trackedX(m,:), '-', 'Color', cmap_markers(m,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('M%d', m));
end
xline(0,        '--','Color',[0.50 0.50 0.50],'LineWidth',1.2,...
    'Label','impact','LabelVerticalAlignment','bottom','FontSize',9);
xline(t_stop_s, '--','Color',[0.80 0.10 0.10],'LineWidth',1.2,...
    'Label','stop',  'LabelVerticalAlignment','bottom','FontSize',9);
grid on;
xlabel('t  (s)'); ylabel('x  (px)');
title('Marker x vs time');
%% 6) SAVE
kinematics = struct( ...
    't_s',           t_s,           'depthRod_cm',  depthRod_cm,  ...
    'z_smooth',      z_smooth,      'v_smooth',     v_smooth,     ...
    'a_smooth',      a_smooth,      'a_plus_g',     a_plus_g,     ...
    'trackedX',      trackedX,      'trackedY',     trackedY,     ...
    'impact_index',  impact_index,  'stopFrame',    stopFrame,    ...
    'postCapFrames', postCapFrames, 'refMarkerID',  refMarkerID,  ...
    'rodBedDist_px', rodBedDist_px, 'sgOrder',      sgOrder,      ...
    'sgWindow',      sgWindow,      'mmPerPx',      mmPerPx,      ...
    'impactDistPx',  impactDistPx,  'bedPoint1',    bedPoint1,    ...
    'bedPoint2',     bedPoint2,     'lineA',        lineA,        ...
    'lineB',         lineB,         'lineC',        lineC);

scalars = struct( ...
    'v0_cm_s',    v0_cm_s,         'd_final_cm', d_final_cm, ...
    't_stop_s',   t_stop_s,        'h_cm',       trialInfo.h_cm, ...
    'H_cm',       trialInfo.h_cm + d_final_cm);

save_trial(trialInfo, kinematics, scalars, [1 2 3 4], ...
    {'kinematics','tracking_qa','phase_space','impact_qa'}, outRoot);