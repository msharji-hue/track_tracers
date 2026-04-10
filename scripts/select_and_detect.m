%% select_and_detect.m
clear; close all; clc;

%% 1) CODE PATH
codeDir = '/Users/muhannadalsharji/Documents/track_tracers';
addpath(fullfile(codeDir, 'src'));

%% 2) USER INPUTS — folders + trial info
framesDir  = uigetdir(pwd,  'Select folder containing PNG frames');
resultsDir = uigetdir(pwd,  'Select folder to save results');
folderName = inputdlg('Name your results subfolder:', 'Output', 1, {'detections'});
detDir     = fullfile(resultsDir, folderName{1});
if ~exist(detDir, 'dir'), mkdir(detDir); end

material    = inputdlg({'Material (GB or CHIN)','Batch','Height label','Trial number'}, ...
                        'Trial Info', 1, {'GB','Batch1','H10','1'});
trialInfo   = struct('material',  material{1}, ...
                     'batchName', material{2}, ...
                     'heightLabel', material{3}, ...
                     'trialNum',  str2double(material{4}), ...
                     'fps_true',  2350);

%% 3) DETECTION PARAMETERS
params = struct();
params.radiusRange       = [9 33];
params.sensitivity       = 0.84;
params.edgeThresh        = 0.10;
params.polarity          = 'bright';
params.alphaG            = 0.50;
params.betaB             = 0.50;
params.doCLAHE           = true;
params.medianK           = 3;
params.showPreviewEveryN = 100;
params.showCenters       = true;
params.heightLabel       = trialInfo.heightLabel;

%% 4) RUN DETECTION + SAVE
detectOut = detect_circles_per_frame(framesDir, params);

saveInfo  = struct('material',    trialInfo.material, ...
                   'batchName',   trialInfo.batchName, ...
                   'heightLabel', trialInfo.heightLabel, ...
                   'trialNum',    trialInfo.trialNum, ...
                   'framesDir',   framesDir, ...
                   'fps_export',  trialInfo.fps_true, ...
                   'detDir',      detDir);

det = save_detections(detectOut, saveInfo);