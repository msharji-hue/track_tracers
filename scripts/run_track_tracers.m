clear; clc; close all;

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot,'src')));

params = struct();
params.projectRoot = "/Users/muhannadalsharji/Library/CloudStorage/GoogleDrive-msharji@umich.edu/My Drive/Research/Videos/SP26/JerboaImpact_VideoExports";
params.material    = "GB";
params.batchName   = "Batch1";
params.heightLabel = "H06";
params.trialNum    = 1;
params.version     = "short";

results = track_tracers(params);

outDir = fullfile(repoRoot,"outputs");
if ~exist(outDir,'dir'); mkdir(outDir); end
save(fullfile(outDir,"results.mat"),"results");
disp("Done.");
