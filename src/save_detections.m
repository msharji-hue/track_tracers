function det = save_detections(detectOut, saveInfo)
% Save detection outputs to MAT and CSV.

trialID  = sprintf('T%02d', saveInfo.trialNum);
baseName = sprintf('%s_%s_long', saveInfo.heightLabel, trialID);

outMat = fullfile(saveInfo.detDir, [baseName '_detections.mat']);
outCsv = fullfile(saveInfo.detDir, [baseName '_detections.csv']);

nSave = detectOut.nFramesReadOK;

centersCell = detectOut.centersCell(1:nSave);
radiiCell   = detectOut.radiiCell(1:nSave);
nDetected   = detectOut.nDetected(1:nSave);
meanRadius  = detectOut.meanRadius(1:nSave);
frame0      = (0:nSave-1)';
iMat        = frame0 + 1;

det = struct();

det.meta.material  = saveInfo.material;
det.meta.batch     = saveInfo.batchName;
det.meta.height    = saveInfo.heightLabel;
det.meta.trialNum  = saveInfo.trialNum;
det.meta.videoFile = saveInfo.videoFile;
det.meta.videoPath = saveInfo.videoPath;

det.export.fps_export = saveInfo.fps_export;
det.export.nFrames    = nSave;

det.frame.frame0 = frame0;
det.frame.iMat   = iMat;

det.detect.centersCell = centersCell;
det.detect.radiiCell   = radiiCell;
det.detect.nDetected   = nDetected;
det.detect.meanRadius  = meanRadius;

save(outMat, 'det', '-v7.3');
fprintf('Saved MAT detections: %s\n', outMat);

fid = fopen(outCsv, 'w');
fprintf(fid, 'frame0,iMat,x,y,r,nDetected,fps_export\n');

for i = 1:nSave
    f0 = frame0(i);
    C  = centersCell{i};
    R  = radiiCell{i};
    nd = nDetected(i);

    if isempty(C)
        fprintf(fid, '%d,%d,NaN,NaN,NaN,%d,%.6f\n', ...
            f0, i, nd, saveInfo.fps_export);
    else
        for k = 1:size(C,1)
            fprintf(fid, '%d,%d,%.6f,%.6f,%.6f,%d,%.6f\n', ...
                f0, i, C(k,1), C(k,2), R(k), nd, saveInfo.fps_export);
        end
    end
end

fclose(fid);
fprintf('Saved CSV detections: %s\n', outCsv);
end