function det = save_detections(detectOut, saveInfo)
% SAVE_DETECTIONS  Save detection results to .mat and .csv.

baseName = sprintf('%s_T%02d', saveInfo.heightLabel, saveInfo.trialNum);
outMat   = fullfile(saveInfo.detDir, [baseName '_detections.mat']);
outCsv   = fullfile(saveInfo.detDir, [baseName '_detections.csv']);

n = detectOut.nFramesReadOK;

% Build output struct
det.meta   = rmfield(saveInfo, 'detDir');
det.frame  = struct('frame0', detectOut.frame0, 'iMat', detectOut.iMat);
det.detect = struct('centersCell', {detectOut.centersCell}, ...
                    'radiiCell',   {detectOut.radiiCell}, ...
                    'nDetected',   detectOut.nDetected, ...
                    'meanRadius',  detectOut.meanRadius);

save(outMat, 'det', '-v7.3');
fprintf('Saved MAT: %s\n', outMat);

% Write CSV
fid = fopen(outCsv, 'w');
fprintf(fid, 'frame0,iMat,x,y,r,nDetected,fps_export\n');
for i = 1:n
    C  = detectOut.centersCell{i};
    R  = detectOut.radiiCell{i};
    nd = detectOut.nDetected(i);
    f0 = detectOut.frame0(i);
    if isempty(C)
        fprintf(fid, '%d,%d,NaN,NaN,NaN,%d,%.6f\n', f0, i, nd, saveInfo.fps_export);
    else
        for k = 1:size(C,1)
            fprintf(fid, '%d,%d,%.6f,%.6f,%.6f,%d,%.6f\n', f0, i, C(k,1), C(k,2), R(k), nd, saveInfo.fps_export);
        end
    end
end
fclose(fid);
fprintf('Saved CSV: %s\n', outCsv);
end
