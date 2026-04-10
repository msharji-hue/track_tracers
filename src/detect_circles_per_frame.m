function detectOut = detect_circles_per_frame(framesDir, params)
% DETECT_CIRCLES_PER_FRAME  Detect red circular markers from a folder of PNGs.
%
%   Inputs:
%       framesDir  - path to folder containing frame_XXXXX.png files
%       params     - detection parameter struct

% Get sorted list of PNGs
files   = dir(fullfile(framesDir, '*.png'));
files   = sort({files.name});
nFrames = numel(files);

% Storage
centersCell      = cell(nFrames,1);
radiiCell        = cell(nFrames,1);
nDetected        = zeros(nFrames,1);
meanRadius       = nan(nFrames,1);
failedFrame_iMat = NaN;
failedMsg        = "";
nFramesReadOK    = 0;

for i = 1:nFrames
    try
        rgb = imread(fullfile(framesDir, files{i}));
    catch ME
        failedFrame_iMat = i;
        failedMsg        = string(ME.message);
        fprintf('\nFAILED reading %s. Saving partial results.\nReason: %s\n', files{i}, failedMsg);
        break
    end

    nFramesReadOK = i;

    % --- Redness image ---
    R = im2double(rgb(:,:,1));
    G = im2double(rgb(:,:,2));
    B = im2double(rgb(:,:,3));
    A = max(min(R - params.alphaG.*G - params.betaB.*B, 1), 0);

    if params.doCLAHE,    A = adapthisteq(A);                              end
    if params.medianK > 1, A = medfilt2(A, [params.medianK params.medianK]); end

    % --- Detect circles ---
    try
        [centers, radii] = imfindcircles(A, params.radiusRange, ...
            'ObjectPolarity', params.polarity, ...
            'Sensitivity',    params.sensitivity, ...
            'EdgeThreshold',  params.edgeThresh);
    catch ME
        fprintf('Detection failed at frame %d: %s\n', i, ME.message);
        centers = []; radii = [];
    end

    centersCell{i} = centers;
    radiiCell{i}   = radii;
    nDetected(i)   = size(centers,1);
    if ~isempty(radii), meanRadius(i) = mean(radii); end

    if mod(i,200)==0 || i==nFrames
        fprintf('Processed %d / %d frames\n', i, nFrames);
    end

    if params.showPreviewEveryN > 0 && (mod(i,params.showPreviewEveryN)==0 || i==1)
        figure(10); clf; imshow(rgb); hold on;
        if ~isempty(centers)
            viscircles(centers, radii, 'Color', 'y');
            if params.showCenters
                plot(centers(:,1), centers(:,2), 'b+', 'LineWidth', 1.0);
            end
        end
        title(sprintf('%s: frame %d / %d  |  detected: %d', ...
            params.heightLabel, i, nFrames, nDetected(i)));
        drawnow;
    end
end

% Trim to successfully read frames
n = nFramesReadOK;
detectOut = struct( ...
    'centersCell',      {centersCell(1:n)}, ...
    'radiiCell',        {radiiCell(1:n)}, ...
    'nDetected',        nDetected(1:n), ...
    'meanRadius',       meanRadius(1:n), ...
    'frame0',           (0:n-1)', ...
    'iMat',             (1:n)', ...
    'nFramesReadOK',    n, ...
    'failedFrame_iMat', failedFrame_iMat, ...
    'failedMsg',        failedMsg);
end
