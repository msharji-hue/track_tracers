function detectOut = detect_circles_per_frame(v, nFrames, params)
% Detect red circular markers frame-by-frame from a video.

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
        rgb = readFrame(v);
    catch ME
        failedFrame_iMat = i;
        failedMsg = string(ME.message);
        fprintf('\nFAILED reading at iMat=%d (frame0=%d). Saving partial detections and moving on.\n', i, i-1);
        fprintf('Reason: %s\n', failedMsg);
        break
    end

    nFramesReadOK = i;

    % --- Convert to channels ---
    R = im2double(rgb(:,:,1));
    G = im2double(rgb(:,:,2));
    B = im2double(rgb(:,:,3));

    % --- Redness image ---
    A = R - params.alphaG .* G - params.betaB .* B;
    A = max(min(A,1),0);

    if params.doCLAHE
        A = adapthisteq(A);
    end
    if params.medianK > 1
        A = medfilt2(A, [params.medianK params.medianK]);
    end

    % --- Detect circles ---
    try
        [centers, radii] = imfindcircles(A, params.radiusRange, ...
            'ObjectPolarity', params.polarity, ...
            'Sensitivity',   params.sensitivity, ...
            'EdgeThreshold', params.edgeThresh);
    catch ME
        fprintf('Detection failed at iMat=%d (frame0=%d): %s\n', i, i-1, ME.message);
        centers = [];
        radii   = [];
    end

    centersCell{i} = centers;
    radiiCell{i}   = radii;
    nDetected(i)   = size(centers,1);

    if ~isempty(radii)
        meanRadius(i) = mean(radii);
    end

    if mod(i,200)==0 || i==nFrames
        fprintf('Detected circles: %d / %d frames\n', i, nFrames);
    end

    if params.showPreviewEveryN > 0 && (mod(i,params.showPreviewEveryN)==0 || i==1)
        figure(10); clf;
        imshow(rgb); hold on;

        if ~isempty(centers)
            viscircles(centers, radii, 'Color', 'y');
            if params.showCenters
                plot(centers(:,1), centers(:,2), 'b+', 'LineWidth', 1.0);
            end
        end

        title(sprintf('%s: frame0=%d (i=%d), nDetected=%d', ...
            params.heightLabel, i-1, i, nDetected(i)));
        drawnow;
    end
end

% Trim to successfully read frames
centersCell = centersCell(1:nFramesReadOK);
radiiCell   = radiiCell(1:nFramesReadOK);
nDetected   = nDetected(1:nFramesReadOK);
meanRadius  = meanRadius(1:nFramesReadOK);
frame0      = (0:nFramesReadOK-1)';
iMat        = frame0 + 1;

% Output struct
detectOut = struct();
detectOut.centersCell      = centersCell;
detectOut.radiiCell        = radiiCell;
detectOut.nDetected        = nDetected;
detectOut.meanRadius       = meanRadius;
detectOut.frame0           = frame0;
detectOut.iMat             = iMat;
detectOut.nFramesReadOK    = nFramesReadOK;
detectOut.failedFrame_iMat = failedFrame_iMat;
detectOut.failedMsg        = failedMsg;
end