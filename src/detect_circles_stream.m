function detectOut = detect_circles_stream(v, startFrame, endFrame, filterType, params)
%DETECT_CIRCLES_STREAM  Detect markers straight from the video, writing no PNGs.
%
%   The default Stage A path. Reads [startFrame, endFrame] sequentially from an
%   open VideoReader, applies the same filter export_frames would have applied,
%   runs the same per-frame detection on the in-memory image, and discards the
%   frame. Nothing is written to disk.
%
%   detectOut = detect_circles_stream(v, startFrame, endFrame, filterType, params)
%     v           open VideoReader (from open_video)
%     startFrame  first ABSOLUTE video frame of the window
%     endFrame    last  ABSOLUTE video frame of the window
%     filterType  as apply_filter: 'none'|'grayscale'|'gaussian'|'sharpen'
%     params      detection parameters (see detect_circles_frame)
%
%   Returns the SAME struct detect_circles_per_frame returns, with indices
%   1..N over the window, so every downstream consumer -- firstValidFrame,
%   tracking, save_detections -- is unchanged. The window's absolute position is
%   carried separately as meta.windowStart.
%
%   EQUIVALENCE WITH THE PNG PATH. export_frames applied apply_filter and wrote
%   a PNG; detect_circles_per_frame read that PNG back and detected on it. PNG
%   is lossless and apply_filter returns a uint8 RGB image, so imread(imwrite(x))
%   is x and detecting on the in-memory filtered frame gives bit-identical
%   results. The only difference is that no file is created.
%
%   READ ORDER. CurrentTime is set once and frames are then read sequentially,
%   rather than seeking per frame. Sequential decode is what video containers
%   are built for; per-frame seeking on a 40k-frame clip is far slower and
%   yields the same frames.

    nAvail     = floor(v.Duration * v.FrameRate);
    startFrame = max(1, round(startFrame));
    endFrame   = min(nAvail, round(endFrame));
    if endFrame < startFrame
        error('detect_circles_stream:emptyWindow', ...
              'Window [%d, %d] is empty (video has %d frames).', ...
              startFrame, endFrame, nAvail);
    end
    nFrames = endFrame - startFrame + 1;

    centersCell      = cell(nFrames,1);
    radiiCell        = cell(nFrames,1);
    nDetected        = zeros(nFrames,1);
    meanRadius       = nan(nFrames,1);
    failedFrame_iMat = NaN;
    failedMsg        = "";
    nFramesReadOK    = 0;

    fprintf('Streaming frames %d-%d (%d) from the video; no PNGs written.\n', ...
            startFrame, endFrame, nFrames);

    v.CurrentTime = (startFrame - 1) / v.FrameRate;

    for i = 1:nFrames
        if ~hasFrame(v)
            failedFrame_iMat = i;
            failedMsg = sprintf(['video ended at window frame %d of %d ' ...
                                 '(absolute %d); saving partial results'], ...
                                 i, nFrames, startFrame + i - 1);
            fprintf('\nSTREAM ENDED EARLY. %s\n', failedMsg);
            break
        end
        try
            rgb = readFrame(v);
        catch ME
            failedFrame_iMat = i;
            failedMsg        = string(ME.message);
            fprintf('\nFAILED reading window frame %d (absolute %d).\nReason: %s\n', ...
                    i, startFrame + i - 1, failedMsg);
            break
        end

        % Same filter the PNG path baked into the file.
        rgb = apply_filter(rgb, filterType);

        nFramesReadOK = i;

        [centers, radii] = detect_circles_frame(rgb, params);

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
            title(sprintf('%s: frame %d / %d (abs %d)  |  detected: %d', ...
                params.heightLabel, i, nFrames, startFrame + i - 1, nDetected(i)));
            drawnow;
        end
    end

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
