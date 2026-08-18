function [centers, radii] = detect_circles_frame(rgb, params)
%DETECT_CIRCLES_FRAME  Detect red circular markers in ONE in-memory RGB frame.
%
%   The per-frame detection core, factored out so the folder-of-PNGs path
%   (detect_circles_per_frame) and the streaming path (detect_circles_stream)
%   run literally the same arithmetic. Any change here reaches both.
%
%   [centers, radii] = detect_circles_frame(rgb, params)
%     rgb     uint8 RGB frame, already filtered (see apply_filter)
%     params  radiusRange, sensitivity, edgeThresh, polarity, alphaG, betaB,
%             doCLAHE, medianK
%
%   A detection failure returns empty centers/radii and prints, rather than
%   throwing, so one bad frame cannot abort a whole trial.

    % --- Redness image ---
    R = im2double(rgb(:,:,1));
    G = im2double(rgb(:,:,2));
    B = im2double(rgb(:,:,3));
    A = max(min(R - params.alphaG.*G - params.betaB.*B, 1), 0);

    if params.doCLAHE,     A = adapthisteq(A);                               end
    if params.medianK > 1, A = medfilt2(A, [params.medianK params.medianK]); end

    % --- Detect circles ---
    try
        [centers, radii] = imfindcircles(A, params.radiusRange, ...
            'ObjectPolarity', params.polarity, ...
            'Sensitivity',    params.sensitivity, ...
            'EdgeThreshold',  params.edgeThresh);
    catch ME
        fprintf('Detection failed on a frame: %s\n', ME.message);
        centers = []; radii = [];
    end
end
