function [t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, impact_index, toeMarkerID, toePx] = toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, dt, mmPerPx, impactDistPx)
% TOE_KINEMATICS  Compute depth, velocity and acceleration of the toe marker.

    [~, toeMarkerID]  = max(mean(trackedX, 2, 'omitnan'));
    bedNorm           = sqrt(lineA^2 + lineB^2);
    toePx             = (lineA.*trackedX(toeMarkerID,:) + lineB.*trackedY(toeMarkerID,:) + lineC) ./ bedNorm;
    [~, impact_index] = min(abs(toePx - impactDistPx));

    nFrames     = size(trackedX, 2);
    t_s         = ((0:nFrames-1).*dt) - (impact_index-1).*dt;
    depthRod_cm = (toePx - toePx(impact_index)) .* mmPerPx ./ 10;

   % Savitzky-Golay parameters
    sgOrder  = 4;
    sgWindow = 111;

    % Apply filter only over valid (non-NaN) region to avoid edge distortion
    validIdx = find(isfinite(depthRod_cm));
    z_valid  = depthRod_cm(validIdx);

    % Ensure window doesn't exceed signal length (must be odd)
    sgWindow = min(sgWindow, numel(z_valid));
    if mod(sgWindow, 2) == 0, sgWindow = sgWindow - 1; end

    z_filt = sgolayfilt(z_valid,                sgOrder, sgWindow);
    v_filt = sgolayfilt(gradient(z_filt, dt),   sgOrder, sgWindow);
    a_filt = sgolayfilt(gradient(v_filt, dt),   sgOrder, sgWindow);

    z_smooth = nan(nFrames, 1);
    v_smooth = nan(nFrames, 1);
    a_smooth = nan(nFrames, 1);
    z_smooth(validIdx) = z_filt;
    v_smooth(validIdx) = v_filt;
    a_smooth(validIdx) = a_filt;

    % Find stop frame (after impact, where velocity first goes negative)
    stopFrame = find(v_smooth(impact_index:end) <= 0, 1, 'first') + impact_index - 1;
    if isempty(stopFrame), stopFrame = nFrames; end
    fprintf('impact: %d  |  stop: %d  |  t_stop: %.4f s\n', impact_index, stopFrame, t_s(stopFrame));

    % NaN pre-impact acceleration and post-stop everything
    a_smooth(1:impact_index-1) = nan;
    a_smooth(stopFrame+1:end)  = nan;
    z_smooth(stopFrame+1:end)  = nan;
    v_smooth(stopFrame+1:end)  = nan;
end