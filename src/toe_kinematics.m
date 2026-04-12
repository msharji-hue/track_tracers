function [t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, impact_index, toeMarkerID, toePx] = toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, dt, mmPerPx, impactDistPx, polyDegree)
% TOE_KINEMATICS  Compute depth, velocity and acceleration of the toe marker.
    if nargin < 9, polyDegree = 4; end

    [~, toeMarkerID]  = max(mean(trackedX, 2, 'omitnan'));
    bedNorm           = sqrt(lineA^2 + lineB^2);
    toePx             = (lineA.*trackedX(toeMarkerID,:) + lineB.*trackedY(toeMarkerID,:) + lineC) ./ bedNorm;
    [~, impact_index] = min(abs(toePx - impactDistPx));

    nFrames     = size(trackedX, 2);
    t_s         = ((0:nFrames-1).*dt) - (impact_index-1).*dt;
    depthRod_cm = (toePx - toePx(impact_index)) .* mmPerPx ./ 10;

    [z_smooth, v_smooth, a_smooth] = poly_kinematics(t_s, depthRod_cm, polyDegree);

    z_smooth = z_smooth(:);
    v_smooth = v_smooth(:);
    a_smooth = a_smooth(:);
    a_smooth(1:impact_index-1) = nan;
end
