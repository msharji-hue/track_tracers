function [t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, impact_index, toeMarkerID, toePx] = toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, dt, mmPerPx, impactDistPx, polyDegree)
% TOE_KINEMATICS  Compute depth, velocity and acceleration of the toe marker.
%
%   Toe marker = rightmost marker (largest mean x = furthest from bed line).
%   z = 0 at impact, positive into the bed. Matches Katsuragi & Durian convention.
%
%   Inputs:
%       trackedX/Y   - M x nFrames tracked positions (px)
%       lineA/B/C    - implicit bed line coefficients (Ax + By + C = 0)
%       dt           - frame interval (s) = 1/fps
%       mmPerPx      - calibration (mm per pixel)
%       impactDistPx - distance from bed line at impact (px)
%       polyDegree   - polynomial degree for fit (default 4)
%
%   Outputs:
%       t_s          - time vector, t=0 at impact (s)
%       depthRod_cm  - depth of toe (cm), +ve into bed
%       z_smooth     - smoothed depth from poly fit (cm)
%       v_smooth     - velocity (cm/s)
%       a_smooth     - acceleration (cm/s²)
%       impact_index - frame index of impact
%       toeMarkerID  - index of toe marker
%       toePx        - signed distance of toe from bed line (px)

    if nargin < 9, polyDegree = 4; end

    [~, toeMarkerID]  = max(mean(trackedX, 2, 'omitnan'));
    bedNorm           = sqrt(lineA^2 + lineB^2);
    toePx             = (lineA.*trackedX(toeMarkerID,:) + lineB.*trackedY(toeMarkerID,:) + lineC) ./ bedNorm;
    [~, impact_index] = min(abs(toePx - impactDistPx));

    nFrames     = size(trackedX, 2);
    t_s         = ((0:nFrames-1).*dt) - (impact_index-1).*dt;
    depthRod_cm = (toePx - toePx(impact_index)) .* mmPerPx ./ 10;

    % Fit from frame 1 to stop (where depth plateaus)
    vTemp     = gradient(depthRod_cm, dt);
    stopFrame = find(vTemp <= 0, 1, 'first');
    if isempty(stopFrame), stopFrame = nFrames; end

    fitRange = 1:stopFrame;
    [z_smooth, v_smooth, a_smooth] = poly_kinematics(t_s(fitRange), depthRod_cm(fitRange), polyDegree);

    % Ensure column vectors before padding
    z_smooth = z_smooth(:);
    v_smooth = v_smooth(:);
    a_smooth = a_smooth(:);

    % Pad post-stop with NaN
    pad      = nan(nFrames - stopFrame, 1);
    z_smooth = [z_smooth; pad];
    v_smooth = [v_smooth; pad];
    a_smooth = [a_smooth; pad];

    % NaN pre-impact acceleration only
    a_smooth(1:impact_index-1) = nan;
end