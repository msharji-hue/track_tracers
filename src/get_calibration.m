function calib = get_calibration()
% GET_CALIBRATION  Single, global calibration set used by every trial.
%
%   Per current agreement, calibration is the SAME for all videos
%   (material / container / batch). If you later confirm it varies, give
%   this function an input (e.g. a metadata struct) and branch on it —
%   every caller already routes through here, so nothing else changes.
%
%   Values mirror the constants previously hard-coded in track_tracers_2.m.

    calib = struct();

    % ── Pixel scale & impact reference ────────────────────────────────────
    calib.mmPerPx      = 0.1429;     % mm per pixel (2 mm = 14 px -> 0.142857;
                                     %  supersedes 0.1079. Confirm via the
                                     %  pre-impact g_eff check, validate_calibration_gcheck.m)
    calib.impactDistPx = -400;       % signed rod-to-bed distance (px) at impact

    % ── Bed line (two points, image px) ──────────────────────────────────
    calib.bedPoint1    = [4,  0];
    calib.bedPoint2    = [4, 36];

    % Implicit line  A*x + B*y + C = 0  through the two bed points
    calib.lineA = calib.bedPoint1(2) - calib.bedPoint2(2);
    calib.lineB = calib.bedPoint2(1) - calib.bedPoint1(1);
    calib.lineC = calib.bedPoint1(1)*calib.bedPoint2(2) - ...
                  calib.bedPoint2(1)*calib.bedPoint1(2);

    % ── Tracking ─────────────────────────────────────────────────────────
    %   Max match distance from a marker's last known position. Must stay
    %   BELOW one inter-marker spacing (~34.6 px measured) or a stale
    %   reference can grab a neighbour's detection. Largest observed real
    %   motion is ~5-6 px/frame, so 25 px leaves ample headroom.
    calib.trackTolerancePx = 25;

    % ── Physics / timing ─────────────────────────────────────────────────
    calib.g_cm_s2   = 980;           % gravity, cm/s^2
    calib.postCapMs = 10;            % ms after t_stop kept in saved tracks
end
