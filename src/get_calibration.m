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
    calib.mmPerPx      = 0.1079;     % mm per pixel
    calib.impactDistPx = -400;       % signed toe-to-bed distance (px) at impact

    % ── Bed line (two points, image px) ──────────────────────────────────
    calib.bedPoint1    = [23,  0];
    calib.bedPoint2    = [23, 32];

    % Implicit line  A*x + B*y + C = 0  through the two bed points
    calib.lineA = calib.bedPoint1(2) - calib.bedPoint2(2);
    calib.lineB = calib.bedPoint2(1) - calib.bedPoint1(1);
    calib.lineC = calib.bedPoint1(1)*calib.bedPoint2(2) - ...
                  calib.bedPoint2(1)*calib.bedPoint1(2);

    % ── Physics / timing ─────────────────────────────────────────────────
    calib.g_cm_s2   = 980;           % gravity, cm/s^2
    calib.postCapMs = 10;            % ms after t_stop kept in saved tracks
end
