function calib = get_calibration(bedX, container, model)
% GET_CALIBRATION  Single global calibration set used by every trial.
%
%   calib = get_calibration()                 % bed x=4, trigger -370
%   calib = get_calibration(10)               % override bed line
%   calib = get_calibration([], 'dense')      % container recorded, no effect
%   calib = get_calibration([], '', 'Tight')  % per-MODEL impact trigger
%
%   impactDistPx is defined RELATIVE to whichever bed line is in use: at contact
%   the reference marker (furthest from the line on the rod side = top of rod)
%   sits this many px from THAT line. It is NOT a fitted parameter, and it is
%   the SAME for every container, but it DOES depend on the foot MODEL: each
%   3D-printed foot has a different rod-to-toe offset, so the top marker sits a
%   different distance from the bed at contact. Values measured before testing:
%       Default  -370      Tight  -376.001      Wide  -409
%
%   The 'container' argument is accepted and recorded for provenance but no
%   longer changes the trigger. (An earlier per-container value of -290 for
%   dense was tested and removed.)
%
%   Consequence: a bed-line shift is not neutral for the ANCHOR (the trigger
%   moves with the line), but it cancels out of DEPTH, which is measured
%   relative to the impact frame: z - z(impact_index).

    if nargin < 3 || isempty(model),     model = 'Default'; end
    if nargin < 2 || isempty(container), container = ''; end
    if nargin < 1 || isempty(bedX),      bedX = 4; end   % x = 4 for ALL conditions

    calib = struct();

    % ── Pixel scale & impact reference ────────────────────────────────────
    calib.mmPerPx      = 0.1079;     % mm per pixel. 2 mm marker = 18.5 px,
                                     %  pxPerMm = 9.27. Settled by direct
                                     %  measurement. Supersedes 0.1429 and a
                                     %  briefly-trialled 0.12245.

    % ── Bed line (two points, image px) ──────────────────────────────────
    calib.bedX         = bedX;
    calib.bedPoint1    = [bedX,  0];
    calib.bedPoint2    = [bedX, 36];

    % ── Impact trigger (relative to the bed line above) ──────────────────
    %   At contact the reference marker (furthest from the bed line on the rod
    %   side = top of the rod) has its centre this far from the line. Defined
    %   RELATIVE to whichever bed line is in use, so it does not change when
    %   bedX changes.
    %
    %   SINGLE VALUE FOR ALL CONDITIONS. An earlier per-container variant used
    %   -290 for dense; that has been retired in favour of one uniform trigger.
    %   Note that dense trials may therefore report ANCHOR_OOR, since the top
    %   marker in those trials often does not reach -370 -- the velocity-peak
    %   refinement still locates impact in that case.
    switch lower(strtrim(char(model)))
        case 'tight', calib.impactDistPx = -376.001;
        case 'wide',  calib.impactDistPx = -409;
        otherwise,    calib.impactDistPx = -370;     % Default geometry
    end
    calib.model        = model;      % geometry that set the trigger
    calib.container    = container;  % recorded for provenance only; no effect

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
