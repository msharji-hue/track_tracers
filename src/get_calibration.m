function calib = get_calibration(bedX, container, model)
% GET_CALIBRATION  Single global calibration set used by every trial.
%
%   calib = get_calibration()                 % bed x=4, trigger -360
%   calib = get_calibration(10)               % override bed line
%   calib = get_calibration([], 'dense')      % container recorded, no effect
%   calib = get_calibration([], '', 'Tight')  % model recorded, same trigger
%
%   impactDistPx is defined RELATIVE to whichever bed line is in use: at contact
%   the reference marker (furthest from the line on the rod side = top of rod)
%   sits roughly this many px from THAT line.
%
%   *** IT IS A SEARCH-WINDOW CENTRE, NOT THE IMPACT FRAME. ***
%   kd_kinematics uses it only to place the geometric anchor
%       [~, geomImpact] = min(abs(rodBedDist_px - calib.impactDistPx));
%   and then searches [geomImpact - 0.5*wMax, geomImpact + 2*wMax] for the PEAK
%   of the smoothed velocity. Impact is that velocity peak. The trigger's only
%   job is to keep the window over the real event so a spurious pre-release
%   maximum cannot win. It therefore needs to be approximately right, not
%   exactly right, and a single standardized value is safe for every model.
%
%   STANDARDIZED: -360 px for ALL models (2026-08 unified campaign). The earlier
%   per-model values (Default -370, Tight -376.001, Wide -409) tracked each
%   3D-printed foot's rod-to-toe offset. Those offsets are well inside the
%   search window, so resolving them separately bought no accuracy while making
%   every downstream comparison depend on which model wrote the trigger.
%
%   RE-VALIDATE IF THE CAMERA FRAMING CHANGES. The window is defined in pixels
%   from the bed line, so a change of lens, working distance, or crop can move
%   the true impact outside it. Run scripts/diag_impact_frame.m on a few trials
%   per model after any framing change and confirm the located impact frame sits
%   at contact. kd_kinematics also warns (kd_kinematics:triggerOutOfRange) when
%   impactDistPx falls outside the reference marker's observed range.
%
%   The 'model' and 'container' arguments are accepted and recorded for
%   provenance but neither changes the trigger. (An earlier per-container value
%   of -290 for dense was tested and removed.)
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

    % ── Impact search-window centre (relative to the bed line above) ─────
    %   At contact the reference marker (furthest from the bed line on the rod
    %   side = top of the rod) has its centre approximately this far from the
    %   line. Defined RELATIVE to whichever bed line is in use, so it does not
    %   change when bedX changes.
    %
    %   SINGLE VALUE FOR ALL MODELS AND ALL CONTAINERS (2026-08 campaign). This
    %   centres the impact SEARCH window; kd_kinematics then locates impact as
    %   the velocity peak inside that window, so the value only has to be close
    %   enough to bracket contact. Retired: per-model -370/-376.001/-409, and an
    %   earlier per-container -290 for dense.
    %
    %   Trials whose reference marker never reaches -360 may report ANCHOR_OOR;
    %   the velocity-peak refinement still locates impact in that case.
    %
    %   Re-validate with scripts/diag_impact_frame.m if the camera framing
    %   changes -- the window is in pixels from the bed line.
    calib.impactDistPx = -360;       % px, all models, all containers
    calib.model        = model;      % recorded for provenance only; no effect
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
