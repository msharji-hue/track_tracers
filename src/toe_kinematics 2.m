% ============================================================================
%  RETIRED BASELINE — not in the pipeline. This is the old Savitzky-Golay
%  position-first, double-differentiated chain. Superseded by kd_kinematics.m
%  (velocity-first, KD 2007). Kept ONLY for the one-time SG vs velocity-first
%  A/B noise comparison for Methods; archive to deprecated/ once that is done.
% ============================================================================
function [t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, ...
          impact_index, refMarkerID, rodBedDist_px, stopFrame, sgOrder, sgWindow] = ...
          toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, ...
                         dt, mmPerPx, impactDistPx, varargin)
% TOE_KINEMATICS  Rod kinematics from tracked marker positions.
%
%   IMPORTANT (naming): the markers are on the RIGID ROD attached to the foot.
%   NONE of them is the foot/toe. Earlier versions of this file called the
%   marker with the largest x the 'toe marker' — that was wrong twice over:
%   (a) no marker is the toe, and (b) with the bed line at low x, max(x) is the
%   marker FARTHEST from the bed. Depth is now the rigid-rod displacement
%   averaged over ALL visible markers (see rod_displacement.m), which is
%   composition-invariant and ~sqrt(N) less noisy than any single marker.
%
%   refMarkerID   : reference marker index (metadata only; NOT 'the toe')
%   rodBedDist_px : signed bed-normal distance of the reference marker (px)
%
%   z_smooth, v_smooth: pre-impact SG fit retained for display.
%   a_smooth: masked to NaN before impact and after stopFrame (meaningless).
%   Analysis functions should index impact_index:stopFrame.
%
%   Optional name-value args:
%       'sgWindow'      - SG filter window (odd integer; default 21; [] = auto)
%       'sgOrder'       - SG filter order (default 4)
%       'postCapFrames' - frames past stopFrame to retain in z/v (default 0)

    p = inputParser;
    addParameter(p, 'sgWindow',      [], @(x) isempty(x) || (isnumeric(x) && x > 0));
    addParameter(p, 'sgOrder',       2, @isnumeric);
    addParameter(p, 'postCapFrames',  0, @isnumeric);
    parse(p, varargin{:});

    sgOrder       = p.Results.sgOrder;
    sgWindow_in   = p.Results.sgWindow;
    postCapFrames = round(p.Results.postCapFrames);

    % ── Reference marker: NEAREST the bed line (metadata / impact timing) ──
    %   Markers are on the rigid rod; none is the toe. We use the marker
    %   closest to the bed purely as the reference for locating impact.
    bedNorm    = sqrt(lineA^2 + lineB^2);
    d_px_all   = (lineA .* trackedX + lineB .* trackedY + lineC) ./ bedNorm;
    firstFrame = find(all(isfinite(trackedX), 1), 1, 'first');
    if isempty(firstFrame)
        firstFrame = find(any(isfinite(trackedX), 1), 1, 'first');
    end
    col = d_px_all(:, firstFrame);
    col(~isfinite(col)) = Inf;
    [~, refMarkerID] = min(abs(col - impactDistPx));   % nearest the impact plane
    rodBedDist_px    = d_px_all(refMarkerID, :);
    [~, impact_index] = min(abs(rodBedDist_px - impactDistPx));

    % ── Depth: rigid-rod displacement averaged over ALL visible markers ────
    %   Composition-invariant (no step when a marker drops out mid-impact) and
    %   ~sqrt(N) less noisy than a single marker. See rod_displacement.m.
    nFrames     = size(trackedX, 2);
    t_s         = ((0:nFrames-1) .* dt) - (impact_index-1) .* dt;
    z_rod_cm    = rod_displacement(trackedX, trackedY, lineA, lineB, lineC, mmPerPx);
    depthRod_cm = z_rod_cm - z_rod_cm(impact_index);   % zero at impact

    % ── SG filter window ──────────────────────────────────────────────────
    validIdx = find(isfinite(depthRod_cm));
    z_valid  = depthRod_cm(validIdx);
    nValid   = numel(z_valid);

    if ~isempty(sgWindow_in)
        sgWindow = sgWindow_in;
        if mod(sgWindow, 2) == 0
            sgWindow = sgWindow + 1;
            warning('sgWindow must be odd — adjusted to %d', sgWindow);
        end
        if sgWindow > nValid
            error('sgWindow (%d) exceeds valid signal length (%d)', sgWindow, nValid);
        end
    else
        sgWindow = max(7, round(0.10 * nValid));
        if mod(sgWindow, 2) == 0, sgWindow = sgWindow + 1; end
    end

    fprintf('SG: order=%d  window=%d frames (%.1f ms)  nValid=%d\n', ...
        sgOrder, sgWindow, sgWindow*dt*1000, nValid);

    % ── Filter ────────────────────────────────────────────────────────────
    z_filt = sgolayfilt(z_valid,               sgOrder, sgWindow);
    v_filt = sgolayfilt(gradient(z_filt, dt),  sgOrder, sgWindow);
    a_filt = sgolayfilt(gradient(v_filt, dt),  sgOrder, sgWindow);

    % ── Initialize all outputs to NaN then fill valid frames ──────────────
    z_smooth = nan(nFrames, 1);
    v_smooth = nan(nFrames, 1);
    a_smooth = nan(nFrames, 1);
    z_smooth(validIdx) = z_filt;
    v_smooth(validIdx) = v_filt;
    a_smooth(validIdx) = a_filt;

    % ── Stop frame ────────────────────────────────────────────────────────
    % (a+g) per unit mass in this depth-sign convention is (-a - g); it is
    % large and positive at impact (strong resistive deceleration) and falls
    % to zero at the end of penetration. Find that first post-impact crossing.
    % NOTE: a_smooth here is still unmasked (masking happens below), which is
    % what the stop search needs. g matches the project convention (980 cm/s^2).
    g_cm_s2       = 980;
    a_plus_g_full = -a_smooth - g_cm_s2;
    postImpact_ag = a_plus_g_full(impact_index:end);
    stopIdx = find(postImpact_ag <= 0, 1, 'first');
    if isempty(stopIdx)
        stopFrame = nFrames;
        warning('Stop frame not found — using last frame.');
    else
        stopFrame = stopIdx + impact_index - 1;
    end

    fprintf('impact: %d  |  stop: %d  |  t_stop: %.4f s\n', ...
        impact_index, stopFrame, t_s(stopFrame));

    % ── Mask ──────────────────────────────────────────────────────────────
    % a: strict — pre-impact and post-stop acceleration is meaningless
    a_smooth(1:impact_index-1) = nan;
    a_smooth(stopFrame+1:end)  = nan;

    % z, v: pre-impact SG fit kept; clip at stop + postCapFrames
    vidEnd = min(stopFrame + postCapFrames, nFrames);
    z_smooth(vidEnd+1:end) = nan;
    v_smooth(vidEnd+1:end) = nan;
end