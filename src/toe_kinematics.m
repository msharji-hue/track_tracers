function [t_s, depthRod_cm, z_smooth, v_smooth, a_smooth, ...
          impact_index, toeMarkerID, toePx, stopFrame, sgOrder, sgWindow] = ...
          toe_kinematics(trackedX, trackedY, lineA, lineB, lineC, ...
                         dt, mmPerPx, impactDistPx, varargin)

    % Optional name-value args
    p = inputParser;
    addParameter(p, 'sgWindow', 71);   % empty = auto
    addParameter(p, 'sgOrder',  3);
    parse(p, varargin{:});
    sgOrder     = p.Results.sgOrder;
    sgWindow_in = p.Results.sgWindow;

    % Toe marker: highest x at first valid frame 
    firstFrame   = find(any(isfinite(trackedX), 1), 1, 'first');
    [~, toeMarkerID] = max(trackedX(:, firstFrame), [], 'omitnan');

    % Project onto bed-normal axis
    bedNorm           = sqrt(lineA^2 + lineB^2);
    toePx             = (lineA .* trackedX(toeMarkerID,:) + ...
                         lineB .* trackedY(toeMarkerID,:) + lineC) ./ bedNorm;
    [~, impact_index] = min(abs(toePx - impactDistPx));

    nFrames     = size(trackedX, 2);
    t_s         = ((0:nFrames-1) .* dt) - (impact_index-1) .* dt;
    depthRod_cm = (toePx - toePx(impact_index)) .* mmPerPx ./ 10;

    % SG window: manual or adaptive
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
        sgWindow = max(7, round(0.15 * nValid));
        if mod(sgWindow, 2) == 0, sgWindow = sgWindow + 1; end
    end

    fprintf('SG: order=%d  window=%d frames (%.1f ms)  nValid=%d\n', ...
        sgOrder, sgWindow, sgWindow*dt*1000, nValid);

    % ── Filter ────────────────────────────────────────────────────────────
    z_filt = sgolayfilt(z_valid,                sgOrder, sgWindow);
    v_filt = sgolayfilt(gradient(z_filt, dt),   sgOrder, sgWindow);
    a_filt = sgolayfilt(gradient(v_filt, dt),   sgOrder, sgWindow);

    z_smooth = nan(nFrames, 1);
    v_smooth = nan(nFrames, 1);
    a_smooth = nan(nFrames, 1);
    z_smooth(validIdx) = z_filt;
    v_smooth(validIdx) = v_filt;
    a_smooth(validIdx) = a_filt;

    % ── Stop frame ────────────────────────────────────────────────────────
    postImpact = v_smooth(impact_index:end);
    stopIdx    = find(postImpact <= 0, 1, 'first');
    if isempty(stopIdx)
        stopFrame = nFrames;
        warning('Stop frame not found — using last frame.');
    else
        stopFrame = stopIdx + impact_index - 1;
    end

    fprintf('impact: %d  |  stop: %d  |  t_stop: %.4f s\n', ...
        impact_index, stopFrame, t_s(stopFrame));

    % ── Mask pre-impact and post-stop ─────────────────────────────────────
    a_smooth(1:impact_index-1) = nan;
    a_smooth(stopFrame+1:end)  = nan;
    z_smooth(stopFrame+1:end)  = nan;
    v_smooth(stopFrame+1:end)  = nan;
end