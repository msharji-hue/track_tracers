function [ang_deg, S] = rod_angle(trackedX, trackedY, varargin)
% ROD_ANGLE  Angle of the marker-centre line vs a reference frame, per frame,
%            plus scalar angle-change summaries for the trial.
%
%   METHOD (deliberately simple):
%     1. Fit ONE straight line through the tracked marker centres in the first
%        valid (all-marker) frame  -> reference line.
%     2. Fit a straight line through the marker centres in every later frame.
%     3. Report the signed angle of that line relative to the reference line.
%
%   The fit is TOTAL LEAST SQUARES (PCA), not polyfit: the rod is near-vertical,
%   where a y = mx + c fit has unbounded slope. TLS is orientation-free.
%
%   WHAT IT MEASURES: a single line through all markers captures ROTATION of the
%   rod axis, not flexure. A rod that bows symmetrically leaves the best-fit line
%   nearly unrotated, so a flat trace shows it did not TILT, not that it stayed
%   straight.
%
%   [ang_deg, S] = rod_angle(trackedX, trackedY, ...)
%     'refFrame'   []   reference frame index (default: first all-marker frame)
%     'minMarkers' 2    minimum valid markers for a frame to be scored
%     'impact'     []   impact frame  \  if both given, scalars are computed over
%     'stop'       []   stop frame    /  impact:stop, and the noise floor from
%                                        the quiet frames after stop
%     'noiseCap'   500  max frames after stop used for the noise estimate
%
%   ang_deg : [1 x nFrames] SIGNED angle vs the reference line, degrees.
%             Sign follows cross(refDir,dir) in pixel coords; because image y
%             increases downward, positive appears clockwise on screen.
%             Reference frame reads exactly 0.
%
%   S (scalars, all degrees):
%     .peak_abs_deg    max |angle| in the window  <- primary "total angle change"
%     .peak_signed_deg the signed value at that peak
%     .net_deg         settled angle at the end of the window (median of last 5)
%     .range_deg       max - min in the window (full excursion)
%     .noise_sd_deg    SD of the angle after stop (rod at rest) = noise floor
%     .peakFrame .refFrame .win .nScored
%
%   Compare peak_abs_deg against noise_sd_deg: a peak within ~2-3 SD of the
%   post-stop scatter is not resolvable rotation.

    p = struct('refFrame', [], 'minMarkers', 2, 'impact', [], 'stop', [], ...
               'noiseCap', 500);
    for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end

    [~, nF] = size(trackedX);
    valid = isfinite(trackedX) & isfinite(trackedY);

    % -- 1) reference line ------------------------------------------------
    if isempty(p.refFrame)
        ff = find(all(valid, 1), 1, 'first');                 % all markers seen
        if isempty(ff)
            ff = find(sum(valid,1) >= p.minMarkers, 1, 'first');
        end
    else
        ff = p.refFrame;
    end
    if isempty(ff)
        error('rod_angle:noReference', ...
              'No frame has at least %d valid markers.', p.minMarkers);
    end
    refU = fit_dir(trackedX(:,ff), trackedY(:,ff));
    if any(~isfinite(refU))
        error('rod_angle:badReference', ...
              'Reference frame %d did not yield a valid line fit.', ff);
    end

    % -- 2-4) per-frame line fit and angle vs reference --------------------
    ang_deg = nan(1, nF);
    for f = 1:nF
        if sum(valid(:,f)) < p.minMarkers, continue; end
        u = fit_dir(trackedX(:,f), trackedY(:,f));
        if any(~isfinite(u)), continue; end
        % PCA eigenvector sign is arbitrary; align with the reference so the
        % angle does not jump by 180 deg between frames.
        if dot(u, refU) < 0, u = -u; end
        % signed angle FROM reference TO current: cross as sine, dot as cosine.
        ang_deg(f) = atan2d(refU(1)*u(2) - refU(2)*u(1), dot(u, refU));
    end

    % -- scalar summaries over the analysis window ------------------------
    %   Default is the whole trial, but when impact/stop are supplied the window
    %   is impact:stop. That matters: max|angle| over a long post-stop tail is a
    %   max over thousands of noise samples and grows with recording length,
    %   which would make trials incomparable.
    if ~isempty(p.impact) && ~isempty(p.stop)
        i1 = max(1, p.impact); i2 = min(nF, p.stop);
    else
        i1 = 1; i2 = nF;
    end
    if i2 < i1, i1 = 1; i2 = nF; end
    win = i1:i2;
    a   = ang_deg(win);
    fin = win(isfinite(a));
    aok = a(isfinite(a));

    S = struct();
    S.refFrame = ff;
    S.win      = [i1 i2];
    S.nScored  = numel(aok);
    if isempty(aok)
        S.peak_abs_deg    = NaN;
        S.peak_signed_deg = NaN;
        S.net_deg         = NaN;
        S.range_deg       = NaN;
        S.peakFrame       = NaN;
    else
        [S.peak_abs_deg, ipk] = max(abs(aok));
        S.peak_signed_deg     = aok(ipk);
        S.peakFrame           = fin(ipk);
        S.range_deg           = max(aok) - min(aok);
        nlast                 = min(5, numel(aok));
        S.net_deg             = median(aok(end-nlast+1:end));
    end

    % -- noise floor from the quiet frames after stop (rod at rest) -------
    S.noise_sd_deg = NaN;
    if ~isempty(p.stop)
        j1 = min(nF, p.stop + 1);
        j2 = min(nF, p.stop + p.noiseCap);
        if j2 > j1
            q = ang_deg(j1:j2); q = q(isfinite(q));
            if numel(q) >= 5, S.noise_sd_deg = std(q); end
        end
    end
end

% -- helper ---------------------------------------------------------------
function u = fit_dir(x, y)
% Principal direction of the point set (total least squares / PCA).
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 2, u = [NaN NaN]; return; end
    xo = x(ok); yo = y(ok);
    P  = [xo(:) - mean(xo), yo(:) - mean(yo)];
    [~, ~, V] = svd(P, 0);
    u = V(:,1).';
end
