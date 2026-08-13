function [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, field, nGrid)
% INTERP_TO_NORM_TIME  Interpolate trials onto normalized time grid.
%   Grid spans pre-impact through stop: tNorm in [-tPre, 1]
%   tNorm = 0 at impact, tNorm = 1 at stop.

    if nargin < 3, nGrid = 300; end

    % Find longest pre-impact window across trials (in normalized units)
    preMax = 0;
    for i = 1:hg.nTrials
        k      = hg.trials(i).kinematics;
        t0     = k.t_s(k.impact_index);
        t1     = k.t_s(k.stopFrame);
        tspan  = t1 - t0;
        tpre   = t0 - k.t_s(1);          % time before impact
        preMax = max(preMax, tpre/tspan); % normalized pre-impact length
    end
    preMax = min(preMax, 0.5);   % cap at half a penetration-time worth

    tNorm_grid = linspace(-preMax, 1, nGrid);
    nT         = hg.nTrials;
    yMat       = nan(nT, nGrid);

    for i = 1:nT
        k  = hg.trials(i).kinematics;
        t  = k.t_s(:);
        y  = k.(field)(:);
        t0 = t(k.impact_index);
        t1 = t(k.stopFrame);

        tNorm = (t - t0) ./ (t1 - t0);

        % Use all frames from start through stop
        mask = tNorm >= -preMax & tNorm <= 1 & isfinite(y);
        if sum(mask) < 3, continue; end

        yMat(i,:) = interp1(tNorm(mask), y(mask), tNorm_grid, 'pchip', nan);
    end

    % Stats where >= 3 trials valid
    nValid    = sum(isfinite(yMat), 1);
    ok        = nValid >= 1;
    yMean     = nan(1, nGrid);
    yStd      = nan(1, nGrid);
    yMean(ok) = mean(yMat(:, ok), 1, 'omitnan');
    yStd(ok)  = std( yMat(:, ok), 0, 1, 'omitnan');

    % Rescale to physical time using mean t_stop
    tPhys = tNorm_grid .* hg.tstop_mean;
end