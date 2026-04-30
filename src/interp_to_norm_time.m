function [tPhys, yMean, yStd, nValid] = interp_to_norm_time(hg, field, nGrid)

    if nargin < 3, nGrid = 200; end

    tNorm_grid = linspace(0, 1, nGrid);
    nT         = hg.nTrials;
    yMat       = nan(nT, nGrid);

    for i = 1:nT
        k  = hg.trials(i).kinematics;

        % Force all to column vectors to prevent size mismatch
        t  = k.t_s(:);
        y  = k.(field)(:);
        t0 = t(k.impact_index);
        t1 = t(k.stopFrame);

        tNorm = (t - t0) ./ (t1 - t0);

        % Only use impact → stop window with finite y
        mask = tNorm >= 0 & tNorm <= 1 & isfinite(y);
        if sum(mask) < 3, continue; end

        yMat(i,:) = interp1(tNorm(mask), y(mask), tNorm_grid, 'pchip', nan);
    end

    % Stats where >= 3 trials valid
    nValid    = sum(isfinite(yMat), 1);
    ok        = nValid >= 3;
    yMean     = nan(1, nGrid);
    yStd      = nan(1, nGrid);
    yMean(ok) = mean(yMat(:, ok), 1, 'omitnan');
    yStd(ok)  = std( yMat(:, ok), 0, 1, 'omitnan');

    % Rescale to physical time using mean t_stop
    tPhys = tNorm_grid .* hg.tstop_mean;
end