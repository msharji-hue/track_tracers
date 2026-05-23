function heights = group_trials_by_height(trials)
% GROUP_TRIALS_BY_HEIGHT  Group loaded trials by h_cm.
%   Snaps h_cm to nearest standard height before grouping
%   to handle minor metadata discrepancies (e.g. 35.36 → 35.56 cm).
%
%   Output fields per height group:
%       h_cm           - drop height (cm)
%       trials         - array of trial structs
%       nTrials        - number of trials
%       v0_mean/std    - impact velocity stats (cm/s)
%       d_mean/std     - final depth stats (cm)
%       tstop_mean/std - stop time stats (s)
%       v0_cm_s        - all v0 values [1 x nTrials]
%       d_final_cm     - all d_final values [1 x nTrials]

    % Standard heights: 6" through 18" in 2" steps
   standardH = [7.62, 10.16, 12.70, 15.24, 17.78, 20.32, 22.86, ...
             25.40, 27.94, 30.48, 33.02, 35.56, 38.10, 40.64, 43.18, 45.72];
    % Snap each trial to nearest standard height
    for i = 1:numel(trials)
        [~, idx]       = min(abs(trials(i).h_cm - standardH));
        trials(i).h_cm = standardH(idx);
    end

    hVals   = unique([trials.h_cm]);
    heights = struct([]);

    for j = 1:numel(hVals)
        h   = hVals(j);
        idx = find(abs([trials.h_cm] - h) < 0.5);

        v0s    = arrayfun(@(i) trials(i).scalars.v0_cm_s,    idx);
        ds     = arrayfun(@(i) trials(i).scalars.d_final_cm, idx);
        tstops = arrayfun(@(i) trials(i).scalars.t_stop_s,   idx);

        heights(j).h_cm       = h;
        heights(j).trials     = trials(idx);
        heights(j).nTrials    = numel(idx);
        heights(j).v0_mean    = mean(v0s);
        heights(j).v0_std     = std(v0s);
        heights(j).d_mean     = mean(ds);
        heights(j).d_std      = std(ds);
        heights(j).tstop_mean = mean(tstops);
        heights(j).tstop_std  = std(tstops);
        heights(j).v0_cm_s    = v0s;
        heights(j).d_final_cm = ds;

        fprintf('h = %5.2f cm | n = %d | v0 = %6.1f ± %4.1f cm/s | d = %.3f ± %.3f cm\n', ...
            h, numel(idx), mean(v0s), std(v0s), mean(ds), std(ds));
    end
end