function [Fz_per_height, Fz_per_height_se, v0_per_height, Fz_pooled] = ...
        compute_fz_intercepts(heights, d1, z_targets, v_min, tol)
% COMPUTE_FZ_INTERCEPTS  Per-height F(zi)/m = mean(a+g - v²/d1) at fixed depths.
%
%   Lightweight — no fitting, no figure. Feeds directly into plot_fz_vs_z.
%   d1 must be pre-computed from plot_ag_vs_v2 (Fig 5).
%
%   Inputs:
%       heights    - struct array from group_trials_by_height
%       d1         - shared inertial length scale [cm] from Fig 5
%       z_targets  - depth targets [cm], e.g. 1.3:0.10:2.2
%       v_min      - minimum velocity cutoff [cm/s], default 40
%       tol        - depth window half-width [cm], default 0.04
%
%   Outputs:
%       Fz_per_height    - [nH x nZ] per-height intercepts [cm s^-2]
%       Fz_per_height_se - [nH x nZ] standard errors
%       v0_per_height    - [nH x 1]  mean impact speed per height [cm/s]
%       Fz_pooled        - [nZ x 1]  pooled mean across all heights

    if nargin < 4, v_min = 20;   end
    if nargin < 5, tol   = 0.04; end

    z_targets = z_targets(:);
    nZ        = numel(z_targets);
    nH        = numel(heights);
    ag_cap = 30000;   % allows v up to ~350 cm/s

    Fz_per_height    = nan(nH, nZ);
    Fz_per_height_se = nan(nH, nZ);
    v0_per_height    = nan(nH, 1);

    for j = 1:nH
        v0_per_height(j) = heights(j).v0_mean;

        for zi = 1:nZ
            v_hj = []; ag_hj = [];

            for i = 1:heights(j).nTrials
                k   = heights(j).trials(i).kinematics;
                idx = k.impact_index:k.stopFrame;
                z_t = k.z_smooth(idx);
                v_t = k.v_smooth(idx);
                ag_t = k.a_plus_g(idx);

                in = abs(z_t - z_targets(zi)) < tol & ...
                v_t >= v_min & isfinite(v_t) & isfinite(ag_t) & ...
                v_t > 0 & ag_t < ag_cap;

                v_hj  = [v_hj;  v_t(in)];
                ag_hj = [ag_hj; ag_t(in)];
            end

            if numel(v_hj) >= 3
                fz_vals              = ag_hj - v_hj.^2 ./ d1;
                Fz_per_height(j,zi)   = mean(fz_vals);
                Fz_per_height_se(j,zi)= std(fz_vals) / sqrt(numel(fz_vals));
            end
        end
    end

    % Pooled mean per depth (mean of per-height means)
    Fz_pooled = nan(nZ, 1);
    for zi = 1:nZ
        vals = Fz_per_height(:,zi);
        ok   = isfinite(vals);
    if nnz(ok) >= 2
        Fz_pooled(zi) = mean(vals(ok));
    end
    end
end
