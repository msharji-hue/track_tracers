function [d_pred, z_traj, v_traj] = run_forward_model(v0, d1, k_over_m, ...
                                                        A_eff_cm2, depth_cm, ...
                                                        A_ref_cm2, dz_cm)
% RUN_FORWARD_MODEL  Numerically integrates the proposed force law:
%
%   F/m = -g + (A_eff(z)/A_ref) * [k/m * z  +  v²/d1]
%
%   This is a GeoPoncelet + Katsuragi hybrid:
%     - Gravity term:   -g             (restored from GeoPoncelet)
%     - Friction term:   k/m * |z|     (Katsuragi depth-linear friction)
%     - Inertial term:   v²/d1         (Katsuragi inertial drag)
%     - Area scaling:    A_eff(z)/A_ref (GeoPoncelet depth-varying geometry)
%
%   Integration uses v² form for numerical stability:
%     v²_{i+1} = v²_i + 2*dz * F_net/m
%
%   Inputs:
%     v0         - impact speed (cm/s)
%     d1         - inertial drag length scale (cm), from plot_ag_vs_v2
%     k_over_m   - friction coefficient (s^-2), from plot_fz_vs_z
%     A_eff_cm2  - effective area profile (cm²), smoothed, from STL extraction
%     depth_cm   - depth vector corresponding to A_eff_cm2 (cm)
%     A_ref_cm2  - reference area at which d1 and k/m were fitted (cm²)
%     dz_cm      - integration step size in cm (default 0.005)
%
%   Outputs:
%     d_pred  - predicted final penetration depth (cm)
%     z_traj  - depth trajectory (cm)
%     v_traj  - velocity trajectory (cm/s)

    if nargin < 7, dz_cm = 0.005; end

    g_cm = 980;   % cm/s²

    % Interpolant for A_eff(z) — clamp extrapolation to endpoint values
    A_interp = @(z) interp1(depth_cm, A_eff_cm2, z, 'pchip', ...
                            A_eff_cm2(end));

    % Initialise
    z  = 0;
    v2 = v0^2;

    z_traj  = z;
    v_traj  = v0;

    max_depth = max(depth_cm) * 0.98;   % safety cap at STL extent

    while v2 > 0 && z < max_depth

        A_z      = max(A_interp(z), 0);
        area_ratio = A_z / A_ref_cm2;

        % Net acceleration (cm/s²)
        F_over_m = -g_cm + area_ratio * (k_over_m * z  +  v2 / d1);

        % v² update
        v2_new = v2 + 2 * dz_cm * F_over_m;

        if v2_new <= 0
            % Interpolate to find exact stopping depth
            frac   = v2 / (v2 - v2_new);
            z_stop = z + frac * dz_cm;
            z_traj(end+1) = z_stop;
            v_traj(end+1) = 0;
            z = z_stop;
            break;
        end

        z  = z  + dz_cm;
        v2 = v2_new;

        z_traj(end+1) = z;
        v_traj(end+1) = sqrt(v2);
    end

    d_pred = z_traj(end);
end
