function d_pred = run_forward_model_simple(v0, d1, k_over_m, dz)
% RUN_FORWARD_MODEL_SIMPLE  Integrate Katsuragi force law to predict
%   final penetration depth.
%
%   F/m = -g + (k/m)*z + v²/d1
%
%   No area scaling — direct analogue of Katsuragi Fig. 2b dotted curve.
%   Parameters d1 and k/m come from kinematic fits (Figs 3a and 3b).
%
%   Inputs:
%       v0        - impact speed [cm/s]
%       d1        - inertial length scale [cm]
%       k_over_m  - friction coefficient [s^-2]
%       dz        - depth step size [cm], default 0.001
%
%   Output:
%       d_pred    - predicted final penetration depth [cm]

    if nargin < 4, dz = 0.001; end

    g  = 980;   % cm/s²
    z  = 0;
    v2 = v0^2;

    while v2 > 0
        F_over_m = g - k_over_m * z - v2 / d1;
        v2_new   = v2 + 2 * F_over_m * dz;
        if v2_new <= 0, break; end
        v2 = v2_new;
        z  = z + dz;
    end

    d_pred = z;
end
