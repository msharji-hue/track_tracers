function [z_smooth, v_smooth, a_smooth] = poly_kinematics(t, z, degree)
% POLY_KINEMATICS  Fit z(t) to a polynomial and return smooth z, v, a
%                  via analytical differentiation — no noise amplification.
%
%   Inputs:
%       t      - time vector (s)
%       z      - depth vector (cm)
%       degree - polynomial degree (default 4)
%
%   Outputs:
%       z_smooth - smoothed depth (cm)
%       v_smooth - velocity dz/dt (cm/s)
%       a_smooth - acceleration d²z/dt² (cm/s²)

    if nargin < 3, degree = 3; end

    t = t(:); z = z(:);
    ok  = isfinite(t) & isfinite(z);
    p   = polyfit(t(ok), z(ok), degree);
    dp  = polyder(p);
    ddp = polyder(dp);

    z_smooth = polyval(p,   t);
    v_smooth = polyval(dp,  t);
    a_smooth = polyval(ddp, t);
end
