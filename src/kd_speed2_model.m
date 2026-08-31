function u = kd_speed2_model(z, v0, k, d1, mass, grav)
%KD_SPEED2_MODEL  Exact v^2(z) solution of the Katsuragi-Durian force law.
%
%   u = kd_speed2_model(z, v0, k, d1, mass, grav)
%
%   The published KD law, in the sign convention this project verified,
%
%       grav - accel = k*depth/mass + speed^2/d1
%
%   is a linear ODE in u = v^2 once written against depth rather than time:
%
%       du/dz + (2/d1)*u = 2*grav - 2*k*z/mass
%
%   Its exact solution with the measured v0 as the boundary condition,
%   u(0) = v0^2, is what this function returns. Writing E = exp(-2*z/d1),
%
%       u(z) = v0^2*E + grav*d1*(1-E) - k*( (d1/mass)*z - (d1^2/(2*mass))*(1-E) )
%
%   Fitting this closed form rather than a finite-difference acceleration
%   keeps the fit on the measured quantity and avoids differentiating noise.
%
%   v0 may be scalar (one trial) or a vector the same size as z (points
%   pooled across trials, each carrying its own impact speed).
%
%   The derivative is kd_dspeed2_dz; scripts/step5_exemplar_fits.m checks the
%   two against each other analytically and numerically before fitting.
%
%   Base MATLAB only.

    Ez = exp(-2*z/d1);
    u = v0(:).^2 .* Ez + grav*d1*(1-Ez) - k*( (d1/mass)*z - (d1^2/(2*mass))*(1-Ez) );
end
