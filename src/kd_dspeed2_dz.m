function du = kd_dspeed2_dz(z, v0, k, d1, mass, grav)
%KD_DSPEED2_DZ  Closed-form d/dz of kd_speed2_model.
%
%   du = kd_dspeed2_dz(z, v0, k, d1, mass, grav)
%
%   Differentiating kd_speed2_model term by term, with E = exp(-2*z/d1) and
%   dE/dz = -(2/d1)*E:
%
%       d/dz [ v0^2*E ]                  = -(2/d1)*v0^2*E
%       d/dz [ grav*d1*(1-E) ]           = +2*grav*E
%       d/dz [ -k*(d1/mass)*z ]          = -k*(d1/mass)
%       d/dz [ k*(d1^2/(2*mass))*(1-E) ] = +k*(d1/mass)*E
%
%   which collects to
%
%       du/dz = -2/d1 * [ v0^2*E - grav*d1*E + k*(d1^2/(2*mass))*(1-E) ]
%
%   and is the same formula as the ODE right-hand side
%   2*grav - 2*u/d1 - 2*k*z/mass, rearranged. scripts/step5_exemplar_fits.m
%   verifies that identity to 1e-9 relative before any fitting runs.
%
%   Base MATLAB only.

    Ez = exp(-2*z/d1);
    du = -2/d1 * ( v0(:).^2 .* Ez - grav*d1*Ez + k*(d1^2/(2*mass))*(1-Ez) );
end
