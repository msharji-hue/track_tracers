function ag = net_accel(kin, calib)
%NET_ACCEL  Net grain acceleration, a + g = g - a, from the RAW a trace.
%
%   ag = net_accel(kin)          % g taken as 980 cm/s^2
%   ag = net_accel(kin, calib)   % g taken from calib.g_cm_s2
%
%   RECOMPUTES rather than reading kin.a_plus_g, and that is the point.
%   _kin.mat files written before the 2026-08 correction carry the old formula
%   (-a - g) in their a_plus_g column: the same quantity offset by a constant
%   -2g with the rest state inverted. The raw a trace those files store is
%   unaffected, so deriving a+g from it here makes every already-processed
%   trial render correctly with no Stage B re-run.
%
%   THE CONVENTION. Depth is positive into the bed, so a free-falling rod has
%   a = +g and the acceleration the grains supply is what remains:
%
%       a + g = g - a
%
%       free fall      a = +g    ->  0             nothing but gravity
%       at rest        a =  0    ->  +g            the bed carries the weight
%       decelerating   a = -|a|  ->  g + |a|       grains resisting, large +ve
%
%   DISPLAY QUANTITY ONLY. No fit reads it: a_stop is fitted to velocity in
%   kd_kinematics/find_stop, v0 comes from v_smooth, d_final from depthRod_cm.
%   Fixing it therefore changes figures and nothing else.
%
%   Inherits kd_kinematics' masking: a is NaN outside [impact_index,
%   stopFrame], so a+g is NaN there too.

    if nargin < 2 || ~isstruct(calib) || ~isfield(calib, 'g_cm_s2')
        g = 980;
    else
        g = calib.g_cm_s2;
    end

    if ~isstruct(kin) || ~isfield(kin, 'a')
        error('net_accel:noRawAccel', ...
            ['kin has no raw ''a'' trace, so a+g cannot be recomputed. Very ' ...
             'old _kin.mat files may predate it; re-run Stage B for that ' ...
             'trial rather than trusting the stored a_plus_g column, which ' ...
             'carries the pre-2026-08 formula.']);
    end

    ag = g - kin.a(:);
end
