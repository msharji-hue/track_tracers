function [k, rss] = kd_inner_k(d1, z, u, v0, mass, grav, w)
%KD_INNER_K  Closed-form (weighted) least squares for k at fixed d1.
%
%   [k, rss] = kd_inner_k(d1, z, u, v0, mass, grav, w)
%
%   At fixed d1 the KD solution splits into a k-independent part and a part
%   proportional to k,
%
%       u(z) = A(z) + k*B(z)
%       A = v0^2*E + grav*d1*(1-E)
%       B = -( (d1/mass)*z - (d1^2/(2*mass))*(1-E) )
%
%   so k drops out of a one-line normal equation and never needs iterating.
%   This is the inner half of kd_fit_nested; it is separate only because
%   MATLAB puts one public function per file.
%
%   k is returned as computed. It is NOT clamped at zero -- a negative k means
%   the trial's speed record wants a force that grows toward the surface, and
%   the caller is expected to flag that rather than have it silently absorbed.
%
%   w is an optional per-point weight vector (default 1).
%
%   Base MATLAB only.

    if nargin < 7 || isempty(w), w = ones(size(z)); end
    z = z(:); u = u(:); w = w(:);

    Ez = exp(-2*z/d1);
    A = v0(:).^2 .* Ez + grav*d1*(1-Ez);             % k-independent part
    B = -( (d1/mass)*z - (d1^2/(2*mass))*(1-Ez) );   % coefficient of k

    k = sum(w .* B .* (u - A)) / sum(w .* B.^2);     % reported, never clamped
    rss = sum(w .* (u - A - k*B).^2);
end
