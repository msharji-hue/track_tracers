function [k_fit, d1_fit, rss, rmse_speed2] = kd_fit_nested(depth, speed2, v0, mass, grav, w)
%KD_FIT_NESTED  Nested least-squares fit of the KD v^2(z) law.
%
%   [k_fit, d1_fit, rss, rmse_speed2] = kd_fit_nested(depth, speed2, v0, mass, grav)
%   [...] = kd_fit_nested(depth, speed2, v0, mass, grav, w)
%
%   kd_speed2_model is LINEAR in k at fixed d1, so the fit nests:
%
%       inner   closed-form (weighted) least squares for k given d1
%       outer   fminbnd over d1 in [0.1, 30] cm minimising the residual sum
%
%   Only d1 is actually searched, which is why a single trial fits in
%   milliseconds and a 1000-sample cluster bootstrap is affordable.
%
%   k is REPORTED, never clamped: a negative k at the optimum is a fact about
%   the trial and the caller flags it rather than this function hiding it.
%
%   INPUTS are already masked -- masking is the caller's decision, not this
%   function's. v0 is scalar (one trial) or per-point (pooled across trials).
%   w is an optional per-point weight vector (default 1). Weighting is how the
%   trial-balanced global estimator gives every trial equal say regardless of
%   how many frames it contributed; with w omitted the arithmetic reduces
%   exactly to the unweighted case.
%
%   Base MATLAB only.

    if nargin < 6 || isempty(w), w = ones(size(depth)); end
    depth = depth(:); speed2 = speed2(:); w = w(:);

    % outer: one-dimensional search over the drag length
    obj = @(d1) local_rss(d1, depth, speed2, v0, mass, grav, w);
    d1_fit = fminbnd(obj, 0.1, 30);

    % inner: the k that minimises the weighted RSS at the winning d1
    [k_fit, rss] = kd_inner_k(d1_fit, depth, speed2, v0, mass, grav, w);

    % weighted root-mean-square residual; equals sqrt(rss/n) when w is all ones
    rmse_speed2 = sqrt(rss / sum(w));
end

function rss = local_rss(d1, z, u, v0, mass, grav, w)
% Objective for the outer search: weighted RSS after the inner k solve.
    [~, rss] = kd_inner_k(d1, z, u, v0, mass, grav, w);
end
