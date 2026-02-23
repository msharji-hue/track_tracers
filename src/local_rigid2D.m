function [R, tvec] = local_rigid2D(X0, X)
% Solve X ≈ R*X0 + t with least squares, R orthonormal.
% X0, X are Nx2 (same IDs). Returns R (2x2) and tvec (2x1).

c0 = mean(X0,1);
c  = mean(X,1);

A = (X0 - c0).';   % 2xN
B = (X  - c ).';   % 2xN

H = B * A.';       % 2x2
[U,~,V] = svd(H);

R = U*V.';
if det(R) < 0
    U(:,2) = -U(:,2);
    R = U*V.';
end

tvec = c.' - R*c0.';
end

