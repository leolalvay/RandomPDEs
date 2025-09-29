function Sigma = matern_covariance(x, nu, rho, sigma2)
% x      : vector of points in [0,1] (row or column)
% nu     : smoothness (0.5, 1.5, 2.5, Inf)
% rho    : correlation length (>0)
% sigma2 : pointwise variance of kappa (>=0)
%
% Sigma(i,j) = C(|x_i - x_j|)

x = x(:);                        % ensure column
r = abs(x - x.');                % pairwise distances (n x n)
switch true
    case nu == 0.5
        Sigma = sigma2 * exp( - r / rho );

    case nu == 1.5
        t = sqrt(3) * r / rho;
        Sigma = sigma2 * (1 + t) .* exp(-t);

    case nu == 2.5
        % Correct Matérn(ν=2.5): (1 + sqrt(5) r/rho + 5 r^2/(3 rho^2)) exp(-sqrt(5) r/rho)
        t = sqrt(5) * r / rho;
        Sigma = sigma2 * (1 + t + (t.^2) * (1/3)) .* exp(-t);   % since t^2 = 5 r^2 / rho^2

    case isinf(nu)   % squared-exponential (Gaussian kernel) as ν→∞
        Sigma = sigma2 * exp( - (r.^2) / (2 * rho^2) );

    otherwise
        error('This helper implements nu in {0.5, 1.5, 2.5, Inf}.');
end
% numerical symmetrization (defensive). To avoid possible non symmetry of
% the Covariance matrix due to rounding operations.
Sigma = 0.5 * (Sigma + Sigma.');
end
