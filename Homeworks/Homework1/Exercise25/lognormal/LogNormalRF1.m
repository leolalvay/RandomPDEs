function [a, kappa, eff_nugget] = LogNormalRF1(x, nu, rho, sigma2, seed)
% Minimal sampler: one log-normal field realization on points x.
% Uses your matern_cholesky_adaptive (and matern_covariance inside it).
%
% Inputs:y si no quiero ponerle un seed com
%   x      : points in [0,1] (column or row)
%   nu     : Matern smoothness (e.g., 0.5, 1.5, 2.5, Inf)
%   rho    : correlation length
%   sigma2 : pointwise variance of kappa
%   seed   : (optional) RNG seed
%
% Outputs:
%   a         : exp(kappa), log-normal field at x
%   kappa     : Gaussian field at x (mean zero, Matern covariance)
%   eff_nugget: nugget actually used by Cholesky (0 if none)

    if nargin >= 5 && ~isempty(seed), rng(seed); end
    [L, eff_nugget] = matern_cholesky_adaptive(x, nu, rho, sigma2);
    z     = randn(numel(x),1);
    kappa = L * z;
    a     = exp(kappa);
end
