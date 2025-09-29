function [L, effective_nugget] = matern_cholesky_adaptive(x, nu, rho, sigma2, base_nugget, max_tries)
% Returns lower-triangular L with (Sigma + effective_nugget*I) = L*L'
% Strategy: try Cholesky WITHOUT nugget first. Only if it fails, add nugget.

    if nargin < 5 || isempty(base_nugget), base_nugget = []; end
    if nargin < 6 || isempty(max_tries),   max_tries   = 6;  end

    % Sigma already symmetrized inside matern_covariance
    Sigma = matern_covariance(x, nu, rho, sigma2);
    n = size(Sigma,1);

    % Try Cholesky WITHOUT nugget
    [L, p] = chol(Sigma, 'lower');
    if p == 0
        effective_nugget = 0.0;
        return;
    end

    % If it failed, define a reasonable base_nugget (scaled to diag level)
    if isempty(base_nugget)
        diag_scale  = mean(diag(Sigma));  % ~ sigma2
        base_nugget = max(1e-12*max(1,diag_scale), eps(diag_scale));
    end

    % Retry growing nugget by x10
    effective_nugget = 0.0;
    for k = 0:max_tries
        effective_nugget = (10^k) * base_nugget;
        [L, p] = chol(Sigma + effective_nugget*eye(n), 'lower');
        if p == 0
            return;
        end
    end

    error('Cholesky failed even after adding nugget up to %.3e. Try fewer points or larger rho.', effective_nugget);
end
