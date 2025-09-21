clear; clc;

K = 3;
mu = [1.3224; 0.7348];                 % your chosen peak
Sigma = [1.3, -0.6; -0.6, 1.3];  % your visual-tuned covariance
N = 2e5;                         % one long run (adjust as you like)
alpha = 0.05;
m = (1:N).';

res = running_convergence_mc_is(K, mu, Sigma, N, alpha);


function res = running_convergence_mc_is(K, mu, Sigma, N, alpha, seed)
%RUNNING_CONVERGENCE_MC_IS  Smooth convergence curves using cumulative stats.
%   res = running_convergence_mc_is(K, mu, Sigma, N, alpha, seed)
%
% INPUTS
%   K     : >0
%   mu    : 2x1 mean for IS proposal N(mu,Sigma)
%   Sigma : 2x2 PD covariance for IS proposal
%   N     : total number of samples in ONE long run (e.g., 2e5 or 1e6)
%   alpha : significance (default 0.05 -> 95% CI)
%   seed  : RNG seed (optional)
%
% OUTPUTS (all length-N column vectors)
%   res.M            : 1:N
%   res.Ihat_MC      : running mean for MC
%   res.HW_MC        : running CLT half-width for MC
%   res.Ihat_IS      : running mean for IS
%   res.HW_IS        : running CLT half-width for IS
%
% Notes:
% - We compute unbiased sample variance via cumulative sums:
%     s^2_m = (sum(x^2) - m*mean_m^2)/(m-1), m>=2.
%   Then half-width = z * sqrt(s^2_m / m).
% - z for two-sided (1-alpha) CI is z = sqrt(2)*erfcinv(alpha),
%   which avoids toolbox calls (erfcinv is base MATLAB).
% - IS weights are computed stably via Cholesky; constants 1/(2π) cancel.

    if nargin < 5 || isempty(alpha), alpha = 0.05; end
    if nargin >= 6 && ~isempty(seed), rng(seed); end

    

    % z_{1-alpha/2} without Statistics Toolbox:
    % norminv(1 - a/2) = sqrt(2)*erfcinv(a)
    z = sqrt(2) * erfcinv(alpha);

    % ----- One long run: MC stream -----
    X = randn(N,2);                                        % X ~ N(0,I)
    gX = max(exp(X(:,1)) + exp(X(:,2)) - K, 0);            % payoff
    c1 = cumsum(gX);
    c2 = cumsum(gX.^2);
    m  = (1:N).';
    mean_MC = c1 ./ m;
    var_MC  = nan(N,1);
    %unbiases variance estimator (1/M-1 * mean of squares)
    var_MC(2:end) = (c2(2:end) - m(2:end).*mean_MC(2:end).^2) ./ (m(2:end)-1);
    HW_MC   = z .* sqrt(var_MC ./ m);

    % --- IS stream (Y ~ N(mu,Sigma))
    L   = chol((Sigma+Sigma.')/2, 'lower');%Sigma = LL^T, L(2x2)
    Z   = randn(N,2);
    Y   = mu.' + (Z * L.'); % (Nx2)Proposition 7, Stochastic Notes (Transpose of eq (16))
    gY  = max(exp(Y(:,1)) + exp(Y(:,2)) - K, 0);

    Ym  = Y - mu.';                          % Y(i,:)=[z1-mu1,z2-mu2]_i
    V    = (L \ Ym.').';                     % solve L*V' = Ym'
    qf   = sum(V.^2,2);                      % (y-mu)'Sigma^{-1}(y-mu)
    t1   = sum(Y.^2,2);                      % ||y||^2
    sqrt_detSigma = prod(diag(L));           % |Sigma|^{1/2}
    w    = sqrt_detSigma .* exp(-0.5*t1 + 0.5*qf);

    Zis  = gY .* w;                          % weighted terms
    c1w  = cumsum(Zis);                      % sum Z
    c2w  = cumsum(Zis.^2);                   % sum Z^2
    mean_IS = c1w ./ m;
    var_IS  = nan(N,1);
    var_IS(2:end) = (c2w(2:end) - m(2:end).*mean_IS(2:end).^2) ./ (m(2:end)-1);
    HW_IS = z .* sqrt(var_IS ./ m);



    % ----- pack results -----
    res = struct('M', m, ...
                 'Ihat_MC', mean_MC, 'HW_MC', HW_MC, ...
                 'Ihat_IS', mean_IS, 'HW_IS', HW_IS);

    % ----- plot: smooth convergence curves (downsample for readability) -----
    % Select ~300 log-spaced indices to avoid plotting N points.
    %k = min(300, N);
    %idx = unique(round(logspace(log10(10), log10(N), k)));  % start at 10 for stable var
    fig = figure('Color','w');
    loglog(m(10:end), HW_MC(10:end), '-', 'LineWidth',1.4); hold on;
    loglog(m(10:end), HW_IS(10:end), '-', 'LineWidth',1.4);
    grid on; grid minor;
    xlabel('M');
    ylabel(sprintf('Error'));
    xlim([1e1, N]);
    title(sprintf('Convergence Rate (K = %g)', K));
    legend('Standard MC','Importance Sampling','Location','southwest');
    exportgraphics(fig, sprintf('convergence_K%d.pdf', K), 'ContentType','image', 'Resolution', 200);

end

