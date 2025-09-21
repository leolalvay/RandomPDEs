function [mu_star, f_star] = optimize_mu(K)
%OPTIMIZE_MU_SHIFT_DILATION  Find mu^* = argmax f(x) for the shift-dilation IS.
%   f(x) = max(exp(x1)+exp(x2)-K, 0) * exp(-(x1^2+x2^2)/2).
%   INPUT:  K > 0 (strike-like parameter)
%   OUTPUT: mu_star (2x1 vector), f_star = f(mu_star)
%
% Notes:
% - We use fminsearch (Nelder–Mead) which MINIMIZES a function, so we
%   minimize -f(x). This is derivative-free (no gradients needed).
% - Start point uses the heuristic x1=x2=log(K/2) which sits near the
%   ridge e^{x1}+e^{x2}≈K. It works well for K≥2; for smaller K it is still OK.

    if K <= 0
        error('K must be positive.');
    end

    % ---------- objective handle (returns negative f) ----------
    % x is a 2-vector [x1; x2]. We negate f for minimization.
    obj = @(x) -f_payoff_times_rho(x(1), x(2), K);

    % ---------- initial guess ----------
    t   = log(K/2);                 % symmetric guess near the active region
    x0  = [t; t];

    % ---------- basic options for readability ----------
    % 'Display' can be 'off'|'iter'; TolX/TolFun control stopping accuracy.
    opts = optimset('Display','off','TolX',1e-8,'TolFun',1e-10,'MaxIter',400);

    % ---------- Nelder–Mead search ----------
    [mu_star, neg_f_star] = fminsearch(obj, x0, opts);

    % Return the actual maximum value
    f_star = -neg_f_star;
end

% ======= helper: f(x) without the constant 1/(2*pi) (doesn't affect argmax) =======
function val = f_payoff_times_rho(x1, x2, K)
% Compute f(x) = max(exp(x1)+exp(x2)-K,0) * exp(-(x1^2+x^2)/2)
% All operations are scalar; this keeps the function easy to read.
    s = exp(x1) + exp(x2) - K;
    if s <= 0
        val = 0.0;
    else
        val = s * exp(-0.5*(x1^2 + x2^2));
    end
end
