clear; clc;

M = 1e4;
N= 10;
ell = 2;
sigma = 0.15;

%============= Estimation of phat====================
% 1) Generate once the same Y for all levels
    Y = unifrnd(-sqrt(3), sqrt(3), [M, N]);

    % 2) Monte Carlo means at each level (paired by Y)
    Qh  = sampleQuPWgivenY(sigma, ell,   Y);
    Qh2 = sampleQuPWgivenY(sigma, ell+1, Y);
    Qh4 = sampleQuPWgivenY(sigma, ell+2, Y);

    mu_h  = mean(Qh);
    mu_h2 = mean(Qh2);
    mu_h4 = mean(Qh4);

    % 3) Order estimate (micro-paso 3)
    Delta1 = mu_h  - mu_h2;
    Delta2 = mu_h2 - mu_h4;
    p_hat  = log2(abs(Delta1) / abs(Delta2));

%========================================================    

% ---- Richardson (mu as h->0) ----
mu_Rich = mu_h2 + (mu_h2 - mu_h) / (2^p_hat - 1);


