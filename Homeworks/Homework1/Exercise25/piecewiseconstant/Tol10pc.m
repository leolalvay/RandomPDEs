clear; clc;

% -------- Parameters you can tune --------
N         = 10;        % number of random blocks in the coefficient
sigma     = 0.15;      % ensure 1 - sigma*sqrt(3) > 0
ell0      = 1;         % starting mesh level (h = 1/(N*2^ell))
M0        = 10;      % initial number of MC samples
tau       = 0.007;      % target relative error (e.g., 10%)
splitBias = 0.7;      % fraction of tau for bias (rest for statistical)

% -------- Run adaptive selection --------
res = choose_hM_piecewise(N, sigma, ell0, M0, tau, splitBias);


% -------- Minimal report --------
h = 1/(N*2^res.ell);
fprintf('Chosen ell = %d  (h = %.6g)\n', res.ell, h);
fprintf('Chosen M   = %d\n', res.M);
fprintf('p_hat      = %.4f\n', res.p_hat);
fprintf('mu_Rich    = %.6e\n', res.mu_Rich);
fprintf('rel bias   = %.3g\n', res.rel_bias);
fprintf('rel stat   = %.3g\n', res.rel_stat);
fprintf('total rel  = %.3g\n', res.rel_bias + res.rel_stat);
fprintf('taub       = %.4f\n', tau*splitBias);
fprintf('taus       = %.4f\n', tau*(1-splitBias))




function out = choose_hM_piecewise(N, sigma, ell0, M0, tau, splitBias)
% Adaptive selection of ell (h=1/(N*2^ell)) and M to meet relative error < tau
% splitBias in (0,1): fraction of tau for bias; rest for statistical error.

    if nargin < 6, splitBias = 0.5; end
    tau_b = tau*splitBias; 
    tau_s = tau*(1-splitBias);
    z = 1.96;                          % 95% CL
    ell = ell0; 
    M   = M0;
    maxIters = 12;

    for it = 1:maxIters
        % --- common randomness for this iteration ---
        rng(202 + it); 
         Y = unifrnd(-sqrt(3), sqrt(3), [M, N]);

        % --- paired MC at three levels (ell, ell+1, ell+2) ---
        Qh  = sampleQuPWgivenY(sigma,ell,   Y);
        Qh2 = sampleQuPWgivenY(sigma,ell+1, Y);
        Qh4 = sampleQuPWgivenY(sigma,ell+2, Y);

        mu_h  = mean(Qh);  mu_h2 = mean(Qh2);  mu_h4 = mean(Qh4);
        Delta1 = mu_h  - mu_h2;
        Delta2 = mu_h2 - mu_h4;
        p_hat  = log2(abs(Delta1)/max(eps,abs(Delta2)));

        % Richardson reference for relative error
        mu_Rich = mu_h2 + (mu_h2 - mu_h) / (2^p_hat - 1);
        D = abs(mu_Rich); 
        if D < 1e-12, D = abs(mu_h2); end   % fallback if near 0

        % Bias component
        bias_h = abs(mu_h - mu_Rich);
        Rb = bias_h / D;

        % Statistical component at level h
        s_h = std(Qh);
        SE  = z * s_h / sqrt(M);
        Rs  = SE / D;

        % Stopping condition
        if (Rb <= tau_b) && (Rs <= tau_s)
            break;
        end

        % Update (prioritize the one that fails)
        if Rb > tau_b
            k = ceil( log2(bias_h/(tau_b*D)) / max(p_hat,1e-6) );
            k = max(k,1);
            ell = ell + k;
            % keep M; re-estimate next iter
        elseif Rs > tau_s
            M_needed = ceil( (z*s_h/(tau_s*D))^2 );
            M = max(M, M_needed);
        end
    end

    out = struct('ell',ell,'M',M,'p_hat',p_hat, ...
                 'mu_h',mu_h,'mu_h2',mu_h2,'mu_h4',mu_h4, ...
                 'mu_Rich',mu_Rich,'bias_abs',bias_h, ...
                 'rel_bias',Rb,'rel_stat',Rs, ...
                 'Delta_check',abs(Delta1)/max(eps,abs(Delta2)));
end
