clc; %clear command window
clear; % Clear worspace
M=1e6;
N=20;
w1=0.5;
u = rand(M,N);
cN=9/N;
Iexact = osc_exact_equal(N,w1);
s = cN * sum(u, 2);                 % Mx1: (9/N)*row-sum
f = cos(2*pi*w1 + s);               % Mx1: integrand at each sample
run_aver = cumsum(f)./((1:M)');
%Ihat   = run_aver(end)


% --- Global sigma (use the same for CLT and B–E) ---
m = (1:M).';

% Prefer MLE for consistency with lambda3's mean-based estimator:
sigma = std(f, 1);   % divide by M (MLE)

% --- CLT 95% band (global sigma) ---
alpha  = 0.05;
zalpha = -sqrt(2)*erfcinv(2*(1 - alpha/2));   % ≈ 1.96
clt_m  = zalpha * sigma ./ sqrt(m);

% --- Berry–Esseen (95%) band (global lambda and same sigma) ---
mu_all    = mean(f);
m3abs_all = mean(abs(f - mu_all).^3);
lambda3   = m3abs_all / max(sigma^3, realmin);   % λ^3

CBE   = 30.51175;
phi_z = exp(-0.5*zalpha^2)/sqrt(2*pi);

Km    = (CBE * lambda3) ./ sqrt(m);                  % penalty term
C0_BE = zalpha + Km ./ (2*phi_z*(1+zalpha)^3);       % 1st-order correction
be_m  = C0_BE .* sigma ./ sqrt(m);                   % USE sigma (not sigma_m!)     % 

% --- Plot from 10^1 --------------------------------
err_m = abs(run_aver - Iexact);

figure;
loglog(m(10:end), max(err_m(10:end), eps), 'LineWidth', 1); hold on;
loglog(m(10:end), max(clt_m(10:end),  eps), 'LineWidth', 1);
loglog(m(10:end), max(be_m(10:end),   eps), 'LineWidth', 1);
xlim([1e1, M]);  grid on; xlabel('M'); ylabel('Statistical error');
legend('$|\epsilon_M|$', 'CLT', 'BE', ...
       'Interpreter', 'latex','Location','best');
title(sprintf('Oscillatory, N = %d', N));



%Exact Solution
function I = osc_exact_equal(N, w1)
cN   = 9/N;
theta = 2*pi*w1 + 9/2;                        % since sum c_n = 9
I     = cos(theta) * (sin(0.5*cN)/(0.5*cN))^N; % product collapses to a power
%I = cos(theta) * (sinc(cN/(2*pi)))^N;
end