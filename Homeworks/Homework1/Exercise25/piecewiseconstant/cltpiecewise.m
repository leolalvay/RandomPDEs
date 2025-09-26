clear; clc;

M=1e4;
sigma = 0.15;
N = [10 20 40];
l = 2;

% --- CLT 95% band (global sigma) ---
alpha  = 0.05;
zalpha = -sqrt(2)*erfcinv(2*(1 - alpha/2));   % ≈ 1.96

step = 10;
Ms   = step:step:M;      % [1000 2000 ... 10000]
Mn = numel(Ms);
m    = Ms(:);                  % columna
n = numel(N);
C = zeros(Mn,n);

for i =1:n
Q = sampleQuPW(M, sigma, N(i), l);

c1 = cumsum(Q);          % sum_{i=1}^m Q_i
c2 = cumsum(Q.^2);       % sum_{i=1}^m Q_i^2

mean_m   = c1(m)./m;               % media en cada m
% Varianza muestral insesgada en cada m (m>=2):
s2_m     = (c2(m) - m.*(mean_m.^2))./(m - 1); 
C(:,i)  = zalpha *sqrt(s2_m) ./ sqrt(m);
end

% --- Plot ---
figure('Color','w','Units','inches','Position',[1 1 9.5 4.0]);
loglog(m, C(:,1), 'LineWidth',1.5); hold on;  
loglog(m, C(:,2), 'LineWidth',1.5);            
loglog(m, C(:,3), 'LineWidth',1.5);            
grid on;
xlabel('M','Interpreter','latex','FontSize',14);
ylabel('Statistical error', 'Interpreter','latex', 'FontSize',14);

legend( ...
    sprintf('$N = %d$', N(1)), ...
    sprintf('$N= %d$', N(2)), ...
     sprintf('$N= %d$', N(3)), ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
title(sprintf('CLT bounds for $Q(u_h)$ (M=%d, $\\sigma$=%.2f)', M, sigma), ...
      'Interpreter','latex','FontSize',14);

exportgraphics(figure(1), 'cltpiecewise.pdf', 'ContentType', 'vector');