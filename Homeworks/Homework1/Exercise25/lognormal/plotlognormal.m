clear; clc; close all;

%% 0) User parameters (edit freely)
I       = 500;          % number of elements (midpoints will be used)
%nu      = 0.5;         % smoothness: 0.5, 1.5, 2.5, Inf
rho     = 0.1;         % correlation length
sigma2  = 2.0;         % pointwise variance of kappa
%seed    = 10;        % RNG seed; set [] to skip seeding

%% 1) Build points in [0,1] (element midpoints for piecewise-constant a)
h  = 1/I;
x  = ((1:I) - 0.5).' * h;   % column vector of midpoints
n  = numel(x);

%% 2) RNG seed (optional)
%if ~isempty(seed), rng(seed); end

%% 3) Plots
y = LogNormalRF1(x, 0.5,rho,sigma2);
z = LogNormalRF1(x, 1.5,rho,sigma2);
v = LogNormalRF1(x, 2.5,rho,sigma2);
w = LogNormalRF1(x, Inf,rho,sigma2);

figure(1);
set(1,'position',[400 200 800 400]);
plot(x,y, x,z,x,v,x,w, 'LineWidth', 2)
xlabel('$x\in(0,1)$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$a(x,\omega)$', 'Interpreter', 'latex', 'FontSize', 16);
title('Realizations of Model 2','Interpreter', 'latex', 'FontSize', 16); grid on;
legend('$\nu = 0.5$', '$\nu = 1.5$','$\nu=2.5$','$\nu\to\infty$', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);

exportgraphics(gcf,'Realizations_Model2.pdf','ContentType','vector');
