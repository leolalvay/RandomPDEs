clear; clc;
%This script plots realizations of the piecewise constant random fields
rng(202);
sigma = 0.25;
Ny = 20;
Nw = 50;
Yy = unifrnd(-sqrt(3), sqrt(3), [Ny 1]);
Yw = unifrnd(-sqrt(3), sqrt(3), [Nw 1]);
x = linspace(0,1,1000);
y = PiecewiseRF(x, sigma, Yy);
w = PiecewiseRF(x, sigma, Yw);
figure(1);
set(1,'position',[400 200 800 400]);
plot(x, y, 'b', x, w, 'r', 'LineWidth', 2)
xlabel('$x\in(0,1)$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$a(x,\omega)$', 'Interpreter', 'latex', 'FontSize', 16);
title( ...
    sprintf('Model 1 with $\\sigma=%.2f$ and $N=%d,\\; N=%d$', sigma, Ny, Nw), ...
    'Interpreter','latex','FontSize',16);
grid on;
legend( ...
    sprintf('$N = %d$', Ny), ...
    sprintf('$N= %d$', Nw), ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);
% export figure 1 in PDF
exportgraphics(figure(1), 'model1_plot.pdf', 'ContentType', 'vector');
