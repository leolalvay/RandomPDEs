%% Visualization of f(x1,x2) = max(exp(x1)+exp(x2)-K,0) * exp(-(x1^2+x2^2)/2)
% This script plots both a surface and a filled contour for two K values.
% Key MATLAB functions used:
% - meshgrid: builds 2D coordinate grids from 1D vectors (for vectorized evals).
% - surf: 3D surface plot.
% - contourf: filled contour (2D heatmap-like) plot.
% - tiledlayout/nexttile: modern layout manager for multiple axes.
% - colorbar: shows the color scale used in the plot.

clear; clc;

Ks = [3, 6];              % You can change or add values here
xlim_plot = [-5, 5];      % Domain for x1, x2 (wide enough to see the action)
gridN = 501;              % Grid resolution (odd -> nice center point)

x = linspace(xlim_plot(1), xlim_plot(2), gridN);
[X1, X2] = meshgrid(x, x);   % Create 2D grid (X1(i,j)=x_j, X2(i,j)=x_i)

% Define f(x1,x2;K) as a vectorized anonymous function.
f = @(X1,X2,K) max(exp(X1) + exp(X2) - K, 0) .* exp(-0.5*(X1.^2 + X2.^2));
% (All operations are elementwise: .*, .^, etc., so it works on matrices.)

tiledlayout(numel(Ks), 2, "TileSpacing", "compact", "Padding", "compact");

for idx = 1:numel(Ks)
    K = Ks(idx);
    F = (1/(2*pi))*f(X1, X2, K);

    % ---- Surface plot ----
    nexttile;
    surf(X1, X2, F, 'EdgeColor', 'none');  % 'EdgeColor','none' removes mesh lines
    view(40, 30);                          % 3D viewing angle (azimuth, elevation)
    camlight headlight; lighting gouraud;  % Simple lighting to see curvature
    xlabel('x_1'); ylabel('x_2'); zlabel('f(x_1,x_2)');
    title(sprintf('f(x_1,x_2), K = %g', K));
    colorbar; grid on; box on;

    % ---- Filled contour + positivity boundary ----
    nexttile;
    contourf(X1, X2, F, 30, 'LineColor', 'none');  % 30 levels, no contour lines
    hold on;
    % Positivity boundary e^{x1}+e^{x2}=K -> x2 = log(K - e^{x1}) for x1 < log K
    x1b = linspace(xlim_plot(1), min(log(K)-1e-8, xlim_plot(2)), 800);
    mask = (exp(x1b) < K);                 % ensure argument of log is positive
    x1b = x1b(mask);
    x2b = log(K - exp(x1b));
    % Only keep boundary points inside the visible window:
    inWin = x2b >= xlim_plot(1) & x2b <= xlim_plot(2);
    plot(x1b(inWin), x2b(inWin), 'k--', 'LineWidth', 1.5);  % dashed black curve
    hold off;

    axis equal tight;                       % equal scales for x and y; fit tightly
    xlabel('x_1'); ylabel('x_2');
    title(sprintf('f(x_1,x_2), K = %g', K));
    colorbar; grid on; box on;
end

% Or just one axes (one tile)
% Whole figure
exportgraphics(gcf, 'plot.pdf', 'ContentType', 'image', 'Resolution', 200);

% Notes:
% - For very large |x|, exp(x) may overflow. On this domain ([-4,4]) it is safe.
% - If you want the normalized density, multiply F by (1/(2*pi)).
% - You can switch to semilogy-like visualization by plotting log(F+eps).
