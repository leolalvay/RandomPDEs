%% Level curves of f(x) with the base density N(0, I2)
% Problem: X ~ N(0, I2), g(x) = max(exp(x1)+exp(x2)-K, 0),
% f(x) = g(x) * rho(x), where rho is the standard bivariate normal pdf.
% We overlay contours of N(0, I2) to visualize the original density.

clear; clc;

Ks     = [3, 8];                  % K values to show
gridN  = 301;                     % grid resolution (bigger = smoother, slower)
x1span = linspace(-5, 8, gridN);  % x1 range (covers the interesting region)
x2span = linspace(-5, 8, gridN);  % x2 range
[X1, X2] = meshgrid(x1span, x2span);

% --- base Gaussian N(0, I2) (normalized; constant doesn't matter for contours) ---
rho = (1/(2*pi)) * exp(-0.5*(X1.^2 + X2.^2));   % standard bivariate normal

% Figure layout
fig = figure('Color','w');
tl  = tiledlayout(fig, 1, numel(Ks), 'TileSpacing','compact','Padding','compact');

for j = 1:numel(Ks)
    K = Ks(j);

    % ---- f(x) = max(exp(x1)+exp(x2)-K, 0) * exp(-0.5*(x1^2+x2^2)) ----
    % (We drop the constant 1/(2*pi) in f since it only rescales.)
    F = max(exp(X1) + exp(X2) - K, 0) .* exp(-0.5*(X1.^2 + X2.^2));

    % ---- plot: filled contours of f + white contours of N(0,I2) ----
    ax = nexttile(tl);
    contourf(ax, X1, X2, F, 20, 'LineColor','none'); hold(ax,'on');  % 20 levels
    levels = max(rho(:)) * linspace(0.1, 0.9, 6);                    % Gaussian levels
    contour(ax, X1, X2, rho, levels, 'LineColor','w', 'LineWidth',1.2);

    % Positivity boundary: e^{x1}+e^{x2} = K  ->  x2 = log(K - e^{x1}) for x1 < log K
    x1b  = linspace(min(x1span), min(log(K)-1e-8, max(x1span)), 800);
    mask = exp(x1b) < K;
    x1b  = x1b(mask);
    x2b  = log(K - exp(x1b));
    in   = (x2b >= min(x2span)) & (x2b <= max(x2span));
    plot(ax, x1b(in), x2b(in), 'k--', 'LineWidth', 1.4);             % dashed boundary

    axis(ax, 'equal'); axis(ax, 'tight'); grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'x_1'); ylabel(ax, 'x_2');
    title(ax, sprintf('f(x) contours  +  N(0, I_2),  K = %g', K));
    colorbar(ax);
    hold(ax,'off');
end

set(fig, 'Renderer','opengl');   % fast interaction/rendering

% Optional: export tightly-cropped PDF (raster inside the PDF for speed)
% exportgraphics(fig, 'base_problem_contours.pdf', 'ContentType','image', 'Resolution', 200);
