function h = plot_f_and_gaussian(K, mu, Sigma, gridN)
%PLOT_F_AND_GAUSSIAN  Visualize f(x) and overlay a N(mu,Sigma) on top of f.
%   h = plot_f_and_gaussian(K, mu, Sigma, gridN)
%
% INPUTS
%   K      : scalar > 0 (parameter in f)
%   mu     : 2x1 vector [mu1; mu2]           (your optimized mean)
%   Sigma  : 2x2 positive-definite covariance (you choose it visually)
%   gridN  : (optional) grid resolution (default 301). Larger = smoother/slower.
%
% OUTPUT
%   h : struct of handles to the figure and axes (useful if you want to export).
%
% WHAT IT DRAWS
%   Left  tile: surface of f(x) = max(exp(x1)+exp(x2)-K,0) * exp(-(x1^2+x2^2)/2)
%   Right tile: filled contour of f, with:
%               - dashed curve  e^{x1}+e^{x2}=K  (positivity boundary),
%               - white contour lines of N(mu,Sigma) overlaid.
%
% NOTES
% * The constant 1/(2*pi) is omitted in f (does not change the shape or argmax).
% * The Gaussian pdf IS normalized (1/(2*pi*sqrt(det(Sigma)))), but only its
%   contour shape matters; lines are overlaid to compare orientation/width.
% * No toolboxes are required: Gaussian pdf is evaluated via Cholesky.
%
% EXAMPLE
%   mu = [log(3/2); log(3/2)];  Sigma = [1.2 -0.5; -0.5 1.2];
%   plot_f_and_gaussian(3, mu, Sigma, 301);

    if nargin < 4, gridN = 301; end
    validate_inputs(K, mu, Sigma);

    % ---- Auto-size a nice window from Sigma and mu (±4 std along major axis) ----
    [V,D] = eig(Sigma);                           % principal axes
    Rmax   = 4*sqrt(max(diag(D)));                % radius along the widest axis
    pad    = 0.5;
    x1 = linspace(mu(1)-Rmax-pad, mu(1)+Rmax+pad, gridN);
    x2 = linspace(mu(2)-Rmax-pad, mu(2)+Rmax+pad, gridN);
    [X1, X2] = meshgrid(x1, x2);

    % ---- f(x) on the grid (no 1/(2*pi) factor on purpose) ----
    F = max(exp(X1) + exp(X2) - K, 0) .* exp(-0.5*(X1.^2 + X2.^2));

    % ---- Gaussian pdf via Cholesky (stable; no toolbox) ----
    G = bvn_pdf_chol(X1, X2, mu, Sigma);

    % ---- Figure with two tiles ----
    h.fig = figure('Color','w');
    tl = tiledlayout(h.fig, 1, 1, 'TileSpacing','compact','Padding','compact');

    % ---- Left: surface of f ----
    %h.ax1 = nexttile(tl);
    %surf(h.ax1, X1, X2, F, 'EdgeColor','none');
    %xlabel(h.ax1, 'x_1'); ylabel(h.ax1, 'x_2'); zlabel(h.ax1, 'f(x)');
    %title(h.ax1, sprintf('Surface of f(x),  K = %g', K));
    %view(h.ax1, 35, 28); camlight(h.ax1, 'headlight'); lighting(h.ax1, 'gouraud');
    %axis(h.ax1, 'tight'); grid(h.ax1, 'on'); box(h.ax1, 'on');

    % ---- Right: f contour + Gaussian contours + positivity boundary ----
    h.ax = nexttile(tl);
    contourf(h.ax, X1, X2, F, 20, 'LineColor','none');   % 20 filled levels
    hold(h.ax, 'on');

    % White contour lines of the Gaussian (levels as fractions of its max)
    levels = max(G(:)) .* linspace(0.1, 0.9, 6);         % 6 light lines
    contour(h.ax, X1, X2, G, levels, 'LineColor','w', 'LineWidth', 1.2);

    % Positivity boundary e^{x1}+e^{x2}=K  ->  x2 = log(K - e^{x1}), for e^{x1}<K
    x1b = linspace(min(x1), min(log(K)-1e-8, max(x1)), 800);
    mask = exp(x1b) < K;   x1b = x1b(mask);
    x2b = log(K - exp(x1b));
    % keep boundary inside window
    in  = x2b >= min(x2) & x2b <= max(x2);
    plot(h.ax, x1b(in), x2b(in), 'k--', 'LineWidth', 1.5);

    hold(h.ax, 'off');
    axis(h.ax, 'equal'); axis(h.ax, 'tight');
    xlabel(h.ax, 'x_1'); ylabel(h.ax, 'x_2');
    title(h.ax, sprintf('f(x) contours  +  N(\\mu,\\Sigma),  K = %g', K));
    colorbar(h.ax); grid(h.ax, 'on'); box(h.ax, 'on');

    % Fast renderer (useful for smooth interaction/exports)
    set(h.fig, 'Renderer', 'opengl');

    % -------- Optional: quick, tightly-cropped export (raster inside PDF) ------
    % exportgraphics(h.fig, 'f_and_gaussian.pdf', 'ContentType','image','Resolution',200);
end

% ================== helpers ==================
function validate_inputs(K, mu, Sigma)
    if ~isscalar(K) || ~(K > 0), error('K must be a positive scalar.'); end
    if ~isequal(size(mu), [2,1]), error('mu must be a 2x1 column vector.'); end
    if ~isequal(size(Sigma), [2,2]), error('Sigma must be 2x2.'); end
    if any(diag(Sigma) <= 0), error('Sigma diagonal must be positive.'); end
    e = eig((Sigma+Sigma.')/2);                 % symmetrize then check PD
    if any(e <= 0), error('Sigma must be positive-definite.'); end
end

function Z = bvn_pdf_chol(X1, X2, mu, Sigma)
%BvN_PDF_CHOL  N(mu,Sigma) pdf on a grid using Cholesky (stable, no toolbox).
% pdf(x) = (2*pi)^(-1) det(Sigma)^(-1/2) * exp(-0.5*(x-mu)'*inv(Sigma)*(x-mu))
% We avoid inv(Sigma): with Sigma = R'*R, solve y = (x-mu)/R, then use ||y||^2.

    [R, p] = chol((Sigma+Sigma.')/2);  %#ok<NASGU>  % ensure symmetry; p>0 if not PD
    if p ~= 0, error('Sigma is not positive-definite (chol failed).'); end

    % center and stack grid points into N-by-2
    d1 = X1(:) - mu(1);
    d2 = X2(:) - mu(2);
    D  = [d1, d2];

    Y = D / R;                % solves y = D * inv(R) (faster and stable)
    q = sum(Y.^2, 2);         % squared Mahalanobis distance

    c = 1/(2*pi)/prod(diag(R));
    Z = c * exp(-0.5*q);
    Z = reshape(Z, size(X1));
end

