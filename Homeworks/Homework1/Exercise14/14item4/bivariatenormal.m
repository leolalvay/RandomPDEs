%% Bivariate Normal Playground (contours + surface)
% Edit the parameters in the block below and re-run.
% Stable, toolbox-free, vectorized. Rotate the 3-D surface with the mouse.

clear; clc;

%% -------------------- Parameters --------------------
% Option A: specify by variances + correlation rho
mu  = [-4.2; 0.0];           % [mu_x1; mu_x2]
sig1 = 1.0;             % std dev of X1  (sigma_x1)
sig2 = 2.0;             % std dev of X2  (sigma_x2)
rho  = -0.7;            % correlation in (-1, 1)
USE_SIGMA_DIRECT = false;

% Option B: (alternative) give Sigma directly (positive-definite 2x2)
% If you set USE_SIGMA_DIRECT=true, this Sigma overrides sig1/sig2/rho.
Sigma_direct = [1.2, -0.5; -0.5, 1.2];

%% -------------------- Build Sigma --------------------
if ~USE_SIGMA_DIRECT
    % Sigma from sigmas and rho
    Sigma = [sig1^2,      rho*sig1*sig2; ...
             rho*sig1*sig2, sig2^2     ];
else
    Sigma = Sigma_direct;
end

% Basic sanity checks
if ~isequal(size(mu),[2,1])
    error('mu must be a 2x1 column vector.');
end
if any(diag(Sigma) <= 0)
    error('Sigma must have positive diagonal entries.');
end
if ~all(eig(Sigma) > 0)
    error('Sigma must be positive-definite.');
end

% Report correlation implied by Sigma (helpful when Sigma_direct is used)
rho_implied = Sigma(1,2) / sqrt(Sigma(1,1)*Sigma(2,2));
fprintf('Using Sigma = \n'); disp(Sigma);
fprintf('Implied correlation rho = %+0.4f\n', rho_implied);

%% -------------------- Grid (auto-sized around mu) --------------------
% Make the window adapt to the scale/orientation (≈ ±4 std along principal axes)
[V, D] = eig(Sigma);                           % principal axes
rad_max = 4 * sqrt(max(diag(D)));             % radius in the largest direction
pad = 0.5;                                    % small padding
x = linspace(mu(1)-rad_max-pad, mu(1)+rad_max+pad, 201);
y = linspace(mu(2)-rad_max-pad, mu(2)+rad_max+pad, 201);
[X1, X2] = meshgrid(x, y);

%% -------------------- Evaluate pdf (stable) --------------------
% pdf(x) = (2*pi)^(-1) * det(Sigma)^(-1/2) * exp(-0.5*(x-mu)'*inv(Sigma)*(x-mu))
% We compute it stably via Cholesky: Sigma = R'*R
Z = bvn_pdf_chol(X1, X2, mu, Sigma);   % defined below

%% -------------------- Plots --------------------
% Contours
fC = figure('Color','w');
contourf(X1, X2, Z, 15, 'LineColor', 'none');  % 15 filled levels
axis equal tight;
colorbar;
xlabel('x_1'); ylabel('x_2');
title('Bivariate Normal Contours');
set(fC,'Renderer','opengl');

% Surface (rotate with mouse)
fS = figure('Color','w');
surf(X1, X2, Z, 'EdgeColor','none');
xlabel('x_1'); ylabel('x_2'); zlabel('pdf');
title('Bivariate Normal Surface');
view(35, 30);
camlight headlight; lighting gouraud;
axis tight;
rotate3d on;
set(fS,'Renderer','opengl');

% ------- Optional: quick, tightly-cropped export (raster inside PDF) ------
% exportgraphics(fC, 'bvn_contours.pdf', 'ContentType','image','Resolution',200);
% exportgraphics(fS, 'bvn_surface.pdf',  'ContentType','image','Resolution',200);

%% -------------------- Local function --------------------
function Z = bvn_pdf_chol(X1, X2, mu, Sigma)
% Bivariate normal pdf using Cholesky (stable; no toolbox).
% X1, X2: meshgrid arrays (same size). mu: 2x1. Sigma: 2x2 PD.

    % Cholesky: Sigma = R'*R with R upper-triangular
    [R, p] = chol(Sigma);
    if p ~= 0
        error('Sigma is not positive-definite (chol failed).');
    end

    % Form centered data, stack into N-by-2
    d1 = X1(:) - mu(1);
    d2 = X2(:) - mu(2);
    D  = [d1, d2];

    % Solve y = D / R, then quadratic form is sum(y.^2,2)
    Y = D / R;
    q = sum(Y.^2, 2);                    % Mahalanobis^2

    % Normalization: (2*pi)^(-1) * 1/prod(diag(R))
    c = 1 / (2*pi) / prod(diag(R));

    Z = c * exp(-0.5 * q);
    Z = reshape(Z, size(X1));
end
