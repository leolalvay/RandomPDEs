clc; clear;

% ------------ Parámetros ------------
M  = 1e6;
w  = 0.5;
Ns = [2 20];      % N=2 y N=20
rng(42);          % reproducible (opcional)

% Tamaño final del PDF (figura completa) en pulgadas
Wtot = 9.5;   % ancho total (p. ej. pensado para \textwidth)
Htot = 4.5;   % alto total

% Evitar que abra ventanas (opcional)
set(groot,'defaultFigureVisible','off');

% ------------ Figura combinada ------------
fig = figure('Units','inches','Position',[1 1 Wtot Htot],'Color','w');
tlo = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

allY = [];   % para unificar Y-lims en ambos subplots

for k = 1:numel(Ns)
    N = Ns(k);
    S = compute_curves(M,N,w);

    ax = nexttile;
    loglog(S.m(10:end), max(S.err(10:end), eps), 'LineWidth', 1.1, 'DisplayName','$|\epsilon_M|$'); hold on;
    loglog(S.m(10:end), max(S.clt(10:end), eps), 'LineWidth', 1.1, 'DisplayName','CLT');
    loglog(S.m(10:end), max(S.be(10:end),  eps), 'LineWidth', 1.1, 'DisplayName','BE');

    grid on; xlim([1e1, M]); xlabel('M');
    if k==1, ylabel('Statistical error'); end
    title(sprintf('Oscillatory, N = %d', N));  % mismo naming que usas

    legend(ax,'Interpreter','latex','Location','best');

    % Guardar datos de Y para límites unificados
    allY = [allY;
            max(S.err(10:end),eps);
            max(S.clt(10:end),eps);
            max(S.be(10:end), eps)];
end

% Unificar límites de Y en ambos ejes (recomendado)
ylims = [min(allY), max(allY)];
axArr = findall(fig,'Type','axes');
set(axArr,'YLim', ylims);

% ------------ Exportar PDF vectorial del tamaño exacto de la figura ------------
outname = 'Oscillatory_N2_N20.pdf';
exportgraphics(fig, outname, 'ContentType','vector', 'BackgroundColor','none');

close(fig);  % opcional
set(groot,'defaultFigureVisible','on');

fprintf('Saved %s\n', outname);

% ========= Helper =========
function S = compute_curves(M,N,w)
    % Tu integrando y estadísticos
    u = rand(M,N);
    cN=9/N;
    Iexact = osc_exact_equal(N,w);
    s = cN * sum(u, 2);                 % Mx1: (9/N)*row-sum
    f = cos(2*pi*w + s);               % Mx1: integrand at each sample

    m = (1:M).';
    run_aver = cumsum(f) ./ m;

    % sigma global (MLE, divide por M)
    sigma = std(f, 1);

    % CLT 95%
    alpha  = 0.05;
    zalpha = -sqrt(2)*erfcinv( 2*(1 - alpha/2) );  % ~1.96
    clt_m  = zalpha * sigma ./ sqrt(m);

    % Berry–Esseen 95%
    mu_all    = mean(f);
    m3abs_all = mean(abs(f - mu_all).^3);
    lambda3   = m3abs_all / max(sigma^3, realmin);

    CBE   = 30.51175;
    phi_z = exp(-0.5*zalpha^2)/sqrt(2*pi);

    Km    = (CBE * lambda3) ./ sqrt(m);
    C0_BE = zalpha + Km ./ (2*phi_z*(1+zalpha)^3);
    be_m  = C0_BE .* sigma ./ sqrt(m);

    % Salida
    S.m   = m;
    S.err = abs(run_aver - Iexact);
    S.clt = clt_m;
    S.be  = be_m;
end

%Exact Solution
function I = osc_exact_equal(N, w1)
cN   = 9/N;
theta = 2*pi*w1 + 9/2;                        % since sum c_n = 9
I     = cos(theta) * (sin(0.5*cN)/(0.5*cN))^N; % product collapses to a power
%I = cos(theta) * (sinc(cN/(2*pi)))^N;
end