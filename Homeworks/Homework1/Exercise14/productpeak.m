clc; clear;

% ------------ Parameters ------------
M  = 1e6;
w  = 0.5;
Ns = [2 20];      % N=2 and N=20

% Final size of the pdf figure
Wtot = 9.5;   % total (e.g. for \textwidth)
Htot = 4.5;   % total height

% Prevent opening windows (optional)
set(groot,'defaultFigureVisible','off');

% ------------ Combined figure ------------
fig = figure('Units','inches','Position',[1 1 Wtot Htot],'Color','w');
tlo = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

allY = [];   % to unify Y-lims across both subplots

for k = 1:numel(Ns)
    N = Ns(k);
    S = compute_curves(M,N,w);

    ax = nexttile;
    loglog(S.m(10:end), max(S.err(10:end), eps), 'LineWidth', 1.1, 'Displayname','RE'); hold on;
    loglog(S.m(10:end), max(S.clt(10:end), eps), 'LineWidth', 1.1, 'DisplayName','CLT');
    loglog(S.m(10:end), max(S.be(10:end),  eps), 'LineWidth', 1.1, 'DisplayName','BE');

    grid on; xlim([1e1, M]); xlabel('M');
    if k==1, ylabel('Relative error'); end
    title(sprintf('Produc Peak, N = %d', N));  % same naming you use

    legend(ax,'Interpreter','latex','Location','best');

    % Store Y data for unified limits
    allY = [allY;
            max(S.err(10:end),eps);
            max(S.clt(10:end),eps);
            max(S.be(10:end), eps)];
end

% Unify Y limits on both axes (recommended)
ylims = [min(allY), max(allY)];
axArr = findall(fig,'Type','axes');
set(axArr,'YLim', ylims);

% ------------ Export vector PDF with exact figure size ------------
outname = 'product_peak_N2_N20.pdf';
exportgraphics(fig, outname, 'ContentType','vector', 'BackgroundColor','none');

close(fig);  % optional
set(groot,'defaultFigureVisible','on');

fprintf('Saved %s\n', outname);

% ========= Helper =========
function S = compute_curves(M,N,w)
    % Your integrand and statistics
    u      = rand(M,N);
    cN     = 7.25/N;
    Iexact = ( 2*cN * atan(cN/2) )^N;
    f      = prod( (cN^2) ./ ( 1 + (cN*(u - w)).^2 ), 2 );

    m = (1:M).';
    run_aver = cumsum(f) ./ m;

    % global sigma (MLE, divides by M)
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

    den  = max(abs(Iexact), realmin('double'));  % avoid division by 0
    S.m   = m;
    S.err = abs(run_aver - Iexact) / den;   % <-- relative error
    S.clt = clt_m / den;                    % <-- relative CLT band
    S.be  = be_m  / den;                    % <-- relative BE band

end
