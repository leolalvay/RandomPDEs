clear; clc;
% this code produces the 3 histograms with M samples for different values of  N
fig = figure('Units','inches','Position',[1 1 9.5 3.8],'Color','w');

M = 1e3;
N = [10 20 40];
% % Parameters for Piecewise model
ell = 2;
sigma = 0.15;
i = 0;
for l = 1:numel(N)
        i = i+1;
        subplot(1,3,i)
        myplotter(M, ell, N(i), sigma)
end


sgtitle(sprintf('Histogram of samples $Q(u_h)$ with $\\sigma=%.2f$, $M=%.0f$', ...
                sigma, M), ...
        'Interpreter','latex','FontWeight','bold','FontSize',16);
exportgraphics(fig, 'histograms.pdf', 'ContentType','vector', 'BackgroundColor','white');

function myplotter(M, ell, N, sigma)
%h = 1/(N*2^ell);
% Monte Carlo estimate of Q with Piecewise model
Qu_hsamples = sampleQuPW(M,sigma,N,ell);
%Qu_hsamples = zeros(M,1);
%for i =1:M
 %   Y_n = unifrnd(-sqrt(3), sqrt(3), [N 1]);
  %  Qu_hsamples(i) = fesolverpiecewise(ell, sigma, Y_n);
%end
% Calculate the mean of the data
EQu_h = mean(Qu_hsamples);
VQu_h = var(Qu_hsamples);
histogram(Qu_hsamples, 'Normalization', 'probability','BinMethod','scott', 'EdgeColor', 'k', 'FaceColor', [0.28 0.38 0.54],'FaceAlpha',0.90);
%xline(EQu_h, 'b', 'LineWidth', 2, 'FontSize', 16);
xlabel(sprintf(' $N=%d$', N), 'Interpreter','latex','FontSize',16);
ylabel('Relative frequency', 'Interpreter', 'latex', 'FontSize', 16);
title([' $Var(Q(u_h))=$',num2str(VQu_h,3)], 'Interpreter', 'latex', 'FontSize', 16);
%legend('Data', 'Mean');
grid on;
%text(EQu_h, 0, sprintf('%.2f', EQu_h), 'FontSize', 12, 'Color', 'b');
end