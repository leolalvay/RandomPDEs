M = 1e6; % Max. number of realization
N = 20; % Dimension of the problem
u = rand(M,N); f = exp(sum(u')');
run_aver = cumsum(f)./(((1:M)')*(exp(1)-1)^N);
figure, plot(1:M, run_aver), xlabel 'M'
figure,plot(1:M,(run_aver-1)), xlabel 'M'
figure,plot(1:M,abs(run_aver-1)), xlabel 'M'
figure,semilogy(1:M,abs(run_aver-1)), xlabel 'M',
