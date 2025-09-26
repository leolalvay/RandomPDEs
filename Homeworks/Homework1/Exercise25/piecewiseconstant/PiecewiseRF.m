% Model 1: Piecewise contanst coefficient
% PiecewiseRF produces a sample of size length(x) of the 
% random field a at the points in vector x
% N is the number of terms in the sum that defines a
function output = PiecewiseRF(x,sigma, Y_n)

% Parameters
%sigma = 0.1; % between 0 and 1/sqrt(3) = 0.5774
%Y = unifrnd(-sqrt(3), sqrt(3), [N 1]);
P = length(x);
N = length(Y_n);
% Using the piecewise function to get the output
output = zeros(P,1);
for k = 1:P
    output(k) = rfa(x(k),N,Y_n,sigma);
end

% Piecewise function
    function value = rfa(x,N,Y_n,sigma)
xhat = (0:1/N:1)';
sum = 0;
for i = 1:N
    sum = sum + Y_n(i)*(x > xhat(i) && x < xhat(i+1));
end
value  = 1 + sigma*sum;
end
end
%Since for the solution of the FE method we will not evaluate a in the 
%coarse nodes, we impose stric inequalities: x > xhat(i) && x < xhat(i+1)