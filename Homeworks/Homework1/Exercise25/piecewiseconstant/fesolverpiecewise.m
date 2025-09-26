function Q = fesolverpiecewise(ell, sigma, Y_n)
% Parameters
N = length(Y_n);
I = N*2^ell;
h = 1/I;
% Uniform grid
x=(0:h:1)';
% Midpoints 
midpoints = (x(1:end-1) + x(2:end)) / 2;
% Sampling the Piecewise random field at the midpoints
realization = PiecewiseRF(midpoints(1:end), sigma, Y_n);
% Building up the matrix
udiag = diag(-realization(2:end-1)/h/h,1); % This is A_{i,i+1}
mdiag = diag((realization(1:end-1) + realization(2:end))/h/h); % This is A_{i,i}
ldiag = diag(-realization(2:end-1)/h/h,-1); % This is A_{i,i-1}

A = udiag + mdiag + ldiag; % Matrix A
F = f(x(2:end-1)); % Vector F
u = A\F; % FEM solution at the internal nodes

Q = h*sum(u); %Estimate of Q(u_h)

function output = f(x)
output = 4*pi^2*cos(2*pi*x);
end
end