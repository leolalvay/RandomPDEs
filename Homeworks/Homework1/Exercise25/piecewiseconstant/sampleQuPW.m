%This function generates samples for Q(u_h) for the piecewise rf

function Q = sampleQuPW(M, sigma, N, l)
% M is the number of samples
% sigma is the sigma parameter in a piecewise constant
% Yn is the uniform vector
%l is the refinement mesh level

Q = zeros(M,1);
for i=1:M
    Yn = unifrnd(-sqrt(3), sqrt(3), [N 1]);
    Q(i) = fesolverpiecewise(l, sigma,Yn);
end

end



