function Q = sampleQuPWgivenY(sigma,ell, Y)
% Y is M-by-N (each row is a realization of [Y1,...,YN])
% Returns Q(m) = Q(u_h) in the level ell for each row of Y

    M = size(Y,1);
    Q = zeros(M,1);
    for m = 1:M
        Yn = Y(m,:).';                 
        Q(m) = fesolverpiecewise(ell, sigma, Yn);
    end
end
