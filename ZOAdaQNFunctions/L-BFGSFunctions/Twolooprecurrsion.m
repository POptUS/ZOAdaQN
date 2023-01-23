function [r] = Twolooprecurrsion(H0, history, grad, k)

    % This function implements the L-BFGS two loop recusion. Tt
    % implements the Algorithm 7.4 given in the book
    % Numerical Optimization by Stephen Wright and Jorge Nocedal
    if k ~= 0
        q = grad;
        y = history.y;
        s = history.s;
        rho = history.rho;
        m = size(y, 1);
        alphai = zeros(m, 1);
        for i = m:-1:1
            alphai(i) = rho(i) * s(i, :) * q;
            q = q - alphai(i) * y(i, :)';
        end
        r = H0 * q;
        for i = 1:m
            beta = rho(i) * y(i, :) * r;
            r = r + s(i, :)' * (alphai(i) - beta);
        end
    else
        r = H0 * grad;
    end

end
