function [gamma] = GetH0(history, k)

    % This function is to get the initial H0 that is to be used in the L-BFGS
    % two loop recursion. This function implements equation 7.20 from the book
    % Numerical Optimization by Stephen Wright and Jorge Nocedal
    if k == 0
        gamma = 1;
    else
        y = history.y;
        s = history.s;
        gamma = s(k, :) * y(k, :)' / (y(k, :) * y(k, :)');
    end
    % H0=gamma*eye(history.n);
end
