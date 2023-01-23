function [Val, g_new, alpha, NLineEvals] = Armijo(Obj, inputvals, batch, f, d, grad, w, c, rho, alpha, tolr, c2)
    % This function is to perfrom armijo linesearch for the Adaptive sampling
    % quasi-Newton algorithm.
    %
    % The inputs are
    % Obj: This is declared as a function which gives the functional value and
    %      the gradient (estimated via ) value for full and batch cases.
    % Obj.func: Gives function and gradient values of the whole function
    % Obj.funcBatch: Gives the sample average objective and sample average
    %                gradient.
    % inputvals: This is the input that need to be given into the Obj for
    %            getting the functional values.
    % batch: Consists of the indices in the set $S_k$
    % f: function value at w
    % d: search direction
    % grad: gradient estimate at the iterate w
    % w: current iterate
    % c: constant in Armijo line search
    % c2:relaxed armijo
    % rho: contraction factor
    % alpha: initial trial alpha
    % tolr: lower bound on alpha
    %
    % The outputs are:
    % Val: function value at w + alpha*d
    % g_new: gradient estimate at w + alpha*d
    % alpha: stepsize value that satisifies armijo condition
    % NLineEvals: Total number of line search evaluations done

    % default settings for alpha, tolr, c2
    if nargin <= 10
        alpha = 1;
    end
    if nargin <= 11
        tolr = 10^-6;
    end

    if nargin < 12
        c2 = 10^-16;
    end
    % Initial estimates
    NLineEvals = 0;
    const = grad' * d * c;
    [Val, g_new] = Obj.funcBatch(w + alpha * d, inputvals, batch);
    NLineEvals = NLineEvals + 1;

    % Checking the loop
    while Val > f + alpha * const + c2
        alpha = alpha * rho;
        [Val, g_new] = Obj.funcBatch(w + alpha * d, inputvals, batch);
        NLineEvals = NLineEvals + 1;
        if alpha <= tolr %|| NLineEvals>=10
            % alpha=10^-6;
            break
        end
    end
end
