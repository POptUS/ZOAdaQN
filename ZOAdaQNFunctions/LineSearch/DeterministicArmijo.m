function [Val, g_new, alpha, NLineEvals] = DeterministicArmijo(Obj, inputvals, f, d, grad, w, c, rho, alpha, tolr)
    % This function perfroms armijo linesearch algorithm: Algorithm 3.1 given in
    % Numerical optimization book by Stephen Wright and Jorge Nocedal
    %
    % The inputs are
    % Obj: This is declared as a function which gives the functional value and
    %      the gradient (estimated via ) value.
    % Obj.func: Gives function and gradient values of the whole function
    % inputvals: This is the input that need to be given into the Obj for
    %            getting the functional values.
    % f: function value at w
    % d: search direction
    % grad: gradient estimate at the iterate w
    % w: current iterate
    % c: constant in Armijo line search
    % rho: contraction factor
    % alpha: initial trial alpha
    % tolr: lower bound on alpha
    %
    % The outputs are:
    % Val: function value at w + alpha*d
    % g_new: gradient estimate at w + alpha*d
    % alpha: stepsize value that satisifies armijo condition
    % NLineEvals: Total number of line search evaluations done

    % Default Values
    if nargin < 9
        alpha = 1;
    end
    if nargin < 10
        tolr = 10^-8;
    end

    % Initialization
    NLineEvals = 0;
    const = grad' * d * c;
    [Val, g_new] = Obj.func(w + alpha * d, inputvals);
    NLineEvals = NLineEvals + 1;

    % Checking for the condition
    while Val > f + alpha * const
        alpha = alpha * rho;
        [Val, g_new] = Obj.func(w + alpha * d, inputvals);
        NLineEvals = NLineEvals + 1;
        if alpha <= tolr %|| NLineEvals>=10
            alpha = tolr;
            break
        end
    end
end
