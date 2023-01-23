function [f, g] = MAD_Expec(w, inputvals)
    % This function computes the objective function and gradient values of the
    % function f(x) = E[\|Ax - b - \zeta\|_1] (zeta is random variable)
    % inputs: w - input parameters
    %        inputvals
    %                A: Data matrix
    %                b: Data vector
    %                Zeta: 'Uni'
    %                U: 'upper bound'
    %                L: 'lower bound'

    A = inputvals.A;
    b = inputvals.b;
    Zeta = inputvals.Zeta;
    U = inputvals.U;
    L = inputvals.L;

    d = length(w);
    m = length(b);

    v = A * w - b;
    switch Zeta
        case 'Uni'
            fi = zeros(m, 1);
            for i = 1:m
                if v(i) > U
                    fi(i) = 1 * v(i);
                elseif v(i) < L
                    fi(i) = -1 * v(i);
                else
                    fi(i) = (v(i)^2 + 1) / 2;
                end
            end
            f = sum(fi);
            if nargout > 1
                gi = zeros(d, m);
                for i = 1:m
                    if v(i) > U
                        gi(:, i) = A(i, :);
                    elseif v(i) < L
                        gi(:, i) = -1 * A(i, :);
                    else
                        gi(:, i) = v(i) * A(i, :);
                    end
                end
                g = sum(gi, 2);
                % DFO gradient
                %             func=@(x)MAD_Expec(x, inputvals);
                %             DFOoptions.Method=inputvals.DFOMethod;
                %             h=Cal_DFOinterval(inputvals,DFOoptions);
                %             DFOoptions.interval=h;
                %             g=DFOGrad(w, func, DFOoptions);
            end
    end
