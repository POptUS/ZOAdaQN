function [f, g] = MAD_BatchDFO_PFLoss(w, inputvals, batch)
    % This function computes the objective function and gradient values of the
    % function f(x) = E[\|Ax - b - \zeta\|_1] (zeta is random variable)
    % inputs: w - input parameters
    %        inputvals
    %                A: Data matrix
    %                b: Data vector
    %                Zeta: 'Uni'
    %                U: 'upper bound'
    %                L: 'lower bound'
    % batch: batch
    A = inputvals.A;
    b = inputvals.b;
    Zeta = inputvals.Zeta;
    U = inputvals.U;
    L = inputvals.L;

    d = length(w);
    m = length(b);
    S = length(batch);

    v = A * w - b;
    randmat = inputvals.randmat;
    randmat = randmat(batch, :);
    switch Zeta
        case 'Uni'
            fi = zeros(m, 1);
            g = zeros(d, 1);
            for i = 1:S
                % rand('seed',batch(i));
                % r=L + rand(m,1)*(U-L);
                r = L + randmat(i, :)' * (U - L);
                fi = fi + abs(v - r);
                if nargout > 1
                    for j = 1:m
                        if v(j) > U
                            g = g + A(j, :)';
                        elseif v(j) < L
                            g = g - A(j, :)';
                        else
                            g = g + (v(j)) * A(j, :)';
                        end
                    end
                end
            end
            f = sum(abs(fi)) / S;
            g = g / S;
            if nargout > 1
                % g=zeros(d,1);
                func = @(x)MAD_BatchDFO(x, inputvals, batch);
                DFOoptions.Method = inputvals.DFOMethod;
                h = Cal_DFOinterval(inputvals, DFOoptions);
                DFOoptions.interval = h;
                if strcmp(DFOoptions.Method, 'GS') == 1 || strcmp(DFOoptions.Method, 'SS') == 1
                    DFOoptions.dirs = inputvals.DFODirs;
                end
                g = DFOGrad(w, func, DFOoptions);
            end
    end
