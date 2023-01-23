function [h] = Cal_DFOinterval(inputvals, Options)
    % This function computes the adaptive DFO interval for problems in cuter
    % dataset and logistic regression

    h = inputvals.DFOMachinePrecision;
    loss = inputvals.loss;
    method = inputvals.DFOMethod;
    factor = inputvals.DFOIntervalFactor;
    interval = inputvals.DFOInterval;
    if strcmp(interval, 'Adaptive') == 1
        switch loss
            case {'CuterDFO', 'Cuter'}
                problem = Options.DFOproblem;
                sigma = Options.DFOsigma;
                adaptive = inputvals.DFOIntervalAdaptive;
                switch adaptive
                    case 'Exact'
                        F = Options.DFOfsmooth;
                    otherwise
                        disp('Error: No function value is given for adaptively choosing finite difference interval \n');
                        disp('And so, h is chosen based on machine precision');
                        return
                end
                switch problem.probtype
                    case 'relnormal'
                        h = h + 2 * sigma * F;
                    case 'absnormal'
                        m = problem.m;
                        h = h + 2 * sigma * sqrt(F) + sqrt(2) * m * sigma^2;
                end
        end
    end

    switch method
        case {'FD', 'GS', 'SS', 'GS-C'} % Forward finite differences, gaussian smoothing, sphere smoothing
            h = sqrt(h) / sqrt(factor);
        case 'CFD' % Central finite-differences
            h = h^(1 / 3) / factor;
    end

end
