function [g, fvals] = DFOGrad(x, func, options)
    % This script estimates the gradient of any function using only function
    % values. The different approaches are:
    % 1. FD (forward)
    % 2. CFD (Central)
    % 3. Gaussian Smoothing
    % 4. Sphere Smoothing
    % 5. Linear Interpolation methods
    % The inputs are
    %   x: variable
    % func: function object whose input is x and output is function values
    % inputvals: input parameters for the function object
    % options:
    %        method: specify 'FD' or 'GS' or 'CFD' or 'LID'
    %                FD: Forward finite differecnes
    %                CFD: Central finite differences
    %                GS: Gaussian smoothing without scaling the directions
    %                SS: sphere smoothing
    %        interval: interval length
    %        directions: # number of directions
    % outputs:
    %   g: gradient estimator
    %   funvals: function values evaluated.

    d = length(x);
    h = options.interval; % *max(1, norm(x));
    switch options.Method
        case 'FD'
            fvals = zeros(d + 1, 1);
            g = zeros(d, 1);
            fvals(1) = func(x);
            for i = 1:d
                e = zeros(d, 1);
                e(i) = 1;
                fvals(i + 1) = func(x + h * e);
                g(i) = (fvals(i + 1) - fvals(1)) / h;
            end
        case 'CFD'
            fvals = zeros(2 * d + 1, 1);
            g = zeros(d, 1);
            fvals(1) = func(x);
            for i = 1:d
                e = zeros(d, 1);
                e(i) = 1;
                fvals(2 * i) = func(x - h * e);
                fvals(2 * i + 1) = func(x + h * e);
                g(i) = (fvals(2 * i + 1) - fvals(2 * i)) / (2 * h);
            end
        case 'GS'
            % working on the code
            dr = options.dirs;
            fvals = zeros(dr, 1);
            fvals(1) = func(x);
            g = zeros(d, 1);
            for i = 1:dr
                v = randn(d, 1);
                v = v / norm(v);
                fvals(i + 1) = func(x + h * v);
                g = g + (fvals(i + 1) - fvals(1)) * v / h;
            end
            g = g / dr;

        case 'SS'
            dr = options.dirs;
            fvals = zeros(dr, 1);
            fvals(1) = func(x);
            g = zeros(d, 1);
            for i = 1:dr
                v = randn(d, 1);
                v = v / norm(v);
                fvals(i + 1) = func(x + h * v);
                g = g + (fvals(i + 1) - fvals(1)) * v / h;
            end
            g = g * d / dr;

        case 'GS-C'
            % working on the code
            fvals(1) = func(x);
            % g=zeros(d,1);
            v = randn(d, 1);
            % v=v/norm(v);
            fvals(2) = func(x + h * v);
            g = (fvals(2) - fvals(1)) * v / h;

    end
