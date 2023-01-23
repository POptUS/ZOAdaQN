function [g, fvals] = DFOGradBatch(x, funcbatch, options)
    % This script estimates the gradient of any function using only function
    % values. The different approaches are:
    % 1. FD
    % 2. CFD
    % 3. Gaussian Smoothing
    % 4. Sphere Smoothing
    % 5. Linear Interpolation methods
    % The inputs are
    % func: Function
    % inputvals: input parameters for the function object
    % options:
    %        method: specify 'FD' or 'GS' or 'CFD' or 'LID'
    %        interval: interval length
    %        directions: # number of directions
    % outputs:
    %   g: gradient estimator
    %   funvals: function values evaluated.
    d = length(x);
    h = options.interval; % *max(1, norm(x));
    batch = options.funcbatch;
    S = length(batch);

    function val = func(x)
        val = 0;
        for ij = 1:S
            val = val + funcbatch(x, batch(ij));
        end
        val = val / S;
    end

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
            % working on the code
            %         fvals(1)=func(x);
            %         %g=zeros(d,1);
            %         v=randn(d,1);
            %         v=v/norm(v);
            %         fvals(2)=func(x + h*v);
            %         g= (fvals(2) - fvals(1))*v/h;
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

    end
end
