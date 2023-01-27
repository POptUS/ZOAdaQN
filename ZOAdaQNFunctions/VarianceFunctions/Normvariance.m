function [var, skew] = Normvariance(batch, mean, Obj, inputvals, w)
    % This function calculates the sample variance and sample skew estimates
    % of the norm condition.
    %
    % The inputs are:
    % batch: indices in S_k^v used for sample estimates
    % mean: average of the sample estimates
    % Obj: This is declared as a function which gives the functional value and
    %      the gradient (estimated via finite-differences) value batch cases.
    % Obj.funcBatch: Gives the sample average objective and sample average
    %                gradient.
    % inputvals: This is the input that need to be given into the Obj for
    %            getting the functional values.
    % w: the current iterate
    %
    %
    % The outputs are:
    %  var: sample variance (2nd moment)
    % skew: sample skew estimate (3rd moment)

    % Problem parameters
    n = length(batch); % number of data points;
    p = length(mean); % number of variables (features);

    % initialization
    var = zeros(p, 1);
    skew = zeros(p, 1);

    for i = 1:n
        [~, g] = Obj.funcBatch(w, inputvals, batch(i));
        var = var + (g - mean).^2;
        if nargout > 1
            % computing the skewness
            skew = skew + (g - mean).^3;
        end
    end

    % New way of implementing it
    % [g]=Obj.funcBatchVar(w,inputvals, batch);
    % var=norm(bsxfun(@minus,g,mean), 'fro')^2;
    var = var / (n - 1);
    if nargout > 1
        skew = skew / n;
    end
end
