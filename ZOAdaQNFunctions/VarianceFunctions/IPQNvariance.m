function [samvar] = IPQNvariance(batch, bias, Obj, inputvals, w, H0, history, g_s, k)
    % This function calculates the sample variance estimate of the
    % inner product quasi-Newton (ipqn) condition.
    %
    % The inputs are:
    % batch: indices in S_k^v used for sample estimates
    % mean: average of the sample estimates
    % Obj: This is declared as a function which gives the functional value and
    %      the gradient (estimated via finit-differences) value batch cases.
    % Obj.funcBatch: Gives the sample average objective and sample average
    %                gradient.
    % inputvals: This is the input that need to be given into the Obj for
    %            getting the functional values.
    % w: the current iterate
    % H0: initial estimate of the quasi-Newton hessian
    % history: consists of curvature pairs history.y, histiry.s
    % k: number of curvature pairs in history
    %
    % The outputs are:
    %  samvar: sample variance (2nd moment) of ipqn condition

    % Problem parameters
    n = length(batch);

    % initialization
    Hip = Twolooprecurrsion(H0, history, bias, k);
    nor = (Hip' * g_s);
    samvar = 0;
    for i = 1:n
        [~, g] = Obj.funcBatch(w, inputvals, batch(i));
        % d=Twolooprecurrsion(H0,history,g,k);
        ip = g' * Hip;
        samvar = samvar + (ip - nor)^2;
    end
    samvar = samvar / (n - 1);
end
