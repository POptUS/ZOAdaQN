function [f, g] = calfun_Deterministic(w, inputvals)
    % This function computes the function and DFO gradient of the cuter
    % problems sets
    problem.m = inputvals.cuter_m;
    problem.nprob = inputvals.cuter_nprob;
    problem.probtype = 'smooth';
    sigma = inputvals.cuter_sigma;
    f = calfun_batch(w, sigma, 0, problem);

    if nargout > 1
        % Jacobian computation is not working; so will try dfo grad
        % [J,fvec]=jacobian(problem.m,length(w),w,problem.nprob);
        % g=sum(J,1)';

        % DFO grad
        options.interval = 10^-8;
        options.Method = 'FD';
        func = @(w)calfun_batch(w, sigma, 0, problem);
        [g, ~] = DFOGrad(w, func, options);
    end

end
