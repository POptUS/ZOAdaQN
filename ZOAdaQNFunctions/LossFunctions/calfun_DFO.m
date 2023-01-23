function [f, g] = calfun_DFO(w, inputvals)
    % This function computes the function and DFO gradient of the cuter
    % problems sets
    problem.m = inputvals.cuter_m;
    problem.nprob = inputvals.cuter_nprob;
    problem.probtype = 'smooth';
    sigma = inputvals.cuter_sigma;
    % f=calfun_batch(w,sigma,0,problem); %Modifying this to be consistent
    % func= @(w)(calfun_batch(w,sigma,0,problem));
    f = calfun_randmat(w, sigma, 0, problem);
    func = @(w)(calfun_randmat(w, sigma, 0, problem));
    DFOoptions.Method = inputvals.DFOMethod;
    DFOoptions.DFOsigma = sigma;
    DFOoptions.DFOproblem = problem;
    DFOoptions.DFOfsmooth = f;
    h = Cal_DFOinterval(inputvals, DFOoptions);
    DFOoptions.interval = h;
    if strcmp(DFOoptions.Method, 'GS') == 1 || strcmp(DFOoptions.Method, 'SS') == 1
        DFOoptions.dirs = inputvals.DFODirs;
    end
    % DFOoptions.interval=inputvals.DFOinterval;
    g = DFOGrad(w, func, DFOoptions);
end
