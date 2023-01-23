function [f, correct] = calfun_Test(w, inputvals)
    % This function computes the function and DFO gradient of the cuter
    % problems sets
    problem.m = inputvals.cuter_m;
    problem.nprob = inputvals.cuter_nprob;
    problem.probtype = 'smooth';
    sigma = inputvals.cuter_sigma;
    % f=calfun_batch(w,sigma,0,problem);  %Modifying this to be consistent
    f = calfun_randmat(w, sigma, 0, problem);
    correct = 0;
    % g=DFOGrad(w, func, DFOoptions);
end
