function [f, g] = calfun_batchDFO_PFLoss(w, inputvals, batch)
    % This function computes the function and DFO gradient of the cuter
    % problems sets
    problem.m = inputvals.cuter_m;
    problem.nprob = inputvals.cuter_nprob;
    problem.probtype = inputvals.cuter_probtype;
    sigma = inputvals.cuter_sigma;
    randmat = inputvals.A;
    randmat = randmat(batch, :);
    [f, fsmooth] = calfun_randmat(w, sigma, randmat, problem);
    func = @(w)(calfun_randmat(w, sigma, randmat, problem));
    DFOoptions.Method = inputvals.DFOMethod;
    DFOoptions.DFOsigma = sigma;
    DFOoptions.DFOproblem = problem;
    DFOoptions.DFOfsmooth = fsmooth;
    h = Cal_DFOinterval(inputvals, DFOoptions);
    DFOoptions.interval = h;
    if strcmp(DFOoptions.Method, 'GS') == 1 || strcmp(DFOoptions.Method, 'SS') == 1
        DFOoptions.dirs = inputvals.DFODirs;
    end
    % DFOoptions.interval=inputvals.DFOinterval;
    g = DFOGrad(w, func, DFOoptions);
end

% [J,~]=jacobian(problem.m,length(w),w,problem.nprob);
% gtrue=sum(J,1)';
% intervals=10.^[-12:0.1:log10(h)+0.1];
% errors=zeros(length(intervals),1);
% for j=1:length(intervals)
%  DFOoptions.interval=intervals(j);
%  gtest=DFOGrad(w, func, DFOoptions);
%  errors(j)=norm(gtrue - gtest);
% end
% semilogx(intervals, errors)
% hold on
% semilogx(h,norm(gtrue - g),'r*')
% title('15 - relnormal')
% xlabel('h-interval')
% ylabel('error')

% [J,~]=jacobian(problem.m,length(w),w,problem.nprob);
% gtrue=sum(J,1)';
% intervals=10^-12:10^-6:h;
% errors=zeros(length(intervals),1);
% for j=1:length(intervals)
%     DFOoptions.interval=intervals(j);
%     gtest=DFOGrad(w, func, DFOoptions);
%     errors(j)=norm(gtrue - gtest);
% end
% semilogy(intervals, errors)
% hold on
% semilogy(h,norm(gtrue - g),'r*')
% title('15 - relnormal')
% xlabel('h-interval')
% ylabel('error')
