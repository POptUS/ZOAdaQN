function [iterations,funvals] = CallPlotExperiments_Cuter_Optimum(data, Options, Expmnt,i, flag1, flag2)
%This function loads the output files that gives the optimum function
%values for problems given in cuter dataset. 
if nargin<4
    i=1;
end

if nargin < 5
    flag1='No';
    flag2='No';
end


curr_path=pwd;
    resultDir = sprintf('%s/Results/%s/%s', curr_path,data,Expmnt);  
    resultString = sprintf('%s/%s_%s_%e_%s_%s_(%s_%f_%e_%e_%e_%s_%s_%i_%s)_(%d_%d_%d_%e_%e_%s_%s)_(%s_%d_%.2f_%.3f_%s)_%d',...
    resultDir,data,Options.Method,Options.lambda,Options.InitialWeights,Options.RandChoice,Options.LineSearch,Options.alpha0, Options.c, Options.rho, Options.tolr,...
    Options.Initial_SteplengthRule,Options.LineSearch_FailSkipRule,Options.LineSearch_Skipparam,Options.step_extra,... 
    Options.MaxIterations,Options.MaxEpochs,Options.MaxStorage, Options.epsilon,Options.StoreInterval,Options.StopTest,Options.Stop_extra,...
    Options.Curvature,Options.m,Options.threshold,Options.skip,Options.curvature_extra,...
    i);
    resultFile=strcat(resultString,'.mat');
%load(resultFile)
resultDir_short = sprintf('%s/Results', curr_path); 
resultString_short = sprintf('%s/%s_%s_%.2e_%d',resultDir_short,data,Options.Method,Options.lambda,Options.MaxEpochs);
resultFile_short=strcat(resultString_short,'.mat');
load(resultFile_short)



