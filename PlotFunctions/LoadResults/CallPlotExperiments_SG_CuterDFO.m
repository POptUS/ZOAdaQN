function [iterations,batchsizes, funvals, Evals,FEvals,Stepsize,norminfgrad,norm2grad,Sample_Grad] = CallPlotExperiments_SG_CuterDFO(data, Options, Expmnt,i, flag1, flag2)
%This function loads the output files corresponding to stochastic gradient
%estimation runs, e.g., Finite-Difference Stochastic Gradient; 
%Sphere Smoothing

if nargin<4
    i=1;
end

if nargin < 5
    flag1='No';
    flag2='No';
end

Folder='Results';

curr_path=pwd;
resultDir = sprintf('%s/%s/%s/%s', curr_path,Folder,data,Expmnt);  
resultString = sprintf('%s/%s_%s_%e_%s_%s_(%s_%f_%e_%e_%e_%s_%s_%i_%s)_(%d_%d_%d_%e_%e_%s_%s)_(%d_%s_%.2f_%d_%s_%d_%d_%s_%s_%f)_(%s_%d_%.2f_%.3f_%s)_(%s_%s_%d)_(%s_%s_%s_%.2f_%s)_%d',...
resultDir,data,Options.Method,Options.lambda,Options.InitialWeights,Options.RandChoice,Options.LineSearch,Options.alpha0, Options.c, Options.rho, Options.tolr,...
Options.Initial_SteplengthRule,Options.LineSearch_FailSkipRule,Options.LineSearch_Skipparam,Options.step_extra,... 
Options.MaxIterations,Options.MaxEpochs,Options.MaxStorage, Options.epsilon,Options.StoreInterval,Options.StopTest,Options.Init_Stor,...
Options.S0, Options.SamplingTest,Options.theta,Options.gamma,Options.Variance,Options.SamVarSize,Options.Maxvarbias,Options.Sample_decrease,Options.Sample_extra,Options.GeometricRate,...
Options.Curvature,Options.m,Options.threshold,Options.skip,Options.curvature_extra,...
Options.ovp_type,Options.ovp_rule,Options.ovp,...
Options.DFOMethod,Options.DFOInterval,Options.DFOIntervalAdaptive, Options.DFOIntervalFactor, Options.DFODirs,...
i);
resultFile=strcat(resultString,'.mat');
%load(resultFile)
resultDir_short = sprintf('%s/Results', curr_path); 
resultString_short = sprintf('%s/%s_%s_%.2e_%i_%.2f_%s_%s_%d',resultDir_short,data,Options.Method,Options.lambda,log2(Options.alpha0),Options.theta,Options.DFOMethod, Options.DFODirs,i);
resultFile_short=strcat(resultString_short,'.mat');
load(resultFile_short)
% Setting the function evaluations to be zeros in constant stepsize sg
% methods
FEvals=zeros(length(iterations),1);
end



