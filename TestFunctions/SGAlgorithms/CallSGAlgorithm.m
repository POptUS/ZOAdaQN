function [] = CallSGAlgorithm(loss, data, Options, Expmnt, i)
% This functions takes in the input values from the model and calls in the
% respective optimization algorithm and stores the data in the respective
% results file with the respective name.
if nargin < 5
    i = 1;
end
if nargin < 4
    Expmnt = 'SG-DFO';
end

% for i=1:Options.RandRuns
    curr_path = pwd;
    resultDir = sprintf('%s/Results/%s/%s', curr_path, data, Expmnt);
    resultString = sprintf(['%s/%s_%s_%e_%s_%s_(%s_%f_%e_%e_%e_%s_%s_%i_%s)_', ...
         '(%d_%d_%d_%e_%e_%s_%s)_(%d_%s_%.2f_%d_%s_%d_%d_%s_%s_%f)_', ...
         '(%s_%d_%.2f_%.3f_%s)_(%s_%s_%d)_(%s_%s_%s_%.2f_%s)_%d'], ...
    resultDir, data, Options.Method, Options.lambda, Options.InitialWeights, ...
    Options.RandChoice, Options.LineSearch, Options.alpha0, Options.c, Options.rho, ...
    Options.tolr, Options.Initial_SteplengthRule, Options.LineSearch_FailSkipRule, ...
    Options.LineSearch_Skipparam, Options.step_extra, ...
    Options.MaxIterations, Options.MaxEpochs, Options.MaxStorage, ...
    Options.epsilon, Options.StoreInterval, Options.StopTest, Options.Init_Stor, ...
    Options.S0, Options.SamplingTest, Options.theta, Options.gamma, ...
    Options.Variance, Options.SamVarSize, Options.Maxvarbias, ...
    Options.Sample_decrease, Options.Sample_extra, Options.GeometricRate, ...
    Options.Curvature, Options.m, Options.threshold, Options.skip, Options.curvature_extra, ...
    Options.ovp_type, Options.ovp_rule, Options.ovp, ...
    Options.DFOMethod, Options.DFOInterval, Options.DFOIntervalAdaptive, ...
    Options.DFOIntervalFactor, Options.DFODirs, i);
    resultFile = strcat(resultString, '.mat');

    resultDir_short = sprintf('%s/Results', curr_path);
    resultString_short = sprintf('%s/%s_%s_%.2e_%i_%.2f_%s_%s_%d', ...
          resultDir_short, data, Options.Method, Options.lambda, log2(Options.alpha0), ...
          Options.theta, Options.DFOMethod, Options.DFODirs, i);
    resultFile_short = strcat(resultString_short, '.mat');

    % if ~exist(resultDir,'dir')
    %    mkdir(resultDir);
    % end
    if ~exist(resultDir_short, 'dir')
        mkdir(resultDir_short);
    end

    if exist(resultFile_short, 'file')
        fprintf('\nResult file exists\n');
    [Obj, inputvals, w] = SettingObjective(loss, data, Options);

    % including this for gaussian smoothing - DFO
    if isempty(Options.DFODirs) == 0
        inputvals.DFODirs = floor(str2num(Options.DFODirs));
    else
        inputvals.DFODirs = 0;
    end
    tic;
    % t=clock;
    [w, output] = SGAlgorithms(Obj, inputvals, Options, w);
    % elaptime=etime(clock,t);
    TotalCPUtime = toc;
    iterations = output.iterations;
    batchsizes = output.batchsizes;
    funvals = output.funvals;
    norm2grad = output.norm2grad;
    norminfgrad = output.norminfgrad;
    Evals = output.Evals;
    Sample_Obj = output.sampleobj;
    Sample_Grad = output.samplegrad;
    Stepsize = output.alpha;
    Descent = output.descent;
    time = output.time;
    testfunvals = output.testfunvals;
    corrclass = output.corrclass;
    lsevals = output.lsEvals;
    save(resultFile_short, 'w', 'batchsizes', 'Evals', 'funvals', 'norm2grad' ...
     , 'norminfgrad', 'Sample_Obj', 'Sample_Grad', 'Stepsize', 'Descent' ...
     , 'time', 'testfunvals', 'corrclass', 'Options', 'lsevals', 'iterations');
    fprintf('Total Cpu time to run the method is %9.4e\n', TotalCPUtime);
    diary off;
    end
end
