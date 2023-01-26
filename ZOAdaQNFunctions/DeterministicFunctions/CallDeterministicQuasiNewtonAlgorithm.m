function [] = CallDeterministicQuasiNewtonAlgorithm(loss, data, Options, Expmnt, i)
    % This functions takes in the input values from the model and calls in the
    % respective optimization algorithm and stores the data in the respective
    % results file with the respective name.
    if nargin < 5
        i = 1;
    end
    if nargin < 4
        Expmnt = 'Optimum-Cuter';
    end

    % for i=1:Options.RandRuns
    curr_path = pwd;
    resultDir = sprintf('%s/Results/%s/%s', curr_path, data, Expmnt);
    resultString = sprintf(['%s/%s_%s_%e_%s_%s_(%s_%f_%e_%e_%e_%s_%s_%i_%s)_', ...
                   '(%d_%d_%d_%e_%e_%s_%s)_(%s_%d_%.2f_%.3f_%s)_%d'], ...
                   resultDir, data, Options.Method, Options.lambda, Options.InitialWeights, Options.RandChoice, ...
                   Options.LineSearch, Options.alpha0, Options.c, Options.rho, Options.tolr, ...
                   Options.Initial_SteplengthRule, Options.LineSearch_FailSkipRule, ...
                   Options.LineSearch_Skipparam, Options.step_extra, Options.MaxIterations, ...
                   Options.MaxEpochs, Options.MaxStorage, Options.epsilon, Options.StoreInterval, ...
                   Options.StopTest, Options.Stop_extra, Options.Curvature, Options.m, ...
                   Options.threshold, Options.skip, Options.curvature_extra, i);
    resultFile = strcat(resultString, '.mat');
    resultDir_short = sprintf('%s/Results', curr_path);
    resultString_short = sprintf('%s/%s_%s_%.2e_%d', resultDir_short, ...
    data, Options.Method, Options.lambda, Options.MaxEpochs);
    resultFile_short = strcat(resultString_short, '.mat');
    %     if ~exist(resultDir,'dir')
    %         mkdir(resultDir);
    %     end
    if ~exist(resultDir_short, 'dir')
        mkdir(resultDir_short);
    end

    if exist(resultFile_short, 'file')
        fprintf('\nResult file exists\n');
    else
        %diary(resultString_short);
        [Obj, inputvals, w] = SettingObjective(loss, data, Options);

        tic;
        % t=clock;
        [w, output] = DeterministicQuasiNewtonAlgorithms(Obj, inputvals, Options, w);
        % elaptime=etime(clock,t);
        TotalCPUtime = toc;
        % Performance Statistics
        iterations = output.iterations;
        funvals = output.funvals;
        norm2grad = output.norm2grad;
        norminfgrad = output.norminfgrad;
        Evals = output.Evals;
        FEvals = output.FEvals;
        time = output.time;
        testfunvals = output.testfunvals;
        corrclass = output.corrclass;

        % Line Search Statistics
        Stepsize = output.alpha;
        lsevals = output.lsEvals;
        zeta_vals = output.zetavals;
        lsFails = output.lsFails;

        Norm_Direction_d = output.normdirection;
        % Curvature Statistics
        curvskips = output.curvskips; % iterations where curvature updates are skipped
        curvature = output.curvature; % y^Ts;
        curvaturecomp = output.curvaturecomp; % s^TBs;

        save(resultFile_short, 'w', 'iterations', 'Evals', 'FEvals', 'funvals' ...
             , 'norm2grad', 'norminfgrad', 'time', 'testfunvals', 'corrclass' ...
             , 'Stepsize', 'lsevals', 'zeta_vals', 'lsFails' ...
             , 'Norm_Direction_d' ...
             , 'curvskips', 'curvature', 'curvaturecomp' ...
             , 'Options', '-v7.3');

        fprintf('Total Cpu time to run the method is %9.4e\n', TotalCPUtime);
        diary off;
    end

end
