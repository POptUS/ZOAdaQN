function [] = CallZOAdaQN(loss, data, Options, Expmnt, i)
    % For license, source, & updates to ZOAdaQN, see https://github.com/POptUS/ZOAdaQN
    %
    % This functions takes in the input values from the model and calls in the
    % respective optimization algorithm and stores the data in the respective
    % results file with the respective name.
    %  Inputs:
    %          loss: (string)'DFO'
    %          data: (string) name of the data file
    %          Options:(class)Options class in the initialize file
    %          Expmnt: (string) name of the experiment
    %               i: random run number

    % Setting defaults
    if nargin < 5   % Default random run
        i = 1;
    end
    if nargin < 4   % Default results folder
        Expmnt = 'Expmnt-Cuter';
    end

    % Creating the result directory and result file
    curr_path = pwd;
    resultDir = sprintf('%s/Results/%s/%s', curr_path, data, Expmnt);
    resultDir_short = sprintf('%s/Results', curr_path);
    resultString_short = sprintf('%s/%s_%s_%.2e_%i_%.2f_%s_%d', resultDir_short, data, ...
         Options.Method, Options.lambda, log2(Options.alpha0), Options.theta, Options.DFOMethod, i);
    resultFile_short = strcat(resultString_short, '.mat');
    if ~exist(resultDir_short, 'dir')
        mkdir(resultDir_short);
    end

    if exist(resultFile_short, 'file')
        fprintf('\nResult file exists\n');
        fprintf('\nRunning the experiment again\n');
    end
    response = 'Y';
    if strcmp(response, 'Y') == 1
        % diary(resultString) % uncomment if you want to store the results in a diary
        [Obj, inputvals, w] = SettingObjective(loss, data, Options);

        tic;
        % t=clock;
        % running the algorithm
        [w, output] = ZOAdaQN(Obj, inputvals, Options, w);
        % elaptime=etime(clock,t);
        TotalCPUtime = toc;

        % Extracting the output
        % Performance Statistics
        iterations = output.iterations;
        funvals = output.funvals;
        norm2grad = output.norm2grad;
        norminfgrad = output.norminfgrad;
        Evals = output.Evals;
        FEvals = output.FEvals;
        time = output.time;
        testfunvals = output.testfunvals; % not for CUTER problem sets
        corrclass = output.corrclass;     % not for CUTER problem sets

        % Line Search Statistics
        Stepsize = output.alpha;
        lsevals = output.lsEvals;
        zeta_vals = output.zetavals;
        lsFails = output.lsFails;

        % Sample Statistics
        batchsizes = output.batchsizes;
        Sample_Obj = output.sampleobj;
        Sample_Grad = output.samplegrad;
        Sample_Var = output.samplevar;
        True_Var = output.truevar;
        varratios = output.varratios;
        Descent = output.descent;
        Norm_Direction_d = output.normdirection;
        DD_gs_g = output.gradtimessamgrad;
        delta_vals = output.deltavals;
        % Curvature Statistics
        curvskips = output.curvskips; % iterations where curvature updates are skipped
        curvature = output.curvature; % y^Ts;
        curvaturecomp = output.curvaturecomp; % s^TBs;

        % Overlap statistics
        S_ovp = output.s_ovp;
        S_bias = output.s_bias;
        S_waste = output.s_waste;

        % Storing the outputs.
        save(resultFile_short, 'w', 'iterations', 'Evals', 'FEvals', 'funvals' ...
             , 'norm2grad', 'norminfgrad', 'time', 'testfunvals', 'corrclass' ...
             , 'Stepsize', 'lsevals', 'zeta_vals', 'lsFails' ...
             , 'batchsizes', 'Sample_Obj', 'Sample_Grad', 'Sample_Var', 'True_Var' ...
             , 'varratios', 'Descent', 'Norm_Direction_d', 'DD_gs_g', 'delta_vals' ...
             , 'curvskips', 'curvature', 'curvaturecomp', 'S_ovp', 'S_bias', 'S_waste' ...
             , 'Options', '-v7.3');

        fprintf('Total Cpu time to run the method is %9.4e\n', TotalCPUtime);
        % diary off
    end
end
