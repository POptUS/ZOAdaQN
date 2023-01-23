function [Obj, inputvals, w] = SettingObjective(loss, data, Options)
    % This function sets the objective functions for a given loss function is given.
    % Usually loss function is CuterDFO
    if nargin < 2 || isempty(data)
        data = '15-absnormal';
    end

    if nargin < 1 || isempty(loss)
        loss = 'CuterDFO';
    end

    if isempty(Options.lambda)
        Options.lambda = [];
    end

    % Printing the problem details
    fprintf(['loss = %s, data = %s, Lambda = %e, maxIter = %d, maxEpochs= %d\n method = %s', ...
            ', InitialSampleSize = %d, stepsizeselection = %s, InitialStepSize = %f\n,theta= %f, Variance = %s'], ...
            loss, data, Options.lambda, Options.MaxIterations, Options.MaxEpochs, ...
            Options.Method, Options.S0, Options.LineSearch, Options.alpha0, ...
            Options.theta, Options.Variance);

    fprintf('\n Loading the data set...');

    inputvals.loss = loss;

    % Depending on data set, we may require addtional adjustments like adding a
    % bias variable etc... The adjustments for the following known data sets
    % is done. If the data set is not specified here, we need to verify them
    % manually.

    % Parameters used in previous codes; ignore these for the current code
    inputvals.A = 0;
    inputvals.y = 0;
    inputvals.At = (inputvals.A)';
    inputvals.ndatapnts = size(inputvals.A, 1);
    inputvals.nvars = size(inputvals.A, 2);
    inputvals.lambda = Options.lambda;

    % Initialization
    if strcmp(Options.InitialWeights, 'Zeros') == 1
        w = zeros(inputvals.nvars, 1); % giving the intitial weights.
    elseif strcmp(Options.InitialWeights, 'Random') == 1
        w = zeros(inputvals.nvars, 1); % giving the intitial weights.
    else
        w = Options.InitialWeights;
    end

    % Setting the Objective Function
    fprintf('\n setting up the Objective function... ');
    switch loss
        % Cuter loss where the gradient are computed using the code provided in the software
        case 'Cuter'
            switch data
                case '20-absnormal'
                    inputvals.cuter_m = 30;
                    inputvals.cuter_n = 20;
                    inputvals.cuter_nprob = 20;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.nvars = 20;
                    w = dfoxs(20, 20, 1);
                case '15-relnormal'
                    inputvals.cuter_m = 45;
                    inputvals.cuter_n = 30;
                    inputvals.cuter_nprob = 15;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.nvars = 30;
                    w = dfoxs(30, 15, 1);
                case '15-absnormal'
                    inputvals.cuter_m = 45;
                    inputvals.cuter_n = 30;
                    inputvals.cuter_nprob = 15;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.nvars = 30;
                    w = dfoxs(30, 15, 1);
                case '16-absnormal'
                    inputvals.cuter_m = 45;
                    inputvals.cuter_n = 30;
                    inputvals.cuter_nprob = 16;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.nvars = 30;
                    w = dfoxs(30, 16, 1);
                case '21-absnormal-50'
                    inputvals.cuter_m = 50;
                    inputvals.cuter_n = 50;
                    inputvals.cuter_nprob = 21;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 21, 1);
                case '22-absnormal'
                    inputvals.cuter_m = 8;
                    inputvals.cuter_n = 8;
                    inputvals.cuter_nprob = 22;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 8;
                    w = dfoxs(8, 22, 1);
                case '19-absnormal-50'
                    inputvals.cuter_m = (50 - 4) * 2; % m=(n-4)*2
                    inputvals.cuter_n = 50; % n>=5;
                    inputvals.cuter_nprob = 19;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 19, 1);
                case '18-absnormal'
                    inputvals.cuter_m = 65;
                    inputvals.cuter_n = 11;
                    inputvals.cuter_nprob = 18;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 11;
                    w = dfoxs(11, 18, 1);
                case '20-relnormal'
                    inputvals.cuter_m = 30;
                    inputvals.cuter_n = 20;
                    inputvals.cuter_nprob = 20;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 20;
                    w = dfoxs(20, 20, 1);
                case '21-relnormal-50'
                    inputvals.cuter_m = 50;
                    inputvals.cuter_n = 50;
                    inputvals.cuter_nprob = 21;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 21, 1);
                case '22-relnormal'
                    inputvals.cuter_m = 8;
                    inputvals.cuter_n = 8;
                    inputvals.cuter_nprob = 22;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 8;
                    w = dfoxs(8, 22, 1);
                case '19-relnormal-50'
                    inputvals.cuter_m = (50 - 4) * 2; % m=(n-4)*2
                    inputvals.cuter_n = 50; % n>=5;
                    inputvals.cuter_nprob = 19;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 19, 1);
                case '18-relnormal'
                    inputvals.cuter_m = 65;
                    inputvals.cuter_n = 11;
                    inputvals.cuter_nprob = 18;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 11;
                    w = dfoxs(11, 18, 1);
                case '205-absnormal'
                    inputvals.cuter_m = 125;
                    inputvals.cuter_n = 125;
                    inputvals.cuter_nprob = 205;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 125;
                    w = dfoxsnew(125, 125, 205);
                case '205-relnormal'
                    inputvals.cuter_m = 125;
                    inputvals.cuter_n = 125;
                    inputvals.cuter_nprob = 205;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 125;
                    w = dfoxsnew(125, 125, 205);
                case '220-absnormal'
                    inputvals.cuter_m = 110;
                    inputvals.cuter_n = 110;
                    inputvals.cuter_nprob = 220;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 110;
                    w = dfoxsnew(110, 110, 220);
                case '220-relnormal'
                    inputvals.cuter_m = 110;
                    inputvals.cuter_n = 110;
                    inputvals.cuter_nprob = 220;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 110;
                    w = dfoxsnew(110, 110, 220);
                case '305-absnormal'
                    inputvals.cuter_m = 198;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 305;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(198, 100, 305);
                case '305-relnormal'
                    inputvals.cuter_m = 198;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 305;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(198, 100, 305);
                case '216-absnormal'
                    inputvals.cuter_m = 108;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 216;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(100, 100, 216);
                case '216-relnormal'
                    inputvals.cuter_m = 108;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 216;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(100, 100, 216);
                case '124-absnormal'
                    inputvals.cuter_m = 200;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 124;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(200, 100, 124);
                case '124-relnormal'
                    inputvals.cuter_m = 200;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 124;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(200, 100, 124);
                case '123-absnormal'
                    inputvals.cuter_m = 101;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 123;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(101, 100, 123);
                case '123-relnormal'
                    inputvals.cuter_m = 101;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 123;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(101, 100, 123);
            end
            Obj.func = @(w, inputvals)calfun_Deterministic(w, inputvals);
            Obj.test = @(inputvals, w)calfun_Test(w, inputvals); % Not used in the analysis

            % Cuter loss where the gradient are computed using function values
        case 'CuterDFO'
            inputvals.DFOMethod = Options.DFOMethod;
            inputvals.DFOMachinePrecision = Options.DFOMachinePrecision;
            inputvals.DFOInterval = Options.DFOInterval;
            inputvals.DFOIntervalAdaptive = Options.DFOIntervalAdaptive;
            inputvals.DFOIntervalFactor = Options.DFOIntervalFactor;
            switch data
                case '20-absnormal'
                    inputvals.cuter_m = 30;
                    inputvals.cuter_n = 20;
                    inputvals.cuter_nprob = 20;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 20;
                    w = dfoxs(20, 20, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '15-relnormal'
                    inputvals.cuter_m = 45;
                    inputvals.cuter_n = 30;
                    inputvals.cuter_nprob = 15;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 30;
                    w = dfoxs(30, 15, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '15-absnormal'
                    inputvals.cuter_m = 45;
                    inputvals.cuter_n = 30;
                    inputvals.cuter_nprob = 15;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 30;
                    w = dfoxs(30, 15, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '15-absnormal-10'
                    inputvals.cuter_m = 450;
                    inputvals.cuter_n = 300;
                    inputvals.cuter_nprob = 15;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 300;
                    w = dfoxs(300, 15, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '16-absnormal'
                    inputvals.cuter_m = 45;
                    inputvals.cuter_n = 30;
                    inputvals.cuter_nprob = 16;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 30;
                    w = dfoxs(30, 16, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '21-absnormal-50'
                    inputvals.cuter_m = 50;
                    inputvals.cuter_n = 50;
                    inputvals.cuter_nprob = 21;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 21, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '22-absnormal'
                    inputvals.cuter_m = 8;
                    inputvals.cuter_n = 8;
                    inputvals.cuter_nprob = 22;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 8;
                    w = dfoxs(8, 22, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '19-absnormal-50'
                    inputvals.cuter_m = (50 - 4) * 2; % m=(n-4)*2
                    inputvals.cuter_n = 50; % n>=5;
                    inputvals.cuter_nprob = 19;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 19, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '18-absnormal'
                    inputvals.cuter_m = 65;
                    inputvals.cuter_n = 11;
                    inputvals.cuter_nprob = 18;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 11;
                    w = dfoxs(11, 18, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '1-absnormal-100'
                    inputvals.cuter_m = 100;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 1;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxs(100, 1, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '20-relnormal'
                    inputvals.cuter_m = 30;
                    inputvals.cuter_n = 20;
                    inputvals.cuter_nprob = 20;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 20;
                    w = dfoxs(20, 20, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '21-relnormal-50'
                    inputvals.cuter_m = 50;
                    inputvals.cuter_n = 50;
                    inputvals.cuter_nprob = 21;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 21, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '22-relnormal'
                    inputvals.cuter_m = 8;
                    inputvals.cuter_n = 8;
                    inputvals.cuter_nprob = 22;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 8;
                    w = dfoxs(8, 22, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '19-relnormal-50'
                    inputvals.cuter_m = (50 - 4) * 2; % m=(n-4)*2
                    inputvals.cuter_n = 50; % n>=5;
                    inputvals.cuter_nprob = 19;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 50;
                    w = dfoxs(50, 19, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '18-relnormal'
                    inputvals.cuter_m = 65;
                    inputvals.cuter_n = 11;
                    inputvals.cuter_nprob = 18;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 11;
                    w = dfoxs(11, 18, 1);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '205-absnormal' % including new datasets after first revision
                    inputvals.cuter_m = 125;
                    inputvals.cuter_n = 125;
                    inputvals.cuter_nprob = 205;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 125;
                    w = dfoxsnew(125, 125, 205);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '205-relnormal' % including new datasets after first revision
                    inputvals.cuter_m = 125;
                    inputvals.cuter_n = 125;
                    inputvals.cuter_nprob = 205;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 125;
                    w = dfoxsnew(125, 125, 205);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '220-absnormal' % including new datasets after first revision
                    inputvals.cuter_m = 110;
                    inputvals.cuter_n = 110;
                    inputvals.cuter_nprob = 220;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 110;
                    w = dfoxsnew(110, 110, 220);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '220-relnormal' % including new datasets after first revision
                    inputvals.cuter_m = 110;
                    inputvals.cuter_n = 110;
                    inputvals.cuter_nprob = 220;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 110;
                    w = dfoxsnew(110, 110, 220);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '305-absnormal' % including new datasets after first revision
                    inputvals.cuter_m = 198;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 305;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(198, 100, 305);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '305-relnormal' % including new datasets after first revision
                    inputvals.cuter_m = 198;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 305;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(198, 100, 305);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '216-absnormal' % including new datasets after first revision
                    inputvals.cuter_m = 100;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 216;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(100, 100, 216);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '216-relnormal' % including new datasets after first revision
                    inputvals.cuter_m = 100;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 216;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(100, 100, 216);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '124-absnormal' % including new datasets after first revision
                    inputvals.cuter_m = 200;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 124;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(200, 100, 124);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '124-relnormal' % including new datasets after first revision
                    inputvals.cuter_m = 200;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 124;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(200, 100, 124);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '123-absnormal' % including new datasets after first revision
                    inputvals.cuter_m = 101;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 123;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'absnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(101, 100, 123);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
                case '123-relnormal' % including new datasets after first revision
                    inputvals.cuter_m = 101;
                    inputvals.cuter_n = 100;
                    inputvals.cuter_nprob = 123;
                    inputvals.cuter_sigma = Options.lambda;
                    inputvals.cuter_probtype = 'relnormal';
                    inputvals.ndatapnts = Options.cuter_maxsample;
                    inputvals.nvars = 100;
                    w = dfoxsnew(101, 100, 123);
                    A = randn(inputvals.ndatapnts, inputvals.cuter_m); % Random matrix consisting of 100000by m
                    inputvals.A = A; % Random matrix consisting of ndatapntsby m
                    Obj.func = @(w, inputvals)calfun_DFO(w, inputvals);
                    Obj.funcBatch = @(w, inputvals, batch)calfun_batchDFO_PFLoss(w, inputvals, batch);
                    Obj.test = @(inputvals, w)calfun_Test(w, inputvals);
            end

        case 'MAD' % Nonsmooth loss
            inputvals.ndatapnts = 100000;
            switch data
                case 'Eye-50-50'
                    m = 50;
                    d = 50;
                    inputvals.A = eye(m);
                    inputvals.xopt = randn(d, 1);
                    inputvals.b = inputvals.A * inputvals.xopt;
                    inputvals.randmat = rand(inputvals.ndatapnts, m);
                    w = zeros(d, 1);
                    inputvals.Zeta = 'Uni';
                    inputvals.U = 1;
                    inputvals.L = -1;
                case 'Rand-50-50'
                    m = 50;
                    d = 50;
                    A = randn(m, m);
                    A = 0.5 * (A + A');
                    A = A + m * eye(m);
                    inputvals.A = A;
                    inputvals.xopt = randn(d, 1);
                    inputvals.b = inputvals.A * inputvals.xopt;
                    inputvals.randmat = rand(inputvals.ndatapnts, m);
                    w = zeros(d, 1);
                    inputvals.Zeta = 'Uni';
                    inputvals.U = 1;
                    inputvals.L = -1;
            end
            inputvals.DFOMethod = Options.DFOMethod;
            % inputvals.DFOinterval=Options.DFOinterval;
            inputvals.DFOMachinePrecision = Options.DFOMachinePrecision;
            inputvals.DFOInterval = Options.DFOInterval;
            inputvals.DFOIntervalAdaptive = Options.DFOIntervalAdaptive;
            inputvals.DFOIntervalFactor = Options.DFOIntervalFactor;
            inputvals.data = data;
            Obj.func = @(w, inputvals)MAD_Expec(w, inputvals);
            Obj.funcBatch = @(w, inputvals, batch)MAD_BatchDFO_PFLoss(w, inputvals, batch);
            % test_seeds=randi([10^6 10^7],35,1);
            Obj.test = @(inputvals, w)(test(inputvals, w));

    end
    fprintf('\nDone\n');
end
