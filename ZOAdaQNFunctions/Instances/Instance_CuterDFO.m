% This file runs the FD-ASQN algorithms with the default parameters given
% in the paper.

% Options is the structue to provide the parameters required in the
% algorithm. Like what line search, sample selection,step lengths etc. As
% of now we are using the basic versions. But, later can incorporate all
% other options

%% Problem Parameters

% Options is the structure that consists of solution parameters
Options = Initialize(); % Initialized with empty/default values

% Sampling Parameters
Options.S0 = 2;   % initial sample size
Options.SamplingTest = 'Yes';
Options.Variance = 'Sample';
Options.gamma = 0.9;
Options.theta = 0.9;
Options.Maxvarbias = 5;   % Maximum sample increase
Options.SamVarSize = 1000; % Samples used in variance estimates
Options.Sample_decrease = [];

% Termination parameters
Options.epsilon = [];
Options.StopTest = [];

% Line Search Parameters
Options.LineSearch = 'Armijo';
Options.alpha0 = 1;
Options.c = 10^-4;
Options.rho = 1 / 2;
Options.tolr = 10^-8;
Options.Initial_SteplengthRule = 'Norm';
% Options.LineSearch_FailSkipRule='Con'; % constant
% Options.LineSearch_Skipparam=2;
Options.LineSearch_c2 = 10^-16;
Options.MaxL = 10^8;

% Curvature L-BFGS parameters
Options.Curvature = 'Yes';
Options.m = 10; % memory size
Options.curvskiprule = 'sp'; % 'wlf' for wolfe; 'sp' for y^Ts > epsilon \|s\|^2
Options.skip = 0.001;
Options.threshold = 0;
Options.ovp_type = 'FO';
Options.ovp_rule = 'Fd';
Options.ovp = 0;

% DFO parameters
datapoints = {100000};
Options.cuter_maxsample = 100000;
Options.MaxEpochs = 1;
Options.DFOMethod = 'FD';
Options.DFOMachinePrecision = 10^-16;
Options.DFOIntervalFactor = 1;
Options.StoreInterval = 1 / 10^12; % Frequency at which output stats should be stored

% Using default data, loss, sigma if nothing is specified.
if ~exist('datas', 'var')
    datas = {'15-absnormal', '15-relnormal'}; % Chebyquad
end
if ~exist('loss', 'var')
    loss = 'CuterDFO';
end
if ~exist('sigmas', 'var')
    sigmas = {10^-3};
end

if ~exist('rand_runs_adamethods', 'var')
    rand_runs_adamethods = 5;
end
% Loss
% loss='CuterDFO';%loss='CuterDFO-PFLoss';
Options.loss = loss;

% Folder Name
switch loss
    case 'CuterDFO'
        ExpmntFolder = 'Expmnt-Cuter';
    case 'MAD'
        ExpmntFolder = 'Expmnt-DFO';
    otherwise
        ExpmntFolder = 'Expmnt-Cuter';
end
% Running for different data and noise levels
for dat = 1:length(datas)
    data = datas{dat};
    % Running for each sigma (noise level)
    for ik = 1:length(sigmas)
        Options.lambda = sigmas{ik};

        % Norm Test;
        for seed = 1:rand_runs_adamethods
            Options.Method = 'Norm-AD1';
            rng(seed, 'twister');
            CallZOAdaQN(loss, data, Options, ExpmntFolder, seed);
        end

        % Inner Prodcut Quasi-Newton Test
        for seed = 1:rand_runs_adamethods
            Options.Method = 'IPQN-AD1';
            rng(seed, 'twister');
            CallZOAdaQN(loss, data, Options, ExpmntFolder, seed);
        end

    end
end
