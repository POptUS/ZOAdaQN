% This file runs the Deterministic algorithms to find the optimal solution
% for Cuter problems.

clear;
clc;

% Options is the structue to provide the parameters required in the
% algorithm.

%% Problem Parameters

Options = Initialize();

Options.RandRuns = 1;
Options.RandChoice = 'WOR';
Options.epsilon = 10^-8;
Options.StopTest = 'InftyNorm';

Options.LineSearch = 'DeterministicArmijo';
Options.alpha0 = 1;
Options.c = 10^-4;
Options.rho = 1 / 2;
Options.tolr = 10^-8;
Options.Initial_SteplengthRule = 'DSS'; % Inactive
Options.LineSearch_FailSkipRule = 'Con'; % constant
Options.LineSearch_Skipparam = 2;

Options.Curvature = 'Yes';
Options.m = 10;
Options.skip = 0.001;
Options.threshold = 0;

%% Cuter DFO settings
datas = {'15-absnormal', '20-absnormal', '15-relnormal', '22-absnormal', ...
       '19-absnormal-50', '18-absnormal', '20-relnormal', '22-relnormal', '19-relnormal-50', ...
       '18-relnormal', '205-absnormal', '216-absnormal', '220-absnormal'...
       '124-absnormal', '123-absnormal', '305-absnormal', ...
       '205-relnormal', '216-relnormal', '220-relnormal', ...
       '124-relnormal', '123-relnormal', '305-relnormal'};
sigmas = {10^-3, 10^-5};
datapoints = {100000};
Options.cuter_maxsample = 100000;
loss = 'Cuter';
dats = 1:length(datas);
Options.MaxEpochs = 2000;
% Options.DFOmethod='FD';
% Options.DFOinterval=sqrt(10^-16);
Options.DFOMethod = 'FD';
Options.DFOMachinePrecision = 10^-16;
Options.DFOIntervalFactor = 1;
Options.StoreInterval = 1 / 10^12;

for dat = 1:length(dats)
    data = datas{dats(dat)};
    for ik = 1:1
        Options.lambda = sigmas{ik};
        for seed = 1:1
            Options.Method = 'Deterministic';
            CallDeterministicQuasiNewtonAlgorithm(loss, data, Options, 'Optimum-Cuter', seed);
        end
    end
end
