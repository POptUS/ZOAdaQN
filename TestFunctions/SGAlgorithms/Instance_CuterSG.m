% This file runs the stochastic finite-difference gradient and sphere
% smoothing algorithms

clear;
clc;

% Options is the structue to provide the parameters required in the
% algorithm.
Options = Initialize();
Options.DFOdirs = [];

Options.S0 = 64;
Options.RandRuns = 1;

Options.LineSearch = 'Constant';

Options.MaxEpochs = 10;
Options.L0 = 1;

Options.epsilon = 10^-6;
Options.StopTest = [];

% Cuter DFO settings
% datas={'15-absnormal', '16-absnormal', '20-absnormal', '15-relnormal','15-absnormal-10'};
datas = {'15-absnormal', '20-absnormal', '15-relnormal', '22-absnormal', ...
    '19-absnormal-50', '18-absnormal', '20-relnormal', '22-relnormal', '19-relnormal-50', ...
    '18-relnormal', '205-absnormal', '205-relnormal', ...
    '220-absnormal', '220-relnormal', '216-absnormal', '216-relnormal', ...
    '305-absnormal', '305-relnormal', '124-absnormal', '124-relnormal', ...
    '123-absnormal', '123-relnormal'};
datapoints = {100000};
Options.cuter_maxsample = 100000;
Options.MaxEpochs = 1;
Options.DFOMethod = 'FD';
Options.DFOMachinePrecision = 10^-16;
Options.DFOInterval = 'Fixed';
Options.DFOIntervalFactor = 1;
Options.StoreInterval = 1 / 1000;
Options.Init_Stor = '1000';
loss = 'CuterDFO';
%% SG Algorithm - Running for sigma=10^-5;
sigma = 10^-5;
Options.S0 = 2;
% Optimal stepsize values for the entries in the datas
steps = {2^-10, 2^-12, 2^-10, 2^-11, 2^-13, 2^-5, 2^-11, ...
         2^-10, 2^-13, 2^-5, 2^-9, 2^-8, 2^-9, 2^-9, 2^-8, ...
         2^-8, 2^-10, 2^-10, 2^-6, 2^-6, 2^-20, 2^-20};
Options.DFODirs = [];
for dat = 1:22
    for seed = 1:5
        data = datas{dat};
        Options.lambda = sigma;
        Options.alpha0 = steps{dat};
        Options.Method = 'SG';
        rng(seed, 'twister');
        t = clock;
        CallSGAlgorithm(loss, data, Options, 'SG-DFO', seed);
    end
end

%% SG Algorithm - Running for sigma=10^-3;
sigma = 10^-3;
Options.S0 = 2;
% Optimal stepsize values
steps = {2^-10, 2^-12, 2^-10, 2^-11, 2^-13, 2^-5, 2^-10, 2^-10, ...
         2^-13, 2^-5, 2^-11, 2^-8, 2^-9, 2^-9, 2^-8, 2^-8, ...
         2^-10, 2^-10, 2^-6, 2^-6, 2^-20, 2^-20};
for dat = 1:22
    for seed = 1:5
        data = datas{dat};
        Options.lambda = sigma;
        Options.alpha0 = steps{dat};
        Options.Method = 'SG';
        rng(seed, 'twister');
        t = clock;
        CallSGAlgorithm(loss, data, Options, 'SG-DFO', seed);
    end
end

%% SG Algorithm - Running for Rand-50-50
loss = 'MAD';
Options.lambda = [];
data = 'Rand-50-50';

% % For identity matrix
% data='Eye-50-50';
% %Optimal stepsize
% Options.alpha0=2^-10; %Eye-50-50 alpha=2^-5 for 100; 2^-10 for 2

% Optimal stepsize
Options.alpha0 = 2^-14; % Rand-50-50 alpha=2^-12 for 100; 2^-14 for 2
Options.S0 = 2;
for seed = 1:5
    CallSGAlgorithm(loss, data, Options, 'SG-DFO', seed);
end

%% Running Sphere Smoothing Algorithms now
datas = {'15-absnormal', '20-absnormal', '15-relnormal', '22-absnormal', ...
    '19-absnormal-50', '18-absnormal', '20-relnormal', '22-relnormal', '19-relnormal-50', ...
    '18-relnormal', '205-absnormal', '205-relnormal', ...
    '220-absnormal', '220-relnormal', '216-absnormal', '216-relnormal', ...
    '305-absnormal', '305-relnormal', '124-absnormal', '124-relnormal', ...
    '123-absnormal', '123-relnormal'};
sigmas = {10^-3};
datapoints = {100000};
Options.cuter_maxsample = 100000;
Options.epsilon = 10^-6;
Options.StopTest = [];
loss = 'CuterDFO';
dats = 1:length(datas);
%% Sphere Smoothing Algorithm- Running for sigma=10^-5;
Options.DFOMethod = 'SS';
Options.DFODirs = '5'; % DFO directions
sigma = 10^-5;
Options.S0 = 2;
steps = {2^-13, 2^-12, 2^-13, 2^-12, 2^-14, 2^-8, 2^-12, ...
         2^-11, 2^-14, 2^-8, 2^-11, 2^-11, 2^-11, 2^-11, 2^-9, ...
         2^-9, 2^-15, 2^-15, 2^-13, 2^-13, 2^-24, 2^-24}; % for ss
% steps={2^-8, 2^-9, 2^-8, 2^-9, 2^-8,2^-4,2^-9,2^-9,2^-8,2^-5}; %for gs
% steps={2^-9, 2^-12, 2^-10, 2^-11, 2^-10,2^-5,2^-11,2^-10,2^-10,2^-5};%for sg
for dat = 1:22
    for seed = 1:5
        data = datas{dat};
        Options.lambda = sigma;
        Options.alpha0 = steps{dat};
        Options.Method = 'SG';
        rng(seed, 'twister');
        t = clock;
        CallSGAlgorithm(loss, data, Options, 'SG-DFO', seed);
    end
end

%% Sphere Smoothing Algorithm - Running for sigma=10^-3;
sigma = 10^-3;
Options.S0 = 2;
variables = {30, 20, 30, 8, 50, 11, 20, 8, 50, 11, 125, ...
             125, 110, 110, 100, 100, 100, 100, 100, 100, ...
             100, 100, 100, 100, 100, 100};
steps = {2^-15, 2^-12, 2^-14, 2^-12, 2^-14, 2^-8, 2^-12, 2^-11, ...
         2^-14, 2^-8, 2^-12, 2^-11, 2^-11, 2^-11, 2^-9, 2^-9, 2^-15, ...
         2^-15, 2^-13, 2^-13, 2^-24, 2^-24};   % For SS
% steps={2^-10, 2^-11, 2^-9, 2^-10,2^-8,2^-4,2^-8,2^-9,2^-8,2^-4};   % For GS
% steps={2^-10, 2^-12, 2^-10, 2^-11,2^-10,2^-5,2^-10,2^-10,2^-12,2^-5}; % for FD-SG
for dat = 1:22
    for seed = 1:5
        data = datas{dat};
        Options.lambda = sigma;
        Options.alpha0 = steps{dat};
        Options.Method = 'SG';
        rng(seed, 'twister');
        t = clock;
        CallSGAlgorithm(loss, data, Options, 'SG-DFO', seed);
    end
end

%% Sphere Smoothing Algorithm - Running for Rand-50-50
loss = 'MAD';
Options.lambda = [];
Options.S0 = 2;
% %For identity matrix data
% data='Eye-50-50';
% %For FD-SG
% Options.alpha0=2^-10; %Eye-50-50 alpha=2^-5 for 100; 2^-10 for 2
% %For GS
% Options.alpha0=2^-8; %Eye-50-50 alpha=2^-5 for 100; 2^-10 for 2

data = 'Rand-50-50';
% Options.alpha0=2^-16; %Rand-50-50 alpha=2^-12 for 100; 2^-14 for 2-FD-SG
Options.alpha0 = 2^-19; % Rand-50-50 alpha=2^-12 for 100; 2^-14 for 2
for seed = 1:5
    rng(seed, 'twister');
    CallSGAlgorithm(loss, data, Options, 'SG-DFO', seed);
end
