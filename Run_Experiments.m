% This is the main file to run the experiments and plot the corresponding
% results given in the paper "Adaptive Sampling Quasi-Newton Methods for
% Zeroth-Order Stochastic Optimization."
% The results are stored in the "Results" folder and plots are stored in
% the "Plots" folder.
%
% We run the expiremts for finite-difference adaptive sampling quasi-Newton
% methods. Base on the user input, we re-run the stochastic gradient based
% experiments for plotting. If the user input is "No" then we will only 
% plot finite-difference adaptive sampling quasi-Newton methods. 
%
% The optimal function values for each dataset are obtained by running
% "Instance_CuterOptimum.m" file in ZOAdaQNFunctions folder. The
% the working directory should be Matlab Code
%
% If one wants to re-run all the stochastic gradient experiments
% then please run the Instance_CuterSG.m file in
% ZOAdaQNFunctions/SGAlgorithms/. The working directory should still be
% Matlab Code. 
% 
% Please note that the optimal step-size for the datasets and
% the noise values (10^-3, 10^-5) have been found by grid search already. 
% If one wants to test noise parameter values other than those given in the 
% paper, please provide the corresponding optimal step-size value in the 
% Instance_CuterSG.m and in the PlotExperiments_DFO.m file
% 

% Please use the names under absloss and relloss in datas variable
% for running the experiments
% Dataset names look-up table:
%   Name        Number  absloss             relloss
%   ChebyQuad   15      15-absnormal        15-relnormal
%   Osborne     18      18-absnormal        18-relnormal
%   Bqdrtic     19      19-absnormal-50     19-relnormal-50
%   Cube        20      20-absnormal        20-relnormal
%   Heart8ls    22      22-absnormal        22-relnormal
%   BRATU3D     205     205-absnormal       205-relnormal
%   EIGENC      220     220-absnormal       220-relnormal
%   ConnBand    216     216-absnormal       216-relnormal
%   ROSENBR     305     305-absnormal       305-relnormal
%   PENALTY2    124     124-absnormal       124-relnormal
%   PENLT1NE    123     123-absnormal       123-relnormal

%% setting up the path
clear all;
clc;

% Adding all the folders in the path.
addpath(genpath(pwd));

%% Running experiments for any dataset (CUTER DFO)
datas = {'15-absnormal'}; % other values:15-relnormal, 18-absnormal, 18-relnormal,19-absnormal-50,
% 19-relnormal-50, 20-absnormal, 20-relnormal, 22-absnormal, 22-relnormal,124-absnormal,
% 124-relnormal,123-absnormal,123-relnormal,305-absnormal, 305-relnormal,
% 216-absnormal, 216-relnormal,205-absnormal,205-relnormal 220-absnormal, 220-relnormal
% datas={'15-absnormal','15-relnormal', '18-absnormal', '18-relnormal','19-absnormal-50',...
% '19-relnormal-50', '20-absnormal', '20-relnormal', '22-absnormal', '22-relnormal'
% '124-absnormal','124-relnormal','123-absnormal','123-relnormal'
% '305-absnormal', 305-relnormal,'305-absnormal', 305-relnormal
% '216-absnormal', '216-relnormal','205-absnormal', '205-relnormal'
% '220-absnormal', '220-relnormal'};
sigmas = {10^-3}; % other values: 10^-5;
loss = 'CuterDFO';
run_SG = input(strcat('Do you want to run the Stochastic Gradient experiments for the dataset \n', ...
    '   with datas = ', datas{:}, ', sigmas = ', num2str(sigmas{:}), ' \n', ...
    '(y/n)? \n', ...
    'If you do not select y, then only Adaptive Sampling Finite-Difference Quasi-Newton methods will be used to generate plots'), 's');
% select number of random runs;
rand_runs_adamethods = 5; % default 5
Instance_CuterDFO;             % Run ZOAdaQN
Instance_CuterOptimum_Single;  % Run Deterministic QN

if strcmp(run_SG, 'y') == 1
    Instance_CuterSG_Single;       % Run SG
    PlotExperiments_DFO            % Plotting SG and Adaptive Sampling Methods
else
    PlotExperiments_DFO_Ada;  % Plotting only Adaptive Sampling Methods       
end


%% Running experiments on nonsmooth loss (MAD loss) for random matrix
datas = {'Rand-50-50'};
sigmas = {[]};
loss = 'MAD';

Instance_CuterDFO;             % Run ZOAdaQN

if strcmp(run_SG, 'y') == 1
    Instance_CuterSG_Single;       % Run SG
    PlotExperiments_DFO            % Plotting SG and Adaptive Sampling Methods
else
    PlotExperiments_DFO_Ada;  % Plotting only Adaptive Sampling Methods       
end

