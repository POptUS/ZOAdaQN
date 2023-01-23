% This is the main file to run the experiments and plot the corresponding
% results given in the paper "Adaptive Sampling Quasi-Newton Methods for 
% Zeroth-Order Stochastic Optimization."
% The results are stored in the "Results" folder and plots are stored in 
% the "Plots" folder. 
% Note that, by default, all the results files exist in the Results 
% folder. Based on the user input, we re-run the experiments again
% for the finite-difference adaptive sampling quasi-Newton methods. We use
% the existing results of stochastic gradient based methods for plotting
% the results. 
%
% If one wants to re-run the stochastic gradient experiments 
% too, then please run the Instance_CuterSG.m file in 
% ZOAdaQNFunctions/SGAlgorithms/. The working directory should still be
% Matlab Code.
%
% Also note the optimum function values in each case are already computed 
% and stored in the results folder. So, if one wants to use noise
% parameter values other than those given in the paper, please run 
% the Instance_CuterOptimum.m file in ZOAdaQNFunctions folder.  Again, 
% the working directory should be Matlab Code 

% Please use the names under absloss and relloss in datas variable 
% for running the experiments 
% Dataset names look-up table:
%   Name        Number  absloss             relloss
%   ChebyQuad	15      15-absnormal        15-relnormal
%   Osborne     18      18-absnormal        18-relnormal
%   Bqdrtic     19      19-absnormal-50     19-relnormal-50
%   Cube        20      20-absnormal        20-relnormal
%   Heart8ls	22      22-absnormal        22-relnormal
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
datas={'15-absnormal'}; %other values:15-relnormal, 18-absnormal, 18-relnormal,19-absnormal-50,
% 19-relnormal-50, 20-absnormal, 20-relnormal, 22-absnormal, 22-relnormal,124-absnormal,
% 124-relnormal,123-absnormal,123-relnormal,305-absnormal, 305-relnormal,
% 216-absnormal, 216-relnormal,205-absnormal,205-relnormal 220-absnormal, 220-relnormal
%datas={'15-absnormal','15-relnormal', '18-absnormal', '18-relnormal','19-absnormal-50',...
%'19-relnormal-50', '20-absnormal', '20-relnormal', '22-absnormal', '22-relnormal'
% '124-absnormal','124-relnormal','123-absnormal','123-relnormal'
% '305-absnormal', 305-relnormal,'305-absnormal', 305-relnormal
% '216-absnormal', '216-relnormal','205-absnormal', '205-relnormal'
% '220-absnormal', '220-relnormal'};
sigmas={10^-3}; %other values: 10^-5;
loss='CuterDFO';
run=input(strcat('Do you want to re-run the experiments for the dataset \n', ...
    '   with datas = ',datas{:},', sigmas = ',num2str(sigmas{:}), ' \n', ...
    '(y/n)? \n', ...
    'If you do not select y, existing results will be used to generate plots'),'s');
if strcmp(run,'y')==1
    % select number of random runs;
    rand_runs_adamethods=5; % default 5
    Instance_CuterDFO; 		% Run ZOAdaQN
    Instance_CuterSG; 		% Run SG
    Instance_CuterOptimum; 	% Run Deterministic QN
end
PlotExperiments_DFO;

%% Running experiments on nonsmooth loss (MAD loss) for random matrix
datas={'Rand-50-50'};
sigmas={[]};
loss='MAD';
run=input(strcat('Do you want to re-run the experiments for the dataset \n', ...
    '   with datas = ',datas{:},', sigmas = ',num2str(sigmas{:}), ' \n', ...
    '(y/n)? \n', ...
    'If you do not select y, existing results will be used to generate plots'),'s');
if strcmp(run,'y')==1
    % select number of random runs;
    rand_runs_adamethods=5; % default 5
    Instance_CuterDFO; 
end
PlotExperiments_DFO;

