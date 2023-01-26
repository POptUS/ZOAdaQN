function [] = SetupExperiment_DFO_Ada(data, lambda, n, theta, S0_SG, S0_Ada, xlimit_FEvals, xlimitwide, gamma, method, rand_runs_adamethods)
% This function sets the parameters to load all the results files
% and generates the plots corresponding to these results.

% Setting default values
if nargin < 4
    S0_SG = 2;
    S0_Ada = 2;
    xlimit_FEvals = 1000000;
    xlimitwide = 10000;
end

d = n; % number of parameters
%% Loading the Optimum Values obtained via L-BFGS
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
Options.Initial_SteplengthRule = 'DSS';
Options.LineSearch_FailSkipRule = 'Con'; % constant
Options.LineSearch_Skipparam = 2;

Options.Curvature = 'Yes';
Options.m = 10;
Options.skip = 0.001;
Options.threshold = 0;

Options.cuter_maxsample = 100000;
Options.MaxEpochs = 2000;
Options.DFOmethod = 'FD';
Options.DFOinterval = sqrt(10^-16);
Options.StoreInterval = 1 / 10^12;
Options.lambda = lambda; % Changing to lambda instead of 10-3 for rerunningthe optimal value scenario. 
Options.Method = 'Deterministic';

switch data
    case {'Rand-50-50', 'Eye-50-50'} % For MAD rand data optimal value is known
        Opt_Val = 25;
    otherwise
        [~, funvals] = CallPlotExperiments_Cuter_Optimum(data, Options, 'Optimum-Cuter', 1);
        Opt_Val = min(funvals);
end

%% Loading the results of FD-ASQN methods
Options = Initialize(); % Initialized with empty/default values

% Sampling Parameters
Options.S0 = S0_Ada;   % initial sample size
Options.SamplingTest = 'Yes';
Options.Variance = 'Sample';
Options.gamma = gamma;
Options.theta = theta;
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
Options.DFOMachinePrecision = 10^-14;
Options.DFOIntervalFactor = 1;
Options.StoreInterval = 1 / 10^12; % Frequency at which output stats should be stored

% Loss
loss = 'CuterDFO'; % loss='CuterDFO-PFLoss';
Options.loss = loss;
Options.lambda = lambda;

Options.Method = "Norm" + method;
Options.DFOMethod = 'FD';
switch data
    case {'Rand-50-50', 'Eye-50-50'}
        [iter3, batch3, funvalsiter3, funvalsepoch3, epoch3, steps3, ExpBudgetMarkers3] = Output_Stats(@CallPlotExperiments_Cuter, data, Options, 'Expmnt-DFO', [1:rand_runs_adamethods], d, d + 1);
    otherwise
        [iter3, batch3, funvalsiter3, funvalsepoch3, epoch3, steps3, ExpBudgetMarkers3] = Output_Stats(@CallPlotExperiments_Cuter, data, Options, 'Expmnt-Cuter', [1:rand_runs_adamethods], d, d + 1);
end

Options.Method = "IPQN" + method;
Options.DFOMethod = 'FD';
switch data
    case {'Rand-50-50', 'Eye-50-50'}
        [iter4, batch4, funvalsiter4, funvalsepoch4, epoch4, steps4, ExpBudgetMarkers4] = Output_Stats(@CallPlotExperiments_Cuter, data, Options, 'Expmnt-DFO', [1:rand_runs_adamethods], d, d + 1);
    otherwise
        [iter4, batch4, funvalsiter4, funvalsepoch4, epoch4, steps4, ExpBudgetMarkers4] = Output_Stats(@CallPlotExperiments_Cuter, data, Options, 'Expmnt-Cuter', [1:rand_runs_adamethods], d, d + 1);
end

%% Plotting the results

% Setting up the plots folder
curr_path = pwd;
Expmnt = 'Experiment_DFO';
% PltDir = sprintf('%s/Plots/%s/%s', curr_path,Expmnt,data);
PltDir = sprintf('%s/Plots/%s', curr_path, data);
if ~exist(PltDir, 'dir')
    mkdir(PltDir);
end
% min_val=Opt_Val;
threshold = 10^-16;
min_val = min([min(funvalsepoch3(:)) min(funvalsepoch4(:)) Opt_Val]) - threshold; % threshold to ignore negative data

%% New changes after revision
% Reporting the best (minimum) average F(x) - F* for each method
BestFDNorm = min(funvalsepoch3(:, 3)) - min_val;
BestFDIPQN = min(funvalsepoch4(:, 3)) - min_val;
fprintf('Best average F(x) - F*  for FD-Norm method is %.4e \n', BestFDNorm);
fprintf('Best average F(x) - F*  for FD-IPQN method is %.4e \n', BestFDIPQN);
%% Figure 1: Function Values vs Effective Gradient Evaluations

options.ylabel = '$F(x) - F^*$';
options.xlabel = 'Total Function Evaluations';
options.legend = {"FD-Norm", "FD-IPQN"};
options.lineStyles = {'-', '-'};
options.markers = {'s', '*'};

options.colors = [0 1 0
    0.9290  0.6940 0.1250
    1 0 1
    ];
options.legendLoc = 'Best'; % 'NorthEast';

str11 = sprintf('%s', data);
str14 = ',   \theta=';
str15 = sprintf('%.2f', theta);
str1 = strcat(str11, str14, str15);
% options.title = str1;
% Uncomment if you want a title
% options.title = str11;

options.logScale = 2;
options.lineWidth = 2;
% options.xlimits=[0 100];
% options.ylimits=[10^-8 10];

% figure;
% options.xlimits=[0 100];
options.ylimits = [max([min([min(funvalsepoch3(:) - min_val) min(funvalsepoch4(:) - min_val)])]) Inf];
% Overwrite just for short version
% options.ylimits=[10^-12.5 10^-1];
% options.ytick = [10.^[-11:3:-2]];

% xdata={Evals1*(d+1),Evals31*(d+1) + FEvals31,Evals3*(d+1) + FEvals3,Evals41*(d+1) + FEvals41,Evals4*(d+1) + FEvals4};
% xstart=1;
% epoch1(1)=xstart;epoch2(1)=xstart;epoch3(1)=xstart;epoch4(1)=xstart;
% xdata={epoch1,epoch2,epoch3,epoch4};
ydata = {funvalsepoch3(:, 3) - min_val, funvalsepoch4(:, 3) - min_val};
errors{1, 1} = funvalsepoch3(:, 1) - min_val;
errors{1, 2} = funvalsepoch3(:, 2) - min_val;
errors{2, 1} = funvalsepoch4(:, 1) - min_val;
errors{2, 2} = funvalsepoch4(:, 2) - min_val;
options.errors = errors;
options.errorFill = 1;
options.markers = [];
options.markerSpacing = [];
options.xlimits = [0 xlimit_FEvals];
options.markers = {'s', '*'};
options.markerSize = 10;
options.markerIndices = {ExpBudgetMarkers3(1:5:end), ExpBudgetMarkers4(1:5:end)};
% prettyPlot(xdata, ydata, options)
% set(gca,'YTick',[10.^[-11:3:-2]],'FontSize', 22)

% str5=sprintf('%s/%s_Experiment1_DFO_Funcvals_%i_%i_%i_%i_%s.pdf',PltDir,data,log10(lambda),S0_SG,S0_Ada,gamma*100,method);
% print(str5, '-dpdf');

% log-log plot
% in log-log xaxis cannot start with 0
xstart = 1;
epoch30 = epoch3; epoch40 = epoch4;
epoch30(1) = xstart; epoch40(1) = xstart;
xdata = {epoch30, epoch40};
figure;
options.logScale = 3;
options.markers = {'s', '*'};
options.markerSize = 10;
options.markerIndices = {ExpBudgetMarkers3(1:5:end), ExpBudgetMarkers4(1:5:end)};
prettyPlot(xdata, ydata, options);
str5 = sprintf('%s/%s_Experiment1_DFO_Funcvals_%i_%i_%i_log_%i_%s_Ada.pdf', PltDir, data, log10(lambda), S0_SG, S0_Ada, gamma * 100, method);
print(str5, '-dpdf');

% Wide plot (not plotting in the production code)
% figure;
options.logScale = 2;
% xdata={epoch1,epoch2,epoch3,epoch4};
options.xlimits = [0 xlimitwide];
% prettyPlot(xdata, ydata, options)
% set(gca,'YTick',[10.^[-11:3:1]],'FontSize', 22)
% str5=sprintf('%s/%s_Experiment1_DFO_Funcvals_%i_%i_%i_wide_%i_%s.pdf',PltDir,data,log10(lambda),S0_SG,S0_Ada,gamma*100,method);
% print(str5, '-dpdf');

% % Function values vs iterations
% figure;
% options.xlabel = 'Iterations';
% options.xlimits=[1 10000];
% xdata={iter1,iter31,iter3,iter41,iter4};
% ydata={funvals1-min_val,funvals31-min_val,funvals3-min_val,funvals41-min_val,funvals4-min_val};
% prettyPlot(xdata, ydata, options)
% set(gca,'YTick',[10.^[-11:3:-2]],'FontSize', 22)
% str5=sprintf('%s/%s_Experiment2_AdaLBFGS_CuterDFO_Funcvals_Iters_%i_%i_%i_%s.pdf',PltDir,data,log10(lambda),S0_SG,S0_Ada,variance);
% print(str5, '-dpdf')

% Batchsizes
options.xlabel = 'Iterations';
options.markers = {'s', '*'};
% options.markerSpacing = [25 1
%     25 1
%     25 1
%   25 1];
% options.colors = [0 1 0
%     0 1 1
%     ];
options.colors = [0 1 0
    0.9290  0.6940 0.1250
    ];
options.xlabel = 'Iterations (k)';
% options.legend = {"FD-Norm" + method,"FD-IPQN" + method};
options.legend = {"FD-Norm", "FD-IPQN"};

options.title = [];
options.markers = {'s', '*'};
options.markerSpacing = [];
options.markerIndices = {[], []};
figure;
options.ylabel = 'Batch Size $\left(|S_k|\right)$';
options.xlimits = [];
options.ylimits = [];
xdata = {iter3, iter4};
ydata = {batch3(:, 3), batch4(:, 3)};
errors = {};
errors{1, 1} = batch3(:, 1);
errors{1, 2} = batch3(:, 2);
errors{2, 1} = batch4(:, 1);
errors{2, 2} = batch4(:, 2);
options.errors = errors;
options.markers = [];
options.markerSpacing = [];
options.markerSize = [];
options.markerIndices = [];
prettyPlot(xdata, ydata, options);
str5 = sprintf('%s/%s_Experiment1-DFO_Batchsizes_%i_%i_%s_Ada.pdf', PltDir, data, log10(lambda), gamma * 100, method);
print(str5, '-dpdf');

% Stepsizes
figure;
options.ylabel = 'Step Length $(\alpha_k)$';
xdata = {iter3, iter4};
ydata = {steps3(:, 3), steps4(:, 3)};
errors = {};
errors{1, 1} = steps3(:, 1);
errors{1, 2} = steps3(:, 2);
errors{2, 1} = steps4(:, 1);
errors{2, 2} = steps4(:, 2);
options.errors = errors;
options.markers = [];
options.markerSpacing = [];
options.markerSize = [];
options.markerIndices = [];
prettyPlot(xdata, ydata, options);
str5 = sprintf('%s/%s_Experiment1-DFO_Stepsizes_%i_%i_%s_Ada.pdf', PltDir, data, log10(lambda), gamma * 100, method);
print(str5, '-dpdf');
