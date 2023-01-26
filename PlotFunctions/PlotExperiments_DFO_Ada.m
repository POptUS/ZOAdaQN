% This file generates plots comparing the results of Finite-Difference
% Adaptive Sampling Quasi-Newton methods 

%%

% Setting up default values
if ~exist('datas', 'var')
    datas = {'15-absnormal'};
end

if ~exist('sigmas', 'var')
    sigmas = {10^-3};
end

if ~exist('loss', 'var')
    loss = 'CuterDFO';
end
if ~exist('rand_runs_adamethods', 'var')
    rand_runs_adamethods = 5;
end

theta = 0.9;
gamma = 0.9;
S0_SG = 2;
S0_Ada = 2;
xlimit = 5 * 10^5;
xlimitwide = 10000;
method = "-AD1";

% Look up table of optimal stepsizes for each case
variables_data = containers.Map;
data_lib = {'15-absnormal', '15-relnormal', '18-absnormal', '18-relnormal', ...
    '19-absnormal-50', '19-relnormal-50', '20-absnormal', '20-relnormal', ...
    '22-absnormal', '22-relnormal', '205-absnormal', '205-relnormal', ...
    '220-absnormal', '220-relnormal', '216-absnormal', '216-relnormal', ...
    '305-absnormal', '305-relnormal', '124-absnormal', '124-relnormal', ...
    '123-absnormal', '123-relnormal'};
variable_lib = {30, 30, 11, 11, 50, 50, 20, 20, 8, 8, 125, 125, 110, 110, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};


for lib_ind = 1:length(data_lib)
    variables_data(data_lib{lib_ind}) = variable_lib{lib_ind};
end
% including one additional key for nonsmooth loss rand-50-50
variables_data('Rand-50-50') = 50;

% Plotting the results for given experimental setting
for dat = 1:length(datas)
    data = datas{dat};
    variables = variables_data(data);
    for ik = 1:length(sigmas)
        lambda = sigmas{ik};
        d = variables;
        SetupExperiment_DFO_Ada(data, lambda, d, theta, S0_SG, S0_Ada, xlimit, xlimitwide, gamma, method, rand_runs_adamethods);
    end
end

