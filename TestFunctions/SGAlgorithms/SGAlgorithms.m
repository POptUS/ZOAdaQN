function [w, output] = SGAlgorithms(Obj, inputvals, Options, w)
% This is the main file for running the comparative optimization algorithms 
% for the ZOAdaQN paper.
% It requires the following inputs as explained below.
%
% Obj: This is declared as a function that gives the functional value and
%      the gradient value too for full and batch caaes. That is
% Obj.func: Gives function and gradient values of the whole function
% Obj.funcBatch: Gives the sample average objective and sample average
%                gradient.
%
% inputvals: This is the input that need to be given into the Obj for
%            getting the functional values. Usually, for logistic
%            regression, it has X,y. This should at least contain:
%               ndatapnts (N): Total number of component functions.
%               nvars(n): Total number of variables.
% w:         if possible w (initial weights) should be given
%            else, initial weight will be taken as 0.
%
% Options: This should contain the parameters required in the algorithm
%          The following are the parameters that need to be specified and
%          their default values, if not specified.
%
%

ndatapnts = inputvals.ndatapnts;
nvars = inputvals.nvars;
alpha = Options.alpha0;
if isempty(w)
    w = zeros(nvars, 1);
end

fprintf('Running Optimization Algorithm..\n');
fprintf('%15s %15s %15s %15s %15s %15s %15s %15s %15s\n', ...
    'Iterations', 'Eff_passes', 'batch_size', 'Sample_Obj', ...
    'Full_Obj', '||Sample_grad||_2', '||grad||_2', '||grad||_inf', 'alpha');

NfunEvals = 0;
S = Options.S0;

switch Options.RandChoice
    case 'WR'
        shuffle = randi(ndatapnts, ndatapnts, 1);
    case {'WOR', 'OP', 'OPSVRG', 'Shu', 'WORC'}
        shuffle = randperm(ndatapnts);
end

% for initial iteration
batch = shuffle(1:S);
shuffle_op = shuffle;
[f_s, g_s] = Obj.funcBatch(w, inputvals, batch);
[f, g] = Obj.func(w, inputvals);
d = -g_s;
NfunEvals = NfunEvals + S;
[ftest, correct] = Obj.test(inputvals, w);
fprintf('%15i %12.3f %15i %15.5e %15.5e %15.5e %15.5e %15.5e %15s\n', ...
         0, 0, S, f_s, f, norm(g_s, 2), norm(g, 2), norm(g, Inf), 'NA');

% For storage
if isempty(Options.MaxStorage)
    max_storage = 100;
else
    max_storage = Options.MaxStorage;
end

batchsizes = zeros(max_storage, 1);
funvals = zeros(max_storage, 1);
norm2grad = zeros(max_storage, 1);
norminfgrad = zeros(max_storage, 1);

Evals = zeros(max_storage, 1);
iterations = zeros(max_storage, 1);
Sample_Obj = zeros(max_storage, 1);
Sample_Grad = zeros(max_storage, 1);
stepsizes = zeros(max_storage, 1);
TestFunc = zeros(max_storage, 1);
CorrClass = zeros(max_storage, 1);
time = zeros(max_storage, 1);
Descent_Condition = zeros(max_storage, 1);
LSEvals = zeros(max_storage, 1);

L = Options.L0;

alphak = 1;
iter = 0;
shuffle = randperm(ndatapnts);

time(1, 1) = 0;
flagtest = 'No';

shuffle = [];

% change Dec 8
varratio = 1;

LineSearch = Options.LineSearch;

% Storage: Store every 1/20th of an epoch and
c_stor = 0; % for storing values
stor_inter = 1 / 5; % store every 1/10^12th of an epoch
stor_inter = Options.StoreInterval;
alphak_old = 1;
S_old = S;
count_ls = 0;
count_stor_eval = NfunEvals;
init_stor = str2num(Options.Init_Stor);

while NfunEvals < Options.MaxEpochs * ndatapnts
    % Printing the values at each iteration
    tic;
    iter = iter + 1;
    if iter == 1 ||  (NfunEvals - count_stor_eval > stor_inter * ndatapnts) || ...
     iter < init_stor % included the storage such for initial iterations
        str_flag = 1;
        count_stor_eval = NfunEvals;
    else
        str_flag = 0;
    end

    if str_flag
        c_stor = c_stor + 1;
        iterations(c_stor, 1) = iter;
        batchsizes(c_stor, 1) = S;
        funvals(c_stor, 1) = f;
        norminfgrad(c_stor, 1) = norm(g, Inf);
        norm2grad(c_stor, 1) = norm(g);
        Evals(c_stor, 1) = NfunEvals;
        Sample_Obj(c_stor, 1) = f_s;
        Sample_Grad(c_stor, 1) = norm(g_s, 2);
        Descent_Condition(c_stor, 1) = -d' * g;
        TestFunc(c_stor) = ftest;
        CorrClass(c_stor) = correct;
        if strcmp(Options.StopTest, 'InftyNorm') == 1
            termination = norminfgrad(c_stor, 1);
        elseif strcmp(Options.StopTest, 'l2Norm') == 1
            termination = norm2grad(c_stor, 1);
        elseif strcmp(Options.StopTest, 'AbsFunc') == 1 && c_stor > 1
            termination = abs(funvals(c_stor, 1) - funvals(c_stor - 1, 1));
        else
            termination = Options.epsilon;
        end
        if termination < Options.epsilon
            % fprintf('norm of gradient is less than tolerance\n')
            fprintf('Optimality conditions are achieved\n');
            break
        end
    end

    d = -g_s;
    w_old = w;
    g_old = g_s;

    switch LineSearch
        case 'Constant'
            alphak = alpha;
            w = w + alphak * d;
            [f_new, g_new] = Obj.funcBatch(w, inputvals, batch);
            NLineEvals = 0;
        case 'Armijo'
            [f_new, g_new, alphak, NLineEvals] = Armijo(Obj, ...
            inputvals, batch, f_s, d, g_s, w, Options.c, Options.rho, Options.alpha0);
            w = w + alphak * d;
    end
    alphak_old = alphak;
    if str_flag
        LSEvals(c_stor, 1) = NLineEvals;
        stepsizes(c_stor, 1) = alphak;
    end
    NfunEvals = NfunEvals + length(batch) * NLineEvals;
    if str_flag
        [ftest, correct] = Obj.test(inputvals, w);
        fprintf('%15i %12.6f %15i %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n', ...
        iter, NfunEvals / ndatapnts, S, f_s, f, Sample_Grad(c_stor, 1), ...
        norm2grad(c_stor, 1), norminfgrad(c_stor, 1), alphak);
        time(c_stor + 1, 1) = time(c_stor, 1) + toc;
        [f, g] = Obj.func(w, inputvals);
    end

    % Computation of the new search direction.
    batchsample = randperm(ndatapnts, S)';
    [f_s, g_s] = Obj.funcBatch(w, inputvals, batchsample);
    NfunEvals = NfunEvals + length(batchsample);
    batch = batchsample;

end
% Adding this for safeguarding of exceeding dimensions
% if iter > max_storage
%    iter=max_storage;
% end
output.sol = w;
output.iterations = iterations(1:c_stor);
output.batchsizes = batchsizes(1:c_stor);
output.funvals = funvals(1:c_stor);
output.norm2grad = norm2grad(1:c_stor);
output.norminfgrad = norminfgrad(1:c_stor);
output.Evals = Evals(1:c_stor);
output.sampleobj = Sample_Obj(1:c_stor);
output.samplegrad = Sample_Grad(1:c_stor);
output.alpha = stepsizes(1:c_stor);
output.descent = Descent_Condition(1:c_stor);
output.time = time(1:c_stor);
output.testfunvals = TestFunc(1:c_stor);
output.corrclass = CorrClass(1:c_stor);

output.lsEvals = LSEvals;
end
