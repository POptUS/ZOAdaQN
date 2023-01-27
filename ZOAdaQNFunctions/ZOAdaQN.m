function [w, output] = ZOAdaQN(Obj, inputvals, Options, w)
    % For license, source, & updates to ZOAdaQN, see https://github.com/POptUS/ZOAdaQN
    % This is the main function file for running the derivative free quasi-
    % Newton optimization algorithms given in the paper by Bollapragada & Wild.
    % The code is written for finite-sum structure; so the original stochastic
    % problem is solved by sampling large (around 1000000) samples from the
    % stochastic distribution.
    %
    % It requires the following inputs as explained below
    %
    % Obj: This is declared as a function which gives the functional value and
    %      the gradient (estimated via ) value for full and batch cases.
    % Obj.func: Gives function and gradient values of the whole function
    % Obj.funcBatch: Gives the sample average objective and sample average
    %                gradient.
    % Obj.test: Gives the test function and correct classification
    % inputvals: This is the input that need to be given into the Obj for
    %            getting the functional values.
    % w:         if possible w (initial weights) should be given
    %            else, initial weight will be taken as 0.
    %
    % Options: This should contain the parameters required in the algorithm
    %          The following are the parameters which need to be specified and
    %          their default values, if not specified.
    %
    %

    % Extracting problem parameters
    ndatapnts = inputvals.ndatapnts; % number of data points
    nvars = inputvals.nvars; % number of variables
    alpha = Options.alpha0; % initial stepsize
    if isempty(w)
        w = zeros(nvars, 1);
    end

    fprintf('Running Optimization Algorithm..\n');
    fprintf('%15s %15s %15s %15s %15s %15s %15s %15s %15s\n', 'Iterations', ...
            'Eff_passes', 'batch_size', 'Sample_Obj', 'Full_Obj', ...
            '||Sample_grad||_2', '||grad||_2', '||grad||_inf', 'alpha');

    % initialization
    NfunEvals = 0;
    NgradEvals = 0;
    S = Options.S0;

    % Random sampling from the finite-sum structure:
    switch Options.RandChoice
        case 'WR' % with replacement
            shuffle = randi(ndatapnts, ndatapnts, 1);
        case {'WOR'} % without replacement
            shuffle = randperm(ndatapnts);
    end

    % for initial iteration
    batch = shuffle(1:S);
    [f_s, g_s] = Obj.funcBatch(w, inputvals, batch);
    [f, g] = Obj.func(w, inputvals);
    % d=-g_s;  % search direction
    NgradEvals = NgradEvals + S; % counting the total gradient evaluations
    [ftest, correct] = Obj.test(inputvals, w); %
    fprintf('%15i %12.3f %15i %15.5e %15.5e %15.5e %15.5e %15.5e %15s\n', ...
    0, 0, S, f_s, f, norm(g_s, 2), norm(g, 2), norm(g, Inf), 'NA');

    % For storage
    if isempty(Options.MaxStorage)
        max_storage = 100;
    else
        max_storage = Options.MaxStorage;
    end

    % Performance Storage
    funvals = zeros(max_storage, 1);
    norm2grad = zeros(max_storage, 1); % l_2 norm of the finite-difference gradient
    norminfgrad = zeros(max_storage, 1); % l_infty norm
    Evals = zeros(max_storage, 1);
    FEvals = zeros(max_storage, 1);
    iterations = zeros(max_storage, 1);
    TestFunc = zeros(max_storage, 1);
    CorrClass = zeros(max_storage, 1);
    time = zeros(max_storage, 1);

    % Sample Storage
    batchsizes = zeros(max_storage, 1);
    Sample_Obj = zeros(max_storage, 1);
    Sample_Grad = zeros(max_storage, 1);
    Sample_Var = zeros(max_storage, 6);
    True_Var = zeros(max_storage, 6);
    varratios = zeros(max_storage, 8);
    Descent_Condition = zeros(max_storage, 1);
    IP_grad_samgrad = zeros(max_storage, 1);
    Norm_Direction = zeros(max_storage, 1);
    delta_vals = zeros(max_storage, 1);

    % Line Search Storage
    stepsizes = zeros(max_storage, 1);
    LSEvals = zeros(max_storage, 1);
    zeta_vals = zeros(max_storage, 1);
    LS_Skip = zeros(max_storage, 1);

    % Overlapping Storage
    ovpsizes = zeros(max_storage, 1);
    biassizes = zeros(max_storage, 1);
    wastesizes = zeros(max_storage, 1);

    % Curvature Storage
    Curvature = zeros(max_storage, 1);
    Curv_Skip = zeros(max_storage, 1);
    Ip_S_BS = zeros(max_storage, 1);

    % Initialization
    alphak = 1;
    zeta1 = 1;
    iter = 0;
    time(1, 1) = 0;
    LineSearch = Options.LineSearch;
    % alphak_old=1;
    S_old = S;

    % Curvature
    history = [];
    history.n = size(g_s, 1);
    history.y = [];
    history.s = [];
    history.rho = [];
    m = Options.m;
    g_new = g_s;
    S_ovp = max(1, floor(Options.ovp * S / 100));

    % Storage:
    c_stor = 0;
    stor_inter = Options.StoreInterval; % store every 1/10^12th of an epoch
    count_stor_eval = NgradEvals;

    % Adaptive Theta
    theta0 = Options.theta;
    theta = theta0;

    gamma = Options.gamma;

    while NgradEvals < Options.MaxEpochs * ndatapnts
        tic;
        iter = iter + 1;
        if iter == 1 || S > S_old || S < S_old || ...
         (NgradEvals - count_stor_eval > stor_inter * ndatapnts || S == ndatapnts)
            str_flag = 1;                 % printing the values
            count_stor_eval = NgradEvals;
        else
            str_flag = 0;
        end

        if str_flag % performance statistics
            c_stor = c_stor + 1;
            iterations(c_stor, 1) = iter;
            batchsizes(c_stor, 1) = S;
            funvals(c_stor, 1) = f;
            norminfgrad(c_stor, 1) = norm(g, Inf);
            norm2grad(c_stor, 1) = norm(g);
            Evals(c_stor, 1) = NgradEvals;
            FEvals(c_stor, 1) = NfunEvals;
            Sample_Obj(c_stor, 1) = f_s;
            Sample_Grad(c_stor, 1) = norm(g_s, 2);
            TestFunc(c_stor) = ftest;
            CorrClass(c_stor) = correct;
            % Inner product of g_s*g:
            IP_grad_samgrad(c_stor, 1) = g_s' * g / (Sample_Grad(c_stor, 1) * norm2grad(c_stor, 1));
            if strcmp(Options.StopTest, 'InftyNorm') == 1
                termination = norminfgrad(c_stor, 1);
            elseif strcmp(Options.StopTest, 'l2Norm') == 1
                termination = norm2grad(c_stor, 1);
            elseif strcmp(Options.StopTest, 'l2_d') == 1
                termination = norm2grad(c_stor, 1) / sqrt(length(w));
            elseif strcmp(Options.StopTest, 'AbsFunc') == 1 && c_stor > 1
                termination = abs(funvals(c_stor, 1) - funvals(c_stor - 1, 1));
            elseif strcmp(Options.StopTest, 'SamAbsFun') == 1 && c_stor > 1
                termination = abs(Sample_Obj(c_stor, 1) - Sample_Obj(c_stor - 1, 1));
            elseif strcmp(Options.StopTest, 'LSFail') == 1 && c_stor > 2
                termination = min(stepsizes(1:c_stor - 1, 1));
            else
                termination = Options.epsilon;
            end
            if termination < Options.epsilon
                fprintf('Termination crietria %s is less than tolerance %.2e\n', Options.StopTest, Options.epsilon);
                % fprintf('Optimality conditions are achieved\n')
                break
            end
        end

        if strcmp(Options.Curvature, 'Yes') == 1 % quasi-Newton directions
            k = size(history.y, 1);
            H0 = GetH0(history, k);
            d = -Twolooprecurrsion(H0, history, g_s, k);
        else
            d = -g_s;
        end

        % Code for Overlap
        if S == ndatapnts
            S_ovp = S;
        else
            switch Options.ovp_type
                case 'MB' % multi-batch approach
                    switch Options.ovp_rule
                        case 'Fd' % Fixed
                            S_ovp = max(1, floor(Options.ovp * S / 100));
                    end
                case 'FO' % Full overlap as given in the paper
                    S_ovp = S;
            end
            % Fool Proof - Ensure correct overlapping samples
            if Options.ovp == 100
                S_ovp = S;
            end
        end

        % Evaluate the overlapping functions and gradient values
        if S_ovp == S
            batchovp = batch;
            [~, g_ovp] = Obj.funcBatch(w, inputvals, batchovp);
            % NgradEvals=NgradEvals + length(batchovp);
        else
            batchovp = batch(randperm(S, S_ovp));
            [~, g_ovp] = Obj.funcBatch(w, inputvals, batchovp);
            % NgradEvals=NgradEvals + length(batchovp);
        end

        % For L-BFGS - store previous iterates
        w_old = w;
        g_old = g_ovp;
        g_old_step = g_s;

        if str_flag % performance statistics
            Descent_Condition(c_stor, 1) = -d' * g;
            Norm_Direction(c_stor, 1) = norm(d);
            ovpsizes(c_stor, 1) = S_ovp;
        end

        % Line Search
        if strcmp(Options.Initial_SteplengthRule, 'Con') == 1
            zeta1 = 1;
        end

        if S == ndatapnts % For full objective reinitialize alpha=1 with no contracting factors
            zeta1 = 1;
            alpha = 1;
        end

        if str_flag % performance statistics
            zeta_vals(c_stor, 1) = zeta1;
        end

        switch LineSearch
            case 'Constant'
                alphak = alpha * zeta1; % Decrease the initial estimate by zeta1
                w = w + alphak * d;
                [f_new, g_new] = Obj.funcBatch(w, inputvals, batchovp);
                NLineEvals = 0;
            case 'Armijo'
                alphak = alpha * zeta1; % Decrease the initial estimate by zeta1
                [f_new, g_new, alphak, NLineEvals] = Armijo(Obj, inputvals, batch, f_s, d, ...
                g_s, w, Options.c, Options.rho, alphak, Options.tolr, Options.LineSearch_c2);
                w = w + alphak * d;
        end
        if S == ndatapnts
            g_ovp_new = g_new;
        else
            switch Options.ovp_type
                case 'MB'
                    [~, g_ovp_new] = Obj.funcBatch(w, inputvals, batchovp);
                case 'FO'
                    g_ovp_new = g_new;
            end
        end
        % alphak_old=alphak;

        if str_flag % Performance Statistics
            LSEvals(c_stor, 1) = NLineEvals;
            stepsizes(c_stor, 1) = alphak;
        end
        NfunEvals = NfunEvals + length(batch) * NLineEvals; % counting function evaluations in linesearch
        if str_flag % Performance Statistics
            [ftest, correct] = Obj.test(inputvals, w);
            fprintf('%15i %12.6f %15i %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n', ...
            iter, NgradEvals / ndatapnts, S, f_s, f, Sample_Grad(c_stor, 1), ...
            norm2grad(c_stor, 1), norminfgrad(c_stor, 1), alphak);
            time(c_stor + 1, 1) = time(c_stor, 1) + toc;
            [f, g] = Obj.func(w, inputvals);
        end

        % Selecting the new gradient
        if S == ndatapnts
            f_s = f_new;
            g_s = g_new;
            batchsample = [1:ndatapnts]';
            NgradEvals = NgradEvals + length(batchsample);
        else
            switch Options.ovp_type
                case 'MB' % multi-batch approach instead of full-overlap
                    batchsample = randperm(ndatapnts, S);
                    [f_s, g_s] = Obj.funcBatch(w, inputvals, batchsample);
                    NgradEvals = NgradEvals + length(batchsample);
                case 'FO'
                    switch Options.ovp_rule
                        case 'Fd' % Fixed
                            % S_bias=max(1, floor(Options.ovp*S/100));
                            S_bias = floor(Options.ovp * S / 100);
                    end
                    % correcting for even one sample overlap
                    if S_bias > 0
                        batchbias = batch(randperm(S, S_bias));
                        batchrem0 = setdiff([1:ndatapnts]', batchbias);
                        batchsample = batchrem0(randperm(length(batchrem0), S - S_bias));
                        batchsample = union(batchsample, batchbias);
                        [f_s, g_s] = Obj.funcBatch(w, inputvals, batchsample);
                        NgradEvals = NgradEvals + length(setdiff(batch, batchsample)) + length(batch);
                    else
                        S_bias = 0;
                        batchsample = randperm(ndatapnts, S);
                        [f_s, g_s] = Obj.funcBatch(w, inputvals, batchsample);
                        NgradEvals = NgradEvals + length(setdiff(batch, batchsample)) + length(batch);
                        % NgradEvals=NgradEvals + length(setdiff(batch,
                        % batchsample)) + length(batch) + length(batchsample);
                        % This is the correct total work for finite-difference
                        % based gradient estimation. We can post process with
                        % the batchsizes computation
                    end
                    if str_flag % Performance Statistics
                        biassizes(c_stor, 1) = S_bias;
                        wastesizes(c_stor, 1) = length(setdiff(batch, batchsample));
                    end
            end
        end
        batch = batchsample;

        % Verifying if we need to increase sample size.
        if (strcmp(Options.SamplingTest, 'Yes') == 1) && S < ndatapnts && S > 1

            if strcmp(Options.Variance, 'Exact') == 1 % Use all the samples in estimating the variance
                batchsample = [1:ndatapnts]';
                [f, g] = Obj.func(w, inputvals);
                mean = g;
            elseif strcmp(Options.Variance, 'Sample') == 1
                mean = g_s;
            else
                batchsample = batch;
            end

            batchsample = batchsample(randperm(length(batchsample), min(length(batchsample), Options.SamVarSize)));

            switch Options.Method
                case {'Norm', 'Norm-AD1', 'Norm-AD2', 'Norm-AD3', 'Norm-AD4', 'Norm-AD5', 'Norm-AD6', 'Norm-AD7'}
                    [samvar, ~] = Normvariance(batchsample, mean, Obj, inputvals, w);
                    samvar = norm(samvar, 1);
                    varbiasratio = norm(samvar, 1) / (S * (theta * norm(mean, 2))^2);
                    varbiasratiot = varbiasratio;
                    if strcmp(Options.Initial_SteplengthRule, 'Norm') == 1
                        zeta = varbiasratio * theta^2 / max(varbiasratio, 1);
                        zeta1 = min(1, 1 / (1 + zeta));
                    elseif strcmp(Options.Initial_SteplengthRule, 'Norm-UB') == 1 % DSS unbiased estimate
                        % zeta=varbiasratio*Options.theta^2/max(varbiasratio,1);
                        % R/(1 - R) for unbiasedness:
                        zeta = (varbiasratio * theta^2) / (1 - (varbiasratio * Options.theta^2));
                        zeta1 = min(1, 1 / (1 + zeta));
                    elseif strcmp(Options.Initial_SteplengthRule, 'Con') == 1
                        zeta1 = 1;
                    end

                    if str_flag
                        varratios(c_stor, 1) = varbiasratio;
                        varratios(c_stor, 7) = varbiasratio;
                        Sample_Var(c_stor, 1) = norm(samvar, 1);
                    end
                case {'IPQN', 'IPQN-AD1', 'IPQN-AD2', 'IPQN-AD3', 'IPQN-AD4', 'IPQN-AD5', 'IPQN-AD6', 'IPQN-AD7'}
                    bias = Twolooprecurrsion(H0, history, mean, k);
                    Hip = Twolooprecurrsion(H0, history, bias, k);
                    [samvaripQN] = IPQNvariance(batchsample, bias, Obj, inputvals, w, H0, history, mean, k);
                    varbiasratio = samvaripQN / (S * (theta^2) * (Hip' * mean)^2);
                    varbiasratiot = varbiasratio;
                    if strcmp(Options.Initial_SteplengthRule, 'Norm') == 1
                        [samvarDSS, ~] = Normvariance(batchsample, mean, Obj, inputvals, w); % working on it
                        samvarDSS = norm(samvarDSS, 1);
                        varbiasratioDSS = norm(samvarDSS, 1) / (S * (theta * norm(mean, 2))^2);
                        zeta = varbiasratioDSS * theta^2 / max(varbiasratio, 1);
                        zeta1 = min(1, 1 / (1 + zeta));
                        if str_flag
                            varratios(c_stor, 1) = varbiasratioDSS;
                            Sample_Var(c_stor, 1) = norm(samvarDSS, 1);
                        end
                    else
                        zeta1 = 1;
                    end

                    if str_flag
                        varratios(c_stor, 4) = varbiasratio;
                        varratios(c_stor, 7) = varbiasratio;
                        Sample_Var(c_stor, 4) = samvaripQN;
                    end
            end

            % Checking if we need to decrease the sample size
            if strcmp(Options.Sample_decrease, 'Yes') == 1
                varbiasratio = max(varbiasratio, 0.5);
                if ceil(S * varbiasratio) < 2 % safeguarding the sample size to be greater than 2
                    varbiasratio = 2 / S;
                end
            else
                varbiasratio = max(varbiasratio, 1);
            end
            varbiasratio = min(varbiasratio, Options.Maxvarbias);

        else
            varbiasratio = 1;

            % Code for constant batch size with adaptive steplength
            if S < ndatapnts && S > 1
                if strcmp(Options.Variance, 'Exact') == 1
                    batchsample = [1:ndatapnts]';
                    [f, g] = Obj.func(w, inputvals);
                    % [f,g]=Obj.func(w,inputvals);
                    mean = g;
                elseif strcmp(Options.Variance, 'Sample') == 1
                    mean = g_s;
                else
                    batchsample = batch;
                end
                batchsample = batchsample(randperm(length(batchsample), min(length(batchsample), Options.SamVarSize)));
                if strcmp(Options.Initial_SteplengthRule, 'Norm') == 1
                    samvarDSS = samplevariance(batchsample, mean, Obj, inputvals, w);
                    varbiasratioDSS = norm(samvarDSS, 1) / (S * (Options.theta * norm(mean, 2))^2);
                    zeta = varbiasratioDSS * Options.theta^2 / max(varbiasratio, 1);
                    zeta1 = min(1, 1 / (1 + zeta));
                    if str_flag
                        varratios(c_stor, 1) = varbiasratioDSS;
                        Sample_Var(c_stor, 1) = norm(samvarDSS, 1);
                    end
                else
                    zeta1 = 1;
                end
            end
        end

        if str_flag % performance statistics
            varratios(c_stor, 8) = varbiasratio;
        end

        S_old = S;
        S = min(ceil(varbiasratio * S), ndatapnts);

        % Increasing Sample Size
        if S > S_old
            batchrem = setdiff([1:ndatapnts]', batch);
            batch_new = batchrem(randperm(length(batchrem), S - S_old));
            [f_temp, g_temp] = Obj.funcBatch(w, inputvals, batch_new);
            NgradEvals = NgradEvals + length(batch_new);
            f_s = (f_s * S_old + f_temp * (S - S_old)) / S;
            g_s = (g_s * S_old + g_temp * (S - S_old)) / S;
            batch = union(batch, batch_new);
        elseif S < S_old % Decreasing sample size
            batch = batch(randperm(S_old, S));
            [f_s, g_s] = Obj.funcBatch(w, inputvals, batch);

        end

        % Adaptive Theta formula
        switch Options.Method
            case {'IPQN-AD1', 'Norm-AD1'}
                if iter > 1 && S <= S_old
                    theta = theta * gamma;
                else
                    theta = theta0;
                end
            case {'IPQN-AD2', 'Norm-AD2'}
                if iter > 1 && S <= S_old
                    theta = theta * gamma;
                else
                    theta = min(theta0, theta / gamma);
                end
            case {'IPQN-AD3', 'Norm-AD3'}
                if iter > 1 && S <= S_old
                    theta = theta * gamma;
                else
                    theta = theta / gamma;
                end
            case {'IPQN-AD4', 'Norm-AD4'}
                if iter > 1 && S <= S_old
                    theta = theta * sqrt(varbiasratiot);
                else
                    theta = theta * sqrt(S / S_old);
                end
            case {'IPQN-AD5', 'Norm-AD5'}
                if iter > 1 && S <= S_old
                    theta = theta * gamma;
                else
                    theta = theta * sqrt(S / S_old);
                end
            case {'IPQN-AD6', 'Norm-AD6'}
                if iter > 1 && S <= S_old
                    theta = theta * gamma;
                else
                    theta = min(theta0, theta * sqrt(S / S_old));
                end
            case {'IPQN-AD7', 'Norm-AD7'}
                if iter > 1 && S <= S_old
                    theta = theta * max(gamma, sqrt(varbiasratiot));
                else
                    theta = min(theta0, theta * sqrt(S / S_old));
                end
        end
        % Curvature Update
        if strcmp(Options.Curvature, 'Yes') == 1 && (S >= Options.threshold * ndatapnts)

            s_new = w - w_old;
            y_new = g_ovp_new - g_old;
            if str_flag
                Curvature(c_stor, 1) = y_new' * s_new;
                Ip_S_BS(c_stor, 1) = alphak^2 * g_old_step' * d;
            end
            switch Options.curvskiprule
                case 'wlf'
                    skp_rule = (y_new' * s_new) / (Options.skip * norm(s_new)^2);
                case 'sp'
                    skp_rule = (y_new' * s_new) / (-Options.skip * alphak^2 * g_old_step' * d);
                otherwise
                    skp_rule = 1;
            end
            % if y_new'*s_new >= Options.skip
            % if y_new'*s_new >= -Options.skip*alphak^2*g_old_step'*d % Skipping update
            % if y_new'*s_new >= Options.skip*norm(y_new)*norm(s_new) % Skipping update
            % if y_new'*s_new > Options.skip*norm(s_new)^2 && norm(y_new) < (Options.MaxL)*norm(s_new)% Skipping update
            if skp_rule >= 1 && norm(y_new) < (Options.MaxL) * norm(s_new) % Skipping update
                if size(history.y, 1) > Options.m
                    y = history.y;
                    s = history.s;
                    rho = history.rho;
                    y(1, :) = [];
                    s(1, :) = [];
                    rho(1) = [];
                    y(m, :) = y_new';
                    s(m, :) = s_new';
                    rho(m) = 1 / (y(m, :) * s(m, :)');
                    history.y = y;
                    history.s = s;
                    history.rho = rho;
                elseif size(history.y, 1) <= Options.m  % && (Curvature(end,1) > 0) && (S>=Options.threshold*ndatapnts)
                    history.y(end + 1, :) = y_new';
                    history.s(end + 1, :) = s_new';
                    history.rho(end + 1) = 1 / (history.y(end, :) * history.s(end, :)');
                end
            else
                if str_flag
                    Curv_Skip(c_stor, 1) = 1;   % If we aren't storing then we won't get exact stat.
                end
                % disp('skip')
            end
        end
    end

    %% Output Storage

    % Performance Statistics
    output.sol = w;
    output.iterations = iterations(1:c_stor);
    output.funvals = funvals(1:c_stor);
    output.norm2grad = norm2grad(1:c_stor);
    output.norminfgrad = norminfgrad(1:c_stor);
    output.Evals = Evals(1:c_stor);
    output.FEvals = FEvals(1:c_stor);
    output.time = time(1:c_stor);
    output.testfunvals = TestFunc(1:c_stor);
    output.corrclass = CorrClass(1:c_stor);

    % Line Search Statistics
    output.alpha = stepsizes(1:c_stor);
    output.lsEvals = LSEvals(1:c_stor);
    output.lsFails = LS_Skip(1:c_stor); % iterations where line search failed
    output.zetavals = zeta_vals(1:c_stor, 1);

    % Sample Statistics
    output.batchsizes = batchsizes(1:c_stor);
    output.sampleobj = Sample_Obj(1:c_stor);
    output.samplegrad = Sample_Grad(1:c_stor);
    output.varratios = varratios(1:c_stor, 1:8);
    output.samplevar = Sample_Var(1:c_stor, 1:6);
    output.truevar = True_Var(1:c_stor, 1:6);
    output.descent = Descent_Condition(1:c_stor);
    output.normdirection = Norm_Direction(1:c_stor);
    output.gradtimessamgrad = IP_grad_samgrad(1:c_stor);
    output.deltavals = delta_vals(1:c_stor, 1);

    % Curvature Statistics
    output.curvskips = Curv_Skip(1:c_stor); % iterations where curvature updates are skipped
    output.curvature = Curvature(1:c_stor); % y^Ts;
    output.curvaturecomp = Ip_S_BS(1:c_stor); % s^TBs;

    % Overlap statistics
    output.s_ovp = ovpsizes(1:c_stor);
    output.s_bias = biassizes(1:c_stor);
    output.s_waste = wastesizes(1:c_stor);

end
