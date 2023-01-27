function [w, output] = DeterministicQuasiNewtonAlgorithms(Obj, inputvals, Options, w)
    % This is the main file for running the Deterministic Optimization algorithms.
    % It requires the following inputs as explained below.
    %
    % Obj: This is declared as a function which gives the functional value and
    %      the gradient value too for full and batch codes. That is
    % Obj.func: Gives function and gradient values of the whole function
    %
    %
    % inputvals: This is the input that need to be given into the Obj for
    %            getting the functional values. Usually, for logistic
    %            regression, it has X,y. This should at least contain:
    %               ndatapnts (N): Total number of component functions.
    %               nvars(n): Total number of variables
    % w:         if possible w (initial weights) should be given
    %            else, initial weight will be taken as 0.
    %
    % Options: This should contain the parameters required in the algorithm
    %          The following are the parameters which need to be specified and
    %          their default values, if not specified.
    %
    %

    nvars = inputvals.nvars;
    alpha = Options.alpha0;
    if isempty(w)
        w = zeros(nvars, 1);
    end

    fprintf('Running Optimization Algorithm..\n');
    fprintf('%15s %15s %15s %15s %15s %15s\n', ...
            'Iterations', 'Eff_passes', 'Full_Obj', '||grad||_2', '||grad||_inf', 'alpha');

    NfunEvals = 0;
    NgradEvals = 0;

    % for initial iteration
    [f, g] = Obj.func(w, inputvals);
    d = -g;
    NgradEvals = NgradEvals + 1;
    [ftest, correct] = Obj.test(inputvals, w);
    fprintf('%15i %12.3f %15.5e %15.5e %15.5e %15s\n', 0, 0, f, norm(g, 2), norm(g, Inf), 'NA');

    % For storage
    if isempty(Options.MaxStorage)
        max_storage = 100;
    else
        max_storage = Options.MaxStorage;
    end

    % Performance Storage
    funvals = zeros(max_storage, 1);
    norm2grad = zeros(max_storage, 1);
    norminfgrad = zeros(max_storage, 1);
    Evals = zeros(max_storage, 1);
    FEvals = zeros(max_storage, 1);
    iterations = zeros(max_storage, 1);
    TestFunc = zeros(max_storage, 1);
    CorrClass = zeros(max_storage, 1);
    time = zeros(max_storage, 1);
    Norm_Direction = zeros(max_storage, 1);

    % Line Search Storage
    stepsizes = zeros(max_storage, 1);
    LSEvals = zeros(max_storage, 1);
    zeta_vals = zeros(max_storage, 1);
    LS_Skip = zeros(max_storage, 1);

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

    % Curvature
    history = [];
    history.n = size(g, 1);
    history.y = [];
    history.s = [];
    history.rho = [];
    m = Options.m;

    % Storage: Store every 1/20th of an epoch and
    c_stor = 0; % for storing values
    stor_inter = Options.StoreInterval; % store every 1/10^12th of an epoch
    count_stor_eval = NgradEvals;

    while NgradEvals < Options.MaxEpochs
        tic;
        iter = iter + 1;
        if iter == 1 || (NgradEvals - count_stor_eval > stor_inter)
            str_flag = 1;                 % printing the values
            count_stor_eval = NgradEvals;
        else
            str_flag = 0;
        end

        if str_flag
            c_stor = c_stor + 1;
            iterations(c_stor, 1) = iter;
            funvals(c_stor, 1) = f;
            norminfgrad(c_stor, 1) = norm(g, Inf);
            norm2grad(c_stor, 1) = norm(g);
            Evals(c_stor, 1) = NgradEvals;
            FEvals(c_stor, 1) = NfunEvals;
            TestFunc(c_stor) = ftest;
            CorrClass(c_stor) = correct;
            if strcmp(Options.StopTest, 'InftyNorm') == 1
                termination = norminfgrad(c_stor, 1);
            elseif strcmp(Options.StopTest, 'l2Norm') == 1
                termination = norm2grad(c_stor, 1);
            elseif strcmp(Options.StopTest, 'l2_d') == 1
                termination = norm2grad(c_stor, 1) / sqrt(length(w));
            else
                termination = Options.epsilon;
            end
            if termination < Options.epsilon
                fprintf('norm of gradient is less than tolerance\n');
                fprintf('Optimality conditions are achieved\n');
                break
            end
        end

        if strcmp(Options.Curvature, 'Yes') == 1
            k = size(history.y, 1);
            H0 = GetH0(history, k);
            d = -Twolooprecurrsion(H0, history, g, k);
        else
            d = -g;
        end
        % For L-BFGS
        w_old = w;
        g_old = g;
        g_old_step = g;

        if str_flag
            Norm_Direction(c_stor, 1) = norm(d);
        end

        % Line Search
        if strcmp(Options.Initial_SteplengthRule, 'Con') == 1
            zeta1 = 1;
        end

        % Reducing the steplength for first iteration
        if iter == 1
            zeta1 = 1 / norm(d);
        else
            zeta1 = 1;
        end
        if str_flag
            zeta_vals(c_stor, 1) = zeta1;
        end

        switch LineSearch
            case 'Constant'
                alphak = alpha * zeta1;
                w = w + alphak * d;
                [f, g] = Obj.func(w, inputvals);
                NLineEvals = 0;
            case 'DeterministicArmijo'
                alphak = alpha * zeta1;
                [f, g, alphak, NLineEvals] = DeterministicArmijo(Obj, inputvals, ...
                f, d, g, w, Options.c, Options.rho, alphak, Options.tolr);
                w = w + alphak * d;
        end
        if str_flag
            LSEvals(c_stor, 1) = NLineEvals;
            stepsizes(c_stor, 1) = alphak;
        end
        NfunEvals = NfunEvals + NLineEvals;
        NgradEvals = NgradEvals + 1;
        if str_flag
            [ftest, correct] = Obj.test(inputvals, w);
            fprintf('%15i %12.6f %15.5e %15.5e %15.5e %15.5e\n', ...
            iter, NgradEvals, f, norm2grad(c_stor, 1), norminfgrad(c_stor, 1), alphak);
            time(c_stor + 1, 1) = time(c_stor, 1) + toc;
        end

        % Curvature Update
        if strcmp(Options.Curvature, 'Yes') == 1

            s_new = w - w_old;
            y_new = g - g_old;
            if str_flag
                Curvature(c_stor, 1) = y_new' * s_new;
                Ip_S_BS(c_stor, 1) = alphak^2 * g_old_step' * d;
            end
            % if y_new'*s_new >= Options.skip
            if y_new' * s_new >= -Options.skip * alphak^2 * g_old_step' * d % Skipping update
                % if y_new'*s_new >= Options.skip*norm(y_new)*norm(s_new) % Skipping update
                % if y_new'*s_new >= Options.skip*norm(s_new)^2 % Skipping update
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
    % Adding this for safeguarding of exceeding dimensions
    % if iter > max_storage
    %    iter=max_storage;
    % end

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
    output.normdirection = Norm_Direction(1:c_stor);

    % Curvature Statistics
    output.curvskips = Curv_Skip(1:c_stor); % iterations where curvature updates are skipped
    output.curvature = Curvature(1:c_stor); % y^Ts;
    output.curvaturecomp = Ip_S_BS(1:c_stor); % s^TBs;

end
