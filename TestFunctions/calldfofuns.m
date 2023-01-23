% 8/27/19
% Sample calling syntax for the dfo funs from
% https://github.com/POptUS/BenDFO
%

clear all;
clf; % careful!

% Define Var array that shows problem specifications (nprob, n, m, xtsrat)
if 1 == 1 % n=30 problems
    Var = [1    30   45    0 % quadratic, m>=n
        2   30   45    0 % quadratic, m>=n
        3   30   45    0 % quadratic, m>=n
        11   30   31    0 % Watson, 31=m>=n
        15   30   45    0 % Chebyquad, m>=n
        16   30   45    0 % Brown, m>=n
        19   30   52    0 % Bdqrtic, n>=5, m = (n-4)*2
        20   30   30    0 % Cube, n>=2; m=n;
        21   30   30    0 % Mancino, n>=2; m=n
        ];
else % Larger problems (no quadratics):
    num = 10;
    Var = [
        15   num * 30   num * 45    0 % Chebyquad, m>=n
        16   num * 30   num * 45    0 % Brown, m>=n
        19   num * 30   (num * 30 - 4) * 2    0 % Bdqrtic, n>=5, m = (n-4)*2
        20   num * 30   num * 30    0 % Cube, n>=2; m=n;
        21   num * 30   num * 30    0 % Mancino, n>=2; m=n
        ];
end

% Need global variables to get correct definitions
global m n nprob probtype

% The following global variables are just to track output f values
global fvals nfev np
% !Should probably initialize fvals
np = 0; % Problem counter
% ! nfev initialized at start of new problem

nrows = size(Var, 1);

protype{1} = 'smooth';
protype{2} = 'reluniform';
protype{3} = 'relnormal';
protype{4} = 'absuniform';
protype{5} = 'absnormal';

S = 1:5;

plotflag = 1; % plot: yes or no?
maxfev = 20; % max f evals to do along direction

for i = 1:nrows

    if plotflag
        % Plotting just to see
        figure(i);
        clf;
        hold on;
        title(num2str(nprob));
    end

    nprob = Var(i, 1);
    n = Var(i, 2);
    m = Var(i, 3);
    factor = 10^(Var(i, 4));
    X0 = dfoxs(n, nprob, factor); % starting point

    % Only for eval example
    rand('state', 1);
    dir = (1e-3) * rand(1, n);

    for ptype = S
        probtype = protype{ptype};

        np = np + 1; % Update problem counter
        nfev = 0; % Reset number of f eval counter

        for j = 1:maxfev
        % for now, overwrite y and let fvals do the work:
            y = calfun_sample(X0 + (j - 1) * dir, 1e-3, 1);
        end

        if plotflag
            % Plotting just to see
            figure(i);
            plot(fvals(:, np));
        end
        % func    [f h] Function handle so that func(x) evaluates f (@calfun)
        % func = @(x) funceval(x,14);
    end
end
