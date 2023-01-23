function [Options] = Initialize()
    % This is an intialization file which contains of all the Options required for
    % running the Adaptive Sampling Project.

    % We are intializing all the parameters here.

    Options.Method = 'DSS';

    Options.lambda = 0;
    Options.InitialWeights = 'Zeros';

    Options.LineSearch = 'Constant';
    Options.Initial_SteplengthRule = [];
    Options.LineSearch_FailSkipRule = [];
    Options.LineSearch_Skipparam = [];
    Options.alpha0 = 1;
    Options.c = [];
    Options.rho = [];
    Options.tolr = [];
    % Options.step_tolr=[]; %not necessary
    Options.step_extra = []; % additional parameter just in case

    Options.RandChoice = 'WOR';
    Options.RandRuns = 1;

    Options.MaxIterations = 150;
    Options.MaxEpochs = 100;
    Options.StoreInterval = 1 / 10^12;
    Options.MaxStorage = 2400000;
    Options.StopTest = [];
    Options.epsilon = [];
    Options.Stop_extra = []; % additional parameter just in case

    Options.S0 = 1;
    Options.SamplingTest = 'No';
    Options.theta = [];
    Options.gamma = [];
    % Options.nu=[];
    Options.Variance = [];
    Options.SamVarSize = [];
    Options.Maxvarbias = [];
    Options.Sample_decrease = []; % whether we allow sample deceases or not...
    Options.Sample_extra = []; % extra parameter just in case

    Options.GeometricRate = [];

    Options.m = [];
    Options.Curvature = 'No';
    Options.threshold = [];
    Options.skip = []; % 0.001 10^-parameter;
    Options.curvature_extra = [];

    Options.ovp_type = [];
    Options.ovp_rule = [];
    Options.ovp = [];

    % DFO initialization
    Options.DFOMethod = [];
    Options.DFOInterval = [];
    Options.DFOIntervalAdaptive = [];
    Options.DFOIntervalfactor = [];
    Options.DFOExtra = [];
end
