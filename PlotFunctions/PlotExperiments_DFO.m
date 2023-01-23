% This file generates plots comparing the results of Finite-Difference
% Adaptive Sampling Quasi-Newton methods with sphere smoothing and
% stochastic finite-difference based gradient estimate methods

% Note that by default the code looks for default optimal stepsizes for
% FD-SG and SS methods in the lookup dictionary (steps_SG_GS) created in 
% this code. 
% If you need to change the stepsizes or use a new dataset please update
% the corresponding value in the ditionary steps_SG_SS
% The keys are 

%%

% Setting up default values
if ~exist('datas','var')
    datas={'15-absnormal'};
end

if ~exist('sigmas','var')
    sigmas={10^-3};
end

if ~exist('loss','var')
    loss='CuterDFO';
end
if ~exist('rand_runs_adamethods','var')
    rand_runs_adamethods=5;
end

theta=0.9;
gamma=0.9;
S0_SG=2;
S0_Ada=2;
xlimit=5*10^5;
xlimitwide=10000;
method="-AD1";

% Look up table of optimal stepsizes for each case
variables_data=containers.Map;
steps_SG_SS=containers.Map;
data_lib={'15-absnormal','15-relnormal','18-absnormal','18-relnormal',...
    '19-absnormal-50','19-relnormal-50','20-absnormal','20-relnormal',...
    '22-absnormal','22-relnormal','205-absnormal','205-relnormal',...
    '220-absnormal','220-relnormal','216-absnormal','216-relnormal',...
    '305-absnormal','305-relnormal','124-absnormal','124-relnormal',...
    '123-absnormal','123-relnormal'};
variable_lib={30,30,11,11,50,50,20,20,8,8,125,125,110,110,100,100,100,100,100,100,100,100,100,100,100,100};
sigma_lib={10^-3,10^-5};

steps_SG_lib=[-10 -10 -5 -5 -13 -13 -12 -10 -11 -10 -11 -8 -9 -9 -8 -8 -10 -10 -6 -6 -20 -20;
    -10 -10 -5 -5 -13 -13 -12 -11 -11 -10 -9 -8 -9 -9 -8 -8 -10 -10 -6 -6 -20 -20];
%Optimal stepsize for SG 1st row corresponds to sigma=10^-3; 2nd row
%corresponds to sigma=10^-5

steps_SS_lib=[-15 -14 -8 -8 -14 -14 -12 -12 -12 -11 -12 -11 -11 -11 -9 -9 -15 -15 -13 -13 -24 -24;
    -13 -13 -8 -8 -14 -14 -12 -12 -12 -11 -11 -11 -11 -11 -9 -9 -15 -15 -13 -13 -24 -24];
%Optimal stepsize for SS 1st row corresponds to sigma=10^-3; 2nd row
%corresponds to sigma=10^-5

for lib_ind=1:length(data_lib)
    variables_data(data_lib{lib_ind})=variable_lib{lib_ind};
    for sig_ind=1:length(sigma_lib)
        key_lib=strcat(data_lib(lib_ind),'-',num2str(sigma_lib{sig_ind},'%5.2e'));
        steps_SG_SS(key_lib{1})=[steps_SG_lib(sig_ind,lib_ind) steps_SS_lib(sig_ind,lib_ind)];
    end
end
%including one additional key for nonsmooth loss rand-50-50
variables_data('Rand-50-50')=50;
steps_SG_SS('Rand-50-50-')=[-14 -19];
%key_set=keys(steps_SG_SS);

% Plotting the results for given experimental setting
for dat=1:length(datas)
    data=datas{dat};
    variables=variables_data(data);
    for ik=1:length(sigmas)
        lambda=sigmas{ik};
        key=strcat(data,'-',num2str(lambda,'%5.2e'));
        if isKey(steps_SG_SS,key)
            steps=steps_SG_SS(key);
            alphaSG=steps(1);
            alphaSS=steps(2);
        else
            disp('optimal stepsizes for FD-SG and FD-SS algorithms are not given');
        end
        d=variables;
        SetupExperiment_DFO(data, lambda,d,alphaSG, alphaSS, theta,S0_SG, S0_Ada,xlimit,xlimitwide,gamma,method,rand_runs_adamethods);
    end
end

% %
% for dat=1:length(datas)
%     data=datas{dat};
%     for ik=1:length(sigmas)
%         lambda=sigmas{ik};
%         switch data
%             case {'15-absnormal','15-relnormal'}
%                 variables=30;
%                     
%             case {'18-absnormal','18-relnormal'}
%                 variables=11;
%             case {'19-absnormal-50','19-relnormal-50'}
%                 variables=50;
%             case {'20-absnormal','20-relnormal'}
%                 variables=20;
%             case {'22-absnormal','22-relnormal'}
%                 variables=8;
%             case 'Rand-50-50'
%                 variables=50;
%                 alphaSG=-14; alphaSS=-19; 
%             otherwise
%                 need to mention a stepsize for FD-SG and SS
%                 alphaSG=steps_SG(dat,ik);
%                 alphaSS=steps_SS(dat,ik);
%         end
%         d=variables; % number of variables in the problem
%         SetupExperiment_DFO(data, lambda,d,alphaSG, alphaSS, theta,S0_SG, S0_Ada,xlimit,xlimitwide,gamma,method)
%     end
% end


%% Old Code
%p = profile('info');
%save myprofiledata p
% %% 
% datas={'15-absnormal','20-absnormal', '15-relnormal','22-absnormal',...
%     '19-absnormal-50','18-absnormal','20-relnormal','22-relnormal','19-relnormal-50',...
%     '18-relnormal'};
% sigmas={10^-3};
% variables={30,20,30,8,50,11,20,8,50,11};
% theta=0.9;
% xlimit=1000000;
% 
% 
% %% Modified by raghu for now
% 
% expvec = [1];
% sigma=10^-5;
% steps={-10, -12, -10, -11,-13,-5,-11,-10,-13,-5}; % 10^-5
% stepsgs={-8, -9, -8, -9,-8,-4,-9,-9,-8,-5};% 10^-5
% steps_ss={-13, -12, -12, -12, -14, -8, -12, -11, -14, -8}; %for
% %sigma=10^-3;
% %steps={-10, -12, -10, -11,-10,-5,-10,-10,-12,-5}; % 10^-3
% %stepsgs={-10, -11, -9, -10,-7,-4,-8,-9,-8,-4}; % 10^-3
% S0_SG=2;
% S0_Ada=2;
% xlimit=5*10^5;
% xlimitwide=10000;
% variance='SAMPLE';
% gamma=0.9;
% method="-AD1";
% for i=expvec
%         data=datas{i};
%         lambda=sigma;
%         alphaSG=steps{i};
%         alphaSS=steps_ss{i};
%         d=variables{i};
%         SetupExperiment_DFO(data, lambda,d,alphaSG, alphaSS, theta,S0_SG, S0_Ada,xlimit,xlimitwide,gamma,method)
% end
% 
% sigma=10^-3;
% steps={-10, -12, -10, -11,-13,-5,-10,-10,-13,-5}; % 10^-3
% stepsgs={-10, -11, -9, -10,-8,-4,-8,-9,-8,-4}; % 10^-3
% 
% for i=expvec
%         data=datas{i};
%         lambda=sigma;
%         alphaSG=steps{i};
%         alphaGS=stepsgs{i};
%         d=variables{i};
%         %SetupExperiment_DFO(data, lambda,d,alphaSG, alphaGS, theta,S0_SG, S0_Ada,xlimit,xlimitwide,variance,gamma,method)
% end
% %%
% datas={'Rand-50-50','Eye-50-50'};
% variable={50,50};
% S0_SG=2;
% S0_Ada=2;
% steps={-14,-10};
% stepsgs={-13,-8};
% sigmas={[]};
% %loss='MAD';
% xlimitwide=10000;
% variance='SAMPLE';
% gamma=0.99;
% method="-AD1";
% loss='MAD-PFLoss';
% for i=1:length(datas)
%         data=datas{i};
%         lambda=sigmas{1};
%         alphaSG=steps{i};
%         alphaGS=stepsgs{i};
%         d=variables{i};
%         %SetupExperiment_DFO(data, lambda,d,alphaSG, alphaGS, theta,S0_SG, S0_Ada,xlimit,xlimitwide,variance,gamma,method)
% end
% 
% % %%sigma=10^-5;
% % steps={-10, -16, -10, -10}; % for sgd with S0=64
% % stepsgs={-9, -18, -8, -10};
% % S0_SG=2;
% % S0_Ada=64;
% % xlimit=5*10^5;
% % xlimitwide=100000;
% % for i=expvec
% %         data=datas{i};
% %         lambda=sigma;
% %         alphaSG=steps{i};
% %         alphaGS=stepsgs{i};
% %         d=variables{i};
% %         SetupExperiment1_AdaLBFGS_CuterDFO(data, lambda,d,alphaSG, alphaGS,theta,S0_SG, S0_Ada,xlimit,xlimitwide)
% % end
