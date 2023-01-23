function [iter, batch_stats, funvalsiter_stats, funvalsepoch_stats,Epoch,steps_stats, ExpBudgetMarkers] = Output_Stats(func, data,Options,folder,seeds,d,interval)
%This functions caluclates the min, max, mean, median for multiple runs of
%algorithms.

% Load all the multiple runs
r=length(seeds);
% initializing the libraries
iter_lib={};batch_lib={};funvals_lib={};Evals_lib={};FEvals_lib={};steps_lib={};FuncEvals_lib={};
min_iter=10^16;
minEpoch=10^16;
maxEpoch=10^16;
% Note that this is not the iteration number for SGD - We have to correct
% it as SGD doesn't store values for each iteration
for run=1:r
    [iter_lib{run},batch_lib{run},funvals_lib{run},Evals_lib{run},FEvals_lib{run},steps_lib{run}]=func(data,Options,folder,run);
    min_iter=min(min_iter,length(iter_lib{run}));
    % Computing total work: function evaluations used in line search + gradient
    % estimates
    switch Options.Method
        case {'Norm-AD1','IPQN-AD1'}
            % correcting the total work in the post processing step for full work 
            FuncEvals_lib{run}=[0; (Evals_lib{run}(1:end-1) + batch_lib{run}(1:end-1))*(d+1) + FEvals_lib{run}(2:end)];
            %FuncEvals={[0; (Evals1(1:end-1) + batch1(1:end-1))*(d+1) + FEvals1(2:end)], [0; (Evals2(1:end-1) + batch2(1:end-1))*(d+1) + FEvals2(2:end)],[0; (Evals3(1:end-1) + batch3(1:end-1))*(d+1) + FEvals3(2:end)],[0; (Evals4(1:end-1) + batch4(1:end-1))*(d+1) + FEvals4(2:end)],[0; (Evals5(1:end-1) + batch5(1:end-1))*(d+1) + FEvals5(2:end)]};
        otherwise
            FuncEvals_lib{run}=[0; (Evals_lib{run}(1:end-1))*(d+1) + FEvals_lib{run}(2:end)];
            %FuncEvals={[0; (Evals1(1:end-1))*(d+1) + FEvals1(2:end)], [0; (Evals2(1:end-1))*(d+1) + FEvals2(2:end)],[0; (Evals3(1:end-1))*(d+1) + FEvals3(2:end)],[0; (Evals4(1:end-1))*(d+1) + FEvals4(2:end)],[0; (Evals5(1:end-1))*(d+1) + FEvals5(2:end)]};
    end
    minEpoch=min(minEpoch,min(FuncEvals_lib{run}));
    maxEpoch=min(maxEpoch, max(FuncEvals_lib{run}));
end
iter=1:min_iter;
batch=zeros(min_iter,r); steps=zeros(min_iter,r); funvals=zeros(min_iter,r);
budget=zeros(min_iter,r);
for run=1:r
    batch(:,run)=batch_lib{run}(iter); steps(:,run)=steps_lib{run}(iter); funvals(:,run)=funvals_lib{run}(iter);
    budget(:,run)=FuncEvals_lib{run}(iter);
end
Evals=Evals_lib;
FEvals=FEvals_lib;
FuncEvals=FuncEvals_lib;

% Computing data at regular epoch intervals
Epoch=minEpoch:interval:maxEpoch;
funvalsepoch=zeros(length(Epoch),r);
funvalscell=funvals_lib;
ind=ones(r,1);
for i=1:length(Epoch)
    eval=Epoch(i);
    for j=1:length(ind)
        if FuncEvals{j}(ind(j)) >= eval
            funvalsepoch(i,j)=funvalscell{j}(ind(j));
        else
            while FuncEvals{j}(ind(j)) < eval
                ind(j)=ind(j) + 1;
            end
            funvalsepoch(i,j)=funvalscell{j}(ind(j)-1);
        end
    end 
end


% Evals=[Evals1(iter) Evals2(iter) Evals3(iter) Evals4(iter) Evals5(iter)];
% FEvals=[FEvals1(iter) FEvals2(iter) FEvals3(iter) FEvals4(iter) FEvals5(iter)];
% FuncEvals=Evals*(d+1) + FEvals;
% Epoch=min(FuncEvals(:)):interval:max(FuncEvals(:));
% funvalsepoch=zeros(length(Epoch),r);
% ind=ones(r,1);
% for i=1:length(Epoch)
%     eval=Epoch(i);
%     for j=1:length(ind)
%         if FuncEvals(ind(j)) >= eval
%             funvalsepoch(i,j)=funvals(ind(j),j);
%         else
%             while FuncEvals(ind(j)) < eval
%                 ind(j)=ind(j) + 1;
%             end
%             funvalsepoch(i,j)=funvals(ind(j)-1,j);
%         end
%     end 
% end


% Computing expected budget at each iteration for including some markers
% around that point. working on the code
% budget=[FuncEvals{1}(iter) FuncEvals{2}(iter) FuncEvals{3}(iter) FuncEvals{4}(iter) FuncEvals{5}(iter)];
budget_mean=sum(budget,2)/r;
ExpBudgetMarkers=zeros(length(budget_mean),1);
ep_ind=1;
for i=1:length(budget_mean)
    while Epoch(ep_ind) < budget_mean(i)
        ep_ind=ep_ind + 1;
        if ep_ind > length(Epoch)
            ep_ind=ep_ind-1;
            break;
        end
    end
    ExpBudgetMarkers(i)=ep_ind;
end

batch_min=min(batch,[],2);
batch_max=max(batch,[],2);
batch_mean=sum(batch,2)/r;
batch_med=median(batch,2);
batch_stats=[batch_min batch_max batch_mean batch_med];

steps_min=min(steps,[],2);
steps_max=max(steps,[],2);
steps_mean=sum(steps,2)/5;
steps_med=median(steps,2);
steps_stats=[steps_min steps_max steps_mean steps_med];

funvalsiter_min=min(funvals,[],2);
funvalsiter_max=max(funvals,[],2);
funvalsiter_mean=sum(funvals,2)/r;
funvalsiter_med=median(funvals,2);
funvalsiter_stats=[funvalsiter_min funvalsiter_max funvalsiter_mean funvalsiter_med];

funvalsepoch_min=min(funvalsepoch,[],2);
funvalsepoch_max=max(funvalsepoch,[],2);
funvalsepoch_mean=sum(funvalsepoch,2)/r;
funvalsepoch_med=median(funvalsepoch,2);
funvalsepoch_stats=[funvalsepoch_min funvalsepoch_max funvalsepoch_mean funvalsepoch_med];

% [iter31,batch31,funvals31,Evals31,FEvals31,steps31]=CallPlotExperiments_AdaLBFGS_CuterDFO(data,Options,'Expmnt1-DFO',1);
% [iter32,batch32,funvals32,Evals32,FEvals32,steps32]=CallPlotExperiments_AdaLBFGS_CuterDFO(data,Options,'Expmnt1-DFO',2);
% [iter33,batch33,funvals33,Evals33,FEvals33,steps33]=CallPlotExperiments_AdaLBFGS_CuterDFO(data,Options,'Expmnt1-DFO',3);
% [iter34,batch34,funvals34,Evals34,FEvals34,steps34]=CallPlotExperiments_AdaLBFGS_CuterDFO(data,Options,'Expmnt1-DFO',4);
% [iter35,batch35,funvals35,Evals35,FEvals35,steps35]=CallPlotExperiments_AdaLBFGS_CuterDFO(data,Options,'Expmnt1-DFO',5);
% 
% iter3=1:min([max(iter31) max(iter32) max(iter33) max(iter34) max(iter35)]);
% batch_dss=[batch31(iter3) batch32(iter3) batch33(iter3) batch34(iter3) batch35(iter3)];
% steps_dss=[steps31(iter3) steps32(iter3) steps33(iter3) steps34(iter3) steps35(iter3)];
% funvals_dss=[funvals31(iter3) funvals32(iter3) funvals33(iter3) funvals34(iter3) funvals35(iter3)];
% Evals_dss=[Evals31(iter31) Evals32(iter31) Evals33(iter31) Evals34(iter31) Evals35(iter31)];
% FEvals_dss=[FEvals31(iter31) FEvals32(iter31) FEvals33(iter31) FEvals34(iter31) FEvals35(iter31)];
% FuncEvals_dss=Evals_dss*(d+1) + FEvals_dss;
% Epoch_dss=min(FuncEvals_dss(:)):d+1:max(FuncEvals_dss);
% funvalsepoch_dss=zeros(length(Epoch_dss),5);
% ind=ones(5,1);
% for i=1:length(Epoch_dss)
%     eval=Epoch_dss(i);
%     for j=1:length(ind)
%         if FuncEvals_dss(ind(j)) >= eval
%             funvalsepoch_dss(i,j)=funvals_dss(ind(j));
%         else
%             while FuncEvals_dss(ind(j)) < eval
%                 ind(j)=ind(j) + 1;
%             end
%             funvalsepoch_dss(i,j)=funvals_dss(ind(j)-1);
%         end
%     end 
% end
% 
% batch3_min=min(batch_dss,2);
% batch3_max=max(batch_dss,2);
% batch3_mean=sum(batch_dss,2)/5;
% batch3_med=median(batch_dss);
% 
% steps3_min=min(steps_dss,2);
% steps3_max=max(steps_dss,2);
% steps3_mean=sum(steps_dss,2)/5;
% steps3_med=median(steps_dss);
% 
% funvalsiter3_min=min(funvals_dss,2);
% funvalsiter3_max=max(funvals_dss,2);
% funvalsiter3_mean=sum(funvals_dss,2)/5;
% funvalsiter3_med=median(funvals_dss);
% 
% funvalsepoch3_min=min(funvalsepoch_dss,2);
% funvalsepoch3_max=max(funvalsepoch_dss,2);
% funvalsepoch3_mean=sum(funvalsepoch_dss,2)/5;
% funvalsepoch3_med=median(funvalsepoch_dss);
 end

