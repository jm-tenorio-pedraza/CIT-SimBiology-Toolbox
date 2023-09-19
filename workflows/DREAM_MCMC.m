%% DREAM MCMC
clear Cell_Field clavijo_counts Clavijo_Names clavijo_unique clavijoData counts data_ext1 ImmuneResp morisada_counts Morisada_Names morisada_unique morisadaData ans burnIn data_ext delta groups_subset immuneResp  doses  groups_subset2 parameters observables stateVar beta cellParamsIndx finalValue i fval_minunc indivParamsIndx lb ncol nrow options_fminsearch p_hat p_hat_cell p_hat_pop popParamsIndx residuals_cell residuals_fn residuals_pop sigma_prior SigmaNames simTime ub 
%% Initial vectors
% finalValues = log([PI.par(:).finalValue]);
d = length(finalValues);
N_samples=26;
% LHS design for prior sampling
thetaLHS_0 = lhsdesign(N_samples,length(PI.H.PopulationParams));
lb=paramValues'*.1;
ub=paramValues'*10;
thetaLHS=thetaLHS_0.*(ub-lb)+lb;
% Generate random noise for var params
Z=randn(N_samples,d-length(lb))*.1;
theta= [log(thetaLHS), Z+repmat(finalValues(length(PI.H.PopulationParams)+1:end),N_samples,1)];
X0 = [finalValues; (theta)];
% Evaluate loglikelihood 
logL=rowfun(obj_fun,table(X0));
[L,I]=sort(logL{:,:});
w0 = X0(I(1:d*2),:);
% priorCov= cov(w0);
% K = chol(priorCov)';
Delta = getGradient(finalValues, obj_fun,'parallel',false,'delta',.001);
epsilon = 1e-3;
X1 =[ finalValues; finalValues+epsilon*(1e2*randn(1e3,length(finalValues))-epsilon*Delta)];
logL1 = rowfun(obj_fun,table(X1));
[L1,I1]=sort(logL1{:,:});
w0 = X1(I1(1:d),:);

%%
i =1;
if i==1
    J=ones(1,3);
    n_id=ones(1,3);
    stepSize=2.38*.25;
else
    J = J1;
    n_id=n_id1;
    stepSize=stepSize1;
    w0 = x1(:,:,end)';
end
%%
% p=parpool('local');
% pctRunOnAll warning off

% tic
% [x, p_x,stepSize, J, n_id,accept] = par_dream_1(w0,prior_fun, likelihood_fun,...
%     N,ceil(.2e4/N), d, 'burnIn',1e6,'stepSize',stepSize,'J',J,'n_id',n_id, 'plotLogL', true,'H', PI.H);
% toc
%% Short run to find MLE
N=3;
T = ceil(.9e5/N);
tic
[x, p_x, stepSize, J, n_id,accept,Z] = par_dream_zs(Z,prior_fun, likelihood_fun,...
    N,T,d, 'burnIn',0,'J',J,'n_id',n_id,'stepSize', stepSize,'plotLogL',true);
toc

%%
x3 =x;
p_x3 =p_x;
%% Save results
N_i='3';  
save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_Z' '.mat'},''), 'Z')
%  save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_w0' '.mat'},''), 'w0')

save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_x_' N_i '.mat'},''), strjoin({'x' N_i},''))
save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_p_x_' N_i '.mat'},''), strjoin({'p_x' N_i},''))
save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_J' '.mat'},''), strjoin({'J' ''},''))
save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_n_id' '.mat'},''), strjoin({'n_id' ''},''))
save(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_stepSize' '.mat'},''), strjoin({'stepSize' ''},''))
%%
load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  PI.model '_Z'  '.mat'},''))
load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  'CIM10_PKPD_MOC1_kin_CD8' '_w0'  '.mat'},''))

load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  PI.model '_x_' N_i '.mat'},''))
load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  PI.model '_p_x_' N_i '.mat'},''))
load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  PI.model '_J' '.mat'},''))
load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  PI.model '_n_id'  '.mat'},''))
load(strjoin({cd '/CIM/PI/CIM01/MCMC/'  PI.model '_stepSize' '.mat'},''))

for i=1:2
    load(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_x_'  num2str(i) '.mat'},''))
    load(strjoin({cd '/CIM/PI/CIM01/MCMC/' PI.model '_p_x_' num2str(i) '.mat'},''))
end
x = cat(3,x1,x2,x3);
p_x=cat(1,p_x1,p_x2,p_x3);
%% clear variables
clearvars -except d J likelihood_fun N obj_fun p prior_fun sim w0 x stepSize n_id paramNames PI finalValues w0 z0 Z