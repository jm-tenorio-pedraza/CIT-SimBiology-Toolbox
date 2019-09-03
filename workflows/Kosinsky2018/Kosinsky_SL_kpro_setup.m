%% Search paths
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky/output/PI_SL_kpro')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky.sbproj');

% Extract model
kosinsky=out.m1;
cs=kosinsky.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)

%% load data and previous results
data_subset=Kosinsky_data([1:2:23 24],:);
PI.data=data_subset;
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]' 'Relative units []'...
    'Relative units []' 'Relative units []'};
load('Kosinsky_data.mat')
load('models_SL_kpro.mat')

%% Create function handle for simulations
% Define parameters to estimate
parameters={'r' 'TV_max' 'e_Td' 'k_LN' 'K_TCD' 'K_pdl' 'S_R' 'S_L' 'k_pro'};

% Define outputs
observables={'TV' 'CD8' 'CD107a' 'DCm' 'ISC' 'PDL1'};
[sim,u]=initializePI(kosinsky,parameters,observables,PI);

%% Optimization setup
% Hierarchical structure
H.PopulationParams=1:7;
H.IndividualParams=struct('name', {'S_L'; 'k_pro'}, 'Index', {10:22; 23:35},...
    'EtaIndex', {8;9}, 'OmegaIndex', {36;37});
H.SigmaParams=38:43;

% Generating PI
% Generating PI
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
sigmaNames={'omega_SL'; 'omega_kpro'; 'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'};
sigma_prior=[1 2 2 2 1 2 2 2 2 repelem(2,26) 1 1 1 1 1 1 1 1]';
PI.par = getParamStruct2(sim,H,13,repelem(0.5,8,1),sigmaNames,'Sigma', sigma_prior);

% Residuals function
residuals_fun=@(p)getResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),4:6);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution3(exp(p),PI,H,'type','normal');

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));

%% MCMC setup
postSample=models_array(:,:,1.8e4:1000:end);
postSample=postSample(:,:)';

w0 = [postSample(:, H.PopulationParams), mean(postSample(:, 8:20),2),mean(postSample(:, 21:33),2),...
    postSample(:,8:33), std(postSample(:,8:20),0,2), std(postSample(:,21:33),0,2), postSample(:,34:39)];

%% Save results
save('PI_SL_kpro.mat', 'PI')
save('models_SL_kpro.mat','models')
save('logP_SL_kpro.mat','logP')
