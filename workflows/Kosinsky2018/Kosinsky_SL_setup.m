%% Search paths
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky/output/PI_SL')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky.sbproj');

% Extract model
kosinsky=out.m1;
cs=kosinsky.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)

%% load data and previous results
load('Kosinsky_data.mat')
data_subset=Kosinsky_data([1:2:23 24],:);
PI.data=data_subset;
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]' 'Relative units []'...
    'Relative units []' 'Relative units []'};

load('models_Kosinsky_MOC1_SL.mat')
load('PI_SL.mat')

%% Create function handle for simulations
% Define parameters to estimate
parameters={'r' 'TV_max' 'e_Td' 'k_LN' 'K_TCD' 'K_pdl' 'S_R' 'k_pro' 'S_L'};

% Define outputs
observables={'TV' 'CD8' 'CD107a' 'DCm' 'ISC' 'PDL1'};

[sim,u]=initializePI(kosinsky,parameters,observables,PI);
%% Optimization setup
% Hierarchical structure
H.PopulationParams=1:9;
H.IndividualParams=struct('name', {'S_L'}, 'Index', {10:22}, 'EtaIndex', 9, 'OmegaIndex', 23);
H.SigmaParams=23:29;

% Generating PI
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
sigmaNames={'omega_kpro'; 'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'};
sigma_prior=[1 2 2 2 1 2 2 2 2 repelem(2,13) 1 1 1 1 1 1 1]';
PI.par = getParamStruct2(sim,H,13,repelem(0.5,7,1),sigmaNames,'Sigma', sigma_prior);

% Residuals function
residuals_fun=@(p)getResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),4:6);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution3(exp(p),PI,H,'type','lognormal');

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));
%% MCMC setup

postSample=models_array(:,:,1.8e4:1000:end);
postSample=postSample(:,:)';

w0 = [postSample(:, H.PopulationParams), mean(postSample(:, 9:21),2),...
    postSample(:,9:21), std(postSample(:,9:21),0,2), postSample(:,22:27)];


%% Save results
save('PI_SL.mat', 'PI')
save('models_Kosinsky_MOC1_SL.mat','models_array')
save('logP_Kosinsky_MOC1_SL.mat','logL')
