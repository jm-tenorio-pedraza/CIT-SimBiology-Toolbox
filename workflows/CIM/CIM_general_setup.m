%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM/output/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)

%% Create function handle for simulations
% Define parameters to estimate
parameters = load('Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/parameters_hat.mat');
parameters = parameters.parameters_hat;
% Define outputs
observables={'TV' 'CD8' 'CD107a' 'Treg' 'DCm' 'ISC' 'PDL1_Tumor_Rel' 'PDL1_Immune_Rel'};
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Percentage [%]' ...
    'Relative units []' 'Relative units []'};
stateVar={'Tumor' 'CD8' 'CD107a' 'Treg' 'DC' 'GMDSC'...
    'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
groups_subset = {'MOC1_Control' 'MOC1_antiPDL1' 'MOC1_Control_Mean', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
doses = {'Dose_antiPDL1', 'Dose_antiCTLA4'};
run('Clavijo_Group_Pre_Processing.m')
[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1');

%% Optimization setup
% Hierarchical structure
H.PopulationParams=1:length(parameters);
H.IndividualParams=struct('name', [],'Index', [], 'EtaIndex', [], 'OmegaIndex', []);
H.SigmaParams=length(parameters)+1:length(parameters)+length(observables);

% Generating PI
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
sigmaNames={ 'b_TV'; 'b_CD8'; 'b_CD107a';'b_Treg'; 'b_DC'; 'b_MDSC'; 'b_PDL1_Tumor'; 'b_PDL1_Immune'};
sigma_prior= [ repelem(2,length(H.PopulationParams), 1); repelem(1, length(H.SigmaParams),1)];
PI.par = getParamStruct2(sim,H,1,repelem(0.5,length(H.SigmaParams),1),sigmaNames,'Sigma', sigma_prior);

% load('PI_CIM.mat')

% Residuals function
residuals_fun=@(p)getResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),7:8);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution3(exp(p),PI,H,'type','lognormal');

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));
tic
obj_fun(log([PI.par(:).startValue]))
toc
%% MCMC setup
postSample=models_array(:,:,1.8e4:1000:end);
postSample=postSample(:,:)';

w0 = [postSample(:, H.PopulationParams), mean(postSample(:, 9:21),2),...
    postSample(:,9:21), std(postSample(:,9:21),0,2), postSample(:,22:27)];


%% Save results
save('PI_CIM.mat', 'PI')
