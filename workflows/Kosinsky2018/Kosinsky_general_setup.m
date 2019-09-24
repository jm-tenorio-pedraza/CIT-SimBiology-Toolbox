%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning off
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kosinsky/output/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)

%% load data and previous results
stateVar={'Tumor' 'CD8_logit' 'CD107a_logit' 'DC_Rel' 'GMDSC_Rel'...
    'Tumor_PDL1_Rel'};

% Create function handle for simulations
% Define parameters to estimate
parameters = load('Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/parameters_hat.mat');
parameters = parameters.parameters_hat;
% Define outputs% Define outputs
observables={'TV' 'CD8_logit' 'CD107a_logit' 'DCm' 'ISC' 'PDL1'};

% Create PI with data
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1'};

doses = {'Dose_antiPDL1'};

PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',stateVar,groups_subset,observables);

[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1');

PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []' 'Relative units []'...
    'Relative units []' 'Relative units []'};


%% Optimization setup
% Hierarchical structure
H.PopulationParams=1:length(parameters);
H.IndividualParams=struct('name', {}, 'Index', { }, 'EtaIndex', {}, 'OmegaIndex', {});
H.SigmaParams=length(parameters)+1:length(parameters)+length(observables);

% Generating PI
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))),{PI.data.dataValue},'UniformOutput',true));
sigmaNames={'b_TV'; 'b_CD8'; 'b_CD107a'; 'b_DC'; 'b_MDSC'; 'b_PDL1'};
sigma_prior=[repelem(2,length(parameters),1);repelem(1,length(observables),1)];
PI.par = getParamStruct2(sim,H,13,repelem(0.5,6,1),sigmaNames,'Sigma', sigma_prior);

% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),4:6);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution3(exp(p),PI,H,'type','uniform');

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));


%% Save results
save('PI.mat', 'PI')
