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
% All parameters 
[name,I] = sort(get(model.Parameters, 'Name'));
value = cell2mat(get(model.Parameters, 'Value'));
value = value(I);
% Define parameters to estimate
parameters=name(value>0);
exclude_parameters = {'CD107a' 'CD8' 'K_D_antiPDL1' 'Total_Cell_Count' 'vol_Tcell' 'vol_Tumor', 'T_0'};
parameters = setdiff(parameters, exclude_parameters);

% Define outputs
observables={'TV' 'CD8_logit' 'CD107a_logit' 'DCm' 'ISC' 'PDL1'};

% Create PI with data
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1'};

doses = {'Dose_antiPDL1'};
PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',stateVar,groups_subset,observables);
dose = {'Dose_antiPDL1'};

[sim,u]=initializePI(model,parameters,observables,PI,doses, 'CT26');

PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []' 'Relative units []'...
    'Relative units []' 'Relative units []'};
%% Save output
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/phat_Control.mat','par_hat')
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/parameters_hat.mat', 'parameters_hat')

%% Load results
phat_Control = load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/phat_Control.mat','par_hat');
