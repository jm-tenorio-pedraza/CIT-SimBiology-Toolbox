%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning off
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj');
% Extract model
model=out.model;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)

%% load data and previous results
stateVar={'Tumor' 'CD8_logit' 'CD107a_logit' 'DC_Rel' 'GMDSC_Rel'...
    'Tumor_PDL1_Rel'};

% Create function handle for simulations
% Define parameters to estimate
if sensitivity
    [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
    % Define parameters to estimate
    parameters=name(value>0);
    exclude_parameters = {'CD107a' 'CD8' 'K_D_antiPDL1' 'Total_Cell_Count' 'vol_Tcell' 'vol_Tumor', 'T_0'};
    parameters = setdiff(parameters, exclude_parameters);
else
    parameters = load(strjoin({cd,'parameters_hat.mat'},'/'));
    parameters = parameters.parameters_hat;
end
% Define outputs% Define outputs
observables={'TV' 'CD8_logit' 'CD107a_logit' 'DCm' 'ISC' 'PDL1'};

% Create PI with data
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1'};

doses = {'Dose_antiPDL1'};

PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',stateVar,groups_subset,observables);

[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1_optimized');

PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []' 'Relative units []'...
    'Relative units []' 'Relative units []'};

normIndx = 4:6;
%% Optimization setup
% Hierarchical structure
H = getHierarchicalStruct(parameters,'n_sigma', length(observables), 'n_rand', 1, 'n_indiv', length(u));
try
    sigmaNames={arrayfun(@(x)strjoin({'Omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false)};
    sigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
catch
    sigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
end
% Generating PI
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
     repelem(1, length([H.IndividualParams(:).Index]),1);...
    repelem(1, length(H.SigmaParams),1)];
PI.par = getParamStruct2(sim,H,size(PI.data,1),repelem(2,length(H.SigmaParams),1),...
    sigmaNames,'Sigma', sigma_prior);

% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),normIndx);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution3(exp(p),PI,H,'type','uniform');

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x')*(-1));

finalValues = log([PI.par(:).startValue]);



%% Save results
save('PI_Kosinsky.mat', 'PI')
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'