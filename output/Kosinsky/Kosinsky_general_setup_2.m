%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning off
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_modified.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)
sensitivity = true;
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
    exclude_parameters = {'CD107a' 'CD8' 'K_D_antiPDL1' 'K_D_antiCTLA4'...
        'Total_Cell_Count' 'vol_Tcell' 'vol_Tumor', 'T_0' 'CD107a_logit'...
        'CD8_logit' 'k_LN' 'ka'};
    parameters = setdiff(parameters, exclude_parameters);
    parameters = [ parameters; 'T_0'];
else
    parameters = load(strjoin({cd,'parameters_hat_2.mat'},'/'));
    parameters = parameters.parameters_hat;
    parameters = ['E_0'; parameters; 'T_0'];

end
% Define outputs% Define outputs
observables={'TV' 'CD8_logit' 'CD107a_logit' 'DCm' 'ISC' 'PDL1'};

% Create PI with data
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1', ...
    'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1', 'MOC2_Control', 'MOC2_Control_Mean',...
    'MOC2_antiPDL1', 'MOC2_antiCTLA4'};

doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};

PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables, 'zeroAction', 'omit','mergePhenotypes', true);

[sim,u]=initializePI(model,parameters,observables,PI,doses, '');

PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []' 'Relative units []'...
    'Relative units []' 'Relative units []'};

normIndx = 4:6;
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1});
% Get initial values
x_0 = getInitialValues([PI.data(:).Group], initialStruct);
close all
%% Optimization setup
% Hierarchical structure
H = getHierarchicalStruct(parameters(1:end-1),'n_sigma', length(observables), 'n_rand', 1, 'n_indiv', length(u));
try
    sigmaNames=arrayfun(@(x)strjoin({'Omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false)';
    sigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
catch
    sigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
end
% Generating PI
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
     repelem(1, length([H.IndividualParams(:).Index]),1);...
    repelem(.1, length(H.SigmaParams),1)];
PI.par = getParamStruct2(sim,H,size(PI.data,1),repelem(0.5,length(H.SigmaParams),1),...
    sigmaNames,'Sigma', sigma_prior);
finalValues =log([PI.par(:).startValue]);

% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1),'initialValue',x_0),(@(x)getCovariance(x,H)),normIndx);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,H,'type','uniform'));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,length(u),'initialValue',x_0),exp(finalValues(end-length(observables)+1:end)),normIndx);


%% Save results
save('PI_Kosinsky_2.mat', 'PI')
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'