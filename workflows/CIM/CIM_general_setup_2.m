%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning off
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_2.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)
sensitivity = true;
%% load data and previous results

% Create function handle for simulations
% Define parameters to estimate
if sensitivity
     [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
     % Define parameters to estimate
    parameters=name(value>0);
    exclude_parameters = {'Avogadro' 'omega' 'vol_Tumor' 'CD25_0' 'CTLA4_CD8_0'...
        'CTLA4_Treg_0' 'KD_CD25' 'KD_antiCTLA4' 'KD_antiPDL1' 'koff_CD25' 'koff_antiCTLA4' ...
        'koff_antiPDL1' 'CD107a' 'CD8_E' 'CD8_N' 'CD4Foxp3' 'Q23' ...
        'PDL1_Immune_0' 'PDL1_Tumor_0' 'k23' 'k32' 'ka'...
        'ks_IFNg' 'Tumor_P' 'kdeg_CTLA4' 'kdeg_PDL1'...
        'K_IFNg' 'f1' 'f2' 'f3' 'T_0'};
    uncertain_parameters = {'kin_DC' 'kin_MDSC'...
        'ks_PDL1_Immune' 'ks_PDL1_Tumor'};
    parameters = setdiff(parameters, [exclude_parameters uncertain_parameters ]);
    parameters = [parameters; 'T_0'];
else
    parameters = load(strjoin({cd,'parameters_hat_2.mat'},'/'));
    parameters = parameters.parameters_hat;
    parameters = [parameters; 'T_0'];

end
% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1', ...
    'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1', 'MOC2_Control', 'MOC2_Control_Mean',...
    'MOC2_antiPDL1', 'MOC2_antiCTLA4'};
observables={'TV' 'CD8_logit' 'CD107a_logit' 'Treg_logit' 'DC_logit' 'MDSC_logit' 'PDL1_Tumor_Rel' 'PDL1_Immune_Rel'};
stateVar={'Tumor' 'CD8_logit' 'CD107a_logit' 'Treg_logit' 'DC_logit' 'GMDSC_logit'...
    'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'input','mergePhenotypes', true);
PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []'  'Logit []' ...
     'Logit [t]'   'Logit []' ...
    'Relative units []' 'Relative units []'};

% Get simulation function
[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1', 'doseUnits', 'mole');
close all


normIndx = 7:8;
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1});
% Get initial values
x_0 = getInitialValues([PI.data(:).Group], initialStruct);
close all
%% Optimization setup
% Hierarchical structure
H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', 23, 'cell_indx', 24, 'n_indiv', length(u));
try
    cellSigmaNames=arrayfun(@(x)strjoin({'lambda', x.name}, '_'),H.CellParams,'UniformOutput',false)';
    indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false)';
    SigmaNames = [cellSigmaNames; indivSigmaNames];
    SigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'sigma', x}, '_'),...
        observables,'UniformOutput', false);
catch
    SigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables','UniformOutput', false);
end
% Generating PI
alpha = [repelem(0.01, length([H.CellParams(:).OmegaIndex]),1);...
    repelem(0.01, length([H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(H.SigmaParams, [H.CellParams(:).OmegaIndex ...
    H.IndividualParams(:).OmegaIndex])),1)];
beta = [repelem(0.1, length([H.CellParams(:).OmegaIndex]),1);...
    repelem(0.1, length([H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(H.SigmaParams, [H.CellParams(:).OmegaIndex ...
    H.IndividualParams(:).OmegaIndex])),1)];
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
    repelem(1, length([H.CellParams(:).Index]),1);
     repelem(1, length([H.IndividualParams(:).Index]),1);...
    alpha];
PI.par = getParamStruct2(sim,H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior);
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,PI.tspan(end),u,PI.tspan),PI,...
    @(x)getPhi2(x,H,size(PI.data,1),'initialValue',x_0),(@(x)getCovariance(x,H)),normIndx);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,H,'type','uniform'));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),u,PI.tspan),PI,...
    @(x)getPhi2(x,H,length(u),'initialValue',x_0),exp(finalValues(end-length(observables)+1:end)),normIndx);

%% Save results
save('PI_CIM_2.mat', 'PI')
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'