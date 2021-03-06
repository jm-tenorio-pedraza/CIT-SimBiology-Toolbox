%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
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
sensitivity = false;
%% load data and previous results

% Create function handle for simulations
% Define parameters to estimate
if sensitivity
     [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
     % Define parameters to estimate
    parameters=name(value>0);
    exclude_parameters = {'vol_Tumor' 'CTLA4_CD8_0' 'CTLA4_E'...
        'KD_antiCTLA4' 'KD_antiPDL1' 'koff_antiCTLA4' ...
        'koff_antiPDL1' 'CD107a' 'CD8_E' 'CD8_N' 'CD4Foxp3' 'GMDSC' 'DC' 'Q23' ...
        'PDL1_Immune_0' 'PDL1_Tumor_0' 'k23' 'k32' 'ka'...
        'ks_IFNg' 'TV' 'kdeg_CTLA4' 'kdeg_PDL1'...
        'K_IFNg' 'f1' 'f2' 'T_0' 'kdep'};
    uncertain_parameters = {'kin_DC' 'kin_MDSC'...
        'ks_PDL1_Immune' 'ks_PDL1_Tumor'};
    parameters = setdiff(parameters, exclude_parameters);
    parameters = [parameters; 'T_0'];
else
    parameters = {'kin_CD8';'K_CD8'; 'KDE_Treg';
        'KDE_MDSC';'kpro_Tumor_0'; 'kill_max'; 'K_el'; 'K_DC';
        'kin_Treg' ; 'K_IFNg';'K_MDSC';'kin_DC';'kin_MDSC'};
    parameters = [parameters; 'T_0'];

end
% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1', ...
    'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1', 'MOC2_Control', 'MOC2_Control_Mean',...
    'MOC2_antiPDL1', 'MOC2_antiCTLA4'};
observables={'TV' 'CD8_logit' 'CD107a_logit' 'Treg_logit' 'DC_logit' 'MDSC_logit' 'PDL1_T' 'PDL1_I'};
stateVar={'Tumor' 'CD8_logit' 'CD107a_logit' 'Treg_logit' 'DC_logit' 'GMDSC_logit'...
    'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    true,'output', 'mean','maxIIV', false);
PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []'  'Logit []' ...
     'Logit [t]'   'Logit []' ...
    'Relative units []' 'Relative units []'};
variants = getvariant(model);

PI.normIndx = 7:8;
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1},...
    'variant', {variants(1); variants(2)});
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1_Optimized','doseUnits', 'mole');

%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', 2, 'cell_indx',[5 6 8 11], 'n_indiv', length(PI.u));
if ~isempty(PI.H.IndividualParams(1).Index)
        indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),PI.H.IndividualParams,'UniformOutput',false)';
else
    indivSigmaNames = [];
end
if ~isempty(PI.H.CellParams(1).Index)
    cellSigmaNames=arrayfun(@(x)strjoin({'lambda', x.name}, '_'),PI.H.CellParams,'UniformOutput',false)';
else
    cellSigmaNames = [];
end
try
SigmaNames = [cellSigmaNames; indivSigmaNames];
SigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'sigma', x}, '_'),...
    observables','UniformOutput', false);
catch
    SigmaNames=cellfun(@(x) strjoin({'sigma', x}, '_'),...
    observables','UniformOutput', false);
end
% Generating PI
alpha = [repelem(0.01, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.01, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
beta = [repelem(0.1, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.1, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
sigma_prior= [ repelem(1,length(PI.H.PopulationParams), 1);...
    repelem(1, length([PI.H.CellParams(:).Index]),1);
     repelem(1, length([PI.H.IndividualParams(:).Index]),1);...
    alpha];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior);
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/uniform'}));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),exp(finalValues(end-length(observables)+1:end)),PI.normIndx);
% Parameter names for plots
paramNames = ['kin_{CD8}' '\eta_{K_{CD8}}' '\eta_{kpro_{Tumor}}' 'kill_{max}' 'KDE_{Treg}'...
    'KDE_{MDSC}' 'K_{pro}' 'K_{el}' 'kpro_{Tumor_{MOC1}}' 'kpro_{Tumor_{MOC2}}'...
    {PI.par([PI.H.IndividualParams(:).Index]).name} '\lambda_{kpro_{Tumor}}'...
    '\omega_{kin_{CD8}}' '\sigma_{TV}' '\sigma_{CD8}' '\sigma_{CD107a}' '\sigma_{Treg}'...
    '\sigma_{DC}' '\sigma_{MDSC}' '\sigma_{PDL1_T}' '\sigma_{PDL1_I}'];

%% Save results
save('PI_CIM_2_1.mat', 'PI')
load(strjoin({cd 'PI_CIM_1.mat'},'/'))

load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'