%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/CIM/output/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_2.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)

%% Create function handle for simulations
% Define parameters to 
if sensitivity
    [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
    % Define parameters to estimate
    parameters=name(value>0);
    exclude_parameters = {'Avogadro' 'omega' 'vol_Tumor' 'CD25_0' 'CTLA4_CD8_0'...
        'CTLA4_Treg_0' 'KD_CD25' 'KD_antiCTLA4' 'KD_antiPDL1' 'koff_CD25' 'koff_antiCTLA4' ...
        'koff_antiPDL1' 'CD107a' 'CD8_E' 'CD8_N'};
    fixed_parameters = {'PDL1_Immune_0' 'PDL1_Tumor_0'...
        'k12' 'k23' 'k32' 'ka' 'ks_IFNg' 'ks_IL2' 'Tumor_P' 'Tumor_P_0'...
        'kdeg_CTLA4' 'kdeg_PDL1' 'kel_Debris'  'K_IFNg' ...
         'f1' 'f2' 'f3'};
    uncertain_parameters = {'K_IFNg' 'K_IL2'  'kpro_Naive_max' ...
        'kel_DC'  'kel_MDSC'  'kel_Tumor'};
    parameters = setdiff(parameters, [exclude_parameters fixed_parameters ]);
else
    parameters = load('Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/parameters_hat.mat');
    parameters = parameters.parameters_hat;
end
% Define outputs
observables={'TV' 'CD8_logit' 'CD107a_logit' 'Treg_logit' 'DC_logit' 'MDSC_logit' 'PDL1_Tumor_Rel' 'PDL1_Immune_Rel'};
stateVar={'Tumor' 'CD8_logit' 'CD107a_logit' 'Treg_logit' 'DC_logit' 'GMDSC_logit'...
    'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
groups_subset = {'MOC1_Control' 'MOC1_antiPDL1' 'MOC1_Control_Mean', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
doses = {'Dose_antiPDL1', 'Dose_antiCTLA4'};

% Obtain PI
PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',stateVar,groups_subset,observables);
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Percentage [%]' ...
    'Relative units []' 'Relative units []'};

% Get simulation function
[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1');
close all
%% Optimization setup
% Hierarchical structure
H = getHierarchicalStruct(parameters,'n_sigma', length(observables), 'n_rand', 1, 'n_indiv', length(u));
try
    sigmaNames=arrayfun(@(x)strjoin({'Omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false);
    sigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
catch
    sigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
end
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
    repelem(1, length(H.IndividualParams.Index),1);...
    repelem(0.5, length(H.SigmaParams),1)];
PI.par = getParamStruct2(sim,H,7,repelem(0.5,length(H.SigmaParams),1),sigmaNames,'Sigma', sigma_prior);

% load('PI_CIM.mat')

% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),7:8);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)createPriorDistribution3(exp(p),PI,H,'type','uniform');

residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,length(u)),exp(finalValues(end-length(observables)+1:end)),7:8);


%% Save results
save('PI_CIM.mat', 'PI')
load('PI_CIM.mat')
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))
%% Save output

save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/parameters_hat.mat', 'parameters_hat')
