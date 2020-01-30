%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM3')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_3.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1},...
    'variant', {variants(1); variants(2)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-12);
set(cs, 'MaximumWallClock', 0.25)
%% Parameter setup
parameters = {'K_CTLA4'; 'kin_CD8';'kill_CD8';...
    'kpro_Tumor';'KDE_MDSC'; 'kin_DC'; 'K_CD8'; 'kill_Treg';...
    'K_IFNg'; 'K_PDL1'};
sensitivity_parameters = {};

parameters = [parameters; 'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC2_Control',...
    'MOC1_antiPDL1','MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1',...
    'MOC2_antiPDL1','MOC2_antiCTLA4', 'MOC2_antiCTLA4_antiPDL1'};
observables={'TV'};
stateVar={'Tumor'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    false,'output', 'mean','maxIIV', true);
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'CIM ICB (K_CD8)';
PI.observablesPlot={'TV'};

% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1',...
    'doseUnits', 'mole', 'parallel', false);

%% Optimization setup
% Hierarchical structure
cell_indx = [3 4];
indiv_indx = [7];
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', indiv_indx, 'cell_indx',cell_indx, 'n_indiv', length(PI.u));
if ~isempty(PI.H.IndividualParams(1).Index)
        indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),PI.H.IndividualParams,'UniformOutput',false)';
else
    indivSigmaNames = [];
end
if ~isempty(PI.H.CellParams(1).Index)
    cellSigmaNames=arrayfun(@(x)strjoin({'psi', x.name}, '_'),PI.H.CellParams,'UniformOutput',false)';
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
alpha = [repelem(0.1, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.1, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.001, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
beta = [repelem(0.1, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.1, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.001, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
sigma_prior= [ repelem(1,length(PI.H.PopulationParams), 1);...
    repelem(1, length([PI.H.CellParams(:).Index]),1);
     repelem(1, length([PI.H.IndividualParams(:).Index]),1);...
    alpha];
lb = ([2.2339   4.1269  0.0029  0.0814  0.7148  100.21  0.0043  1.2e-6  54.7    2348.2]');
ub = ([52.332   24.963  1.2283  2.2747  4.8938  377.62  8.3409  0.0086  702.51  8.9e5]');

lb = ([2    3   0.002   0.05    1  90   0.07    1e-6  35   1600]');
ub = ([50   20  4       4       6  300  10      0.02  552  8.6e5]');
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones', 'LB',lb,'UB', ub);
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',...
    {'uniform/normal/inverse gamma/inverse gamma'}));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),exp(finalValues(end-length(observables)+1:end)),PI.normIndx);
paramNames = getParamNames(PI,sim, observables);
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc

PI.AIC = 2*length(PI.par)-2*obj_fun(finalValues)*(-1);
%% Save results
save('PI_CIM_ICB2_2.mat', 'PI')

load(strjoin({cd 'PI_CIM_ICB_1.mat'},'/'))

load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'