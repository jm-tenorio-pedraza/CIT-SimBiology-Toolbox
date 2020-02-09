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
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1; 1},...
    'variant', {variants(1); variants(2); variants(3)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-10);
set(cs, 'MaximumWallClock', 0.25)
%% Parameter setup
parameters = {'kin_CD8';'KDE_MDSC';'K_CTLA4'; ...
    'kpro_Tumor'; 'kill_CD8'; 'kin_Treg' ; ...
     'K_PDL1'; 'kill_Treg'; 'K_CD8'};
parameters = [parameters; 'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC2_Control'...
    'MOC1_antiPDL1' 'MOC2_antiPDL1' 'MOC1_antiCTLA4' 'MOC2_antiCTLA4'...
    'MOC1_antiCTLA4_antiPDL1' 'MOC2_antiCTLA4_antiPDL1'};
observables={'TV' };
stateVar={'Tumor'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    true,'output', 'mean','maxIIV', false);

PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Percentage [%]' ...
    'Relative units []' 'Relative units []'};
PI.normIndx =[];
PI.model = 'CIM ICB';
PI.observablesPlot={'TV' };
% simulationNames = {'MOC1 (Morisada 2017)'; 'MOC1 (Morisada 2017)'; 'MOC1 (Clavijo, 2017)'; 'MOC2 (Clavijo, 2017)'};
% [PI.data(:).Name] = simulationNames{:,:};
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1','doseUnits', 'mole');

%% Optimization setup
% Hierarchical structure
cell_indx = [3 4];
indiv_indx = [1];

PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', indiv_indx, 'cell_indx',cell_indx, 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
% Define lower and upper bounds for population params
lb=([3     9e-3    4        0.10    0.01	1.9      1   1e-6    1e-6])';
ub=([250   1.58    116.3    0.86    0.05	116      1e5 1e6     1e6])';
% Generate parameter estimation structure
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones','LB', lb, 'UB', ub);
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));
% Generate parameter names for plots
paramNames = getParamNames(PI,sim, observables);
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc
%% Parameter selection of inter-cell line varying params

w = arrayfun(@(x) (finalValues(x.Index)), PI.H.CellParams,'UniformOutput', false);
w = (std(cell2mat(w'),[], 2));
[w,cell_indx] =sort(w,'descend');
table({PI.H.CellParams(cell_indx).name}', w)


%% Save results
save('PI_CIM3_red2.mat', 'PI')
load(strjoin({cd 'PI_CIM_Control_3_full.mat'},'/'),'PI')

load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'