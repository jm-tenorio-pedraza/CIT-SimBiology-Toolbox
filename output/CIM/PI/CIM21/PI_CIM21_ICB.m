%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM21')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_4.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1;},...
    'variant', {variants(5); variants(6)});

cs=model.getconfigset;
set(cs.SolverType, 'ode15s');

set(cs.Sol, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 0.25)
sbioaccelerate(model, cs)
%% Parameter setup
parameters = {'kpro_Tumor';'kpro_Tumor_Linear';'kill_CD8'; 'S_L';...
    'S_R'; 'K_CTLA4';'K_PDL1'; 'kill_Treg'};
parameters = [parameters; 'kin_MDSC';'kin_TIC';  'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_antiCTLA4', 'MOC1_antiPDL1', 'MOC1_antiCTLA4_antiPDL1' 'MOC2_antiCTLA4',...
    'MOC2_antiPDL1' 'MOC2_antiCTLA4_antiPDL1'};
observables={'TV'};
stateVar={'Tumor'  };
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData3('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,'output', 'mean','responseGrouping', true, 'kineticGrouping', true);
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Relative units []' ...
    'Relative units []' 'Relative units []'};
PI.normIndx = [];
PI.model = 'CIM21 ICB';
PI.observablesPlot={'TV'};
plotData(PI,PI.observablesPlot,'responseGrouping',true, 'kineticGrouping',true)
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct, 'parameters', {'kin_MDSC' 'kin_TIC'});

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1_optimized','doseUnits', 'mole');

%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-3),PI,'n_sigma', length(observables),...
    'rand_indx', [1:3] , 'cell_indx',[], 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([1 1 .001], [1, 1 0.001], [1 1 1], PI);

lb=([.07   1.6  150    1e-7     1e-7   1e0   1e3    1e-4])';
ub=([.3    10    900    1        1      1e4   1e7     1e3])';
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones','LB', lb, 'UB', ub);

PI = assignPrior(PI);
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun_MCMC=@(p)getPriorPDFMCMC2(exp(p),PI);

paramNames = getParamNames(PI,sim, observables);
%% Objective function
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun_MCMC(x)*(-1));
tic
obj_fun((finalValues))
toc
clearvars ans beta doses groups_subset initialStruct lb observables out sigma_prior SigmaNames stateVar ub variants
%% Parameter selection of inter-cell line varying params

w = arrayfun(@(x) (finalValues(x.Index)), PI.H.CellParams,'UniformOutput', false);
w = (std(cell2mat(w'),[], 2));
[w,cell_indx] =sort(w,'descend');
cell_params = [ {PI.H.CellParams(:).name}'];

z = arrayfun(@(x) (finalValues(x.Index)), PI.H.IndividualParams,'UniformOutput', false);
z = (std(cell2mat(z'),[], 2));
[z,ind_indx] =sort(z,'descend');
ind_params = [{PI.H.IndividualParams(:).name}'];

table([cell_params(cell_indx); ind_params(ind_indx)], [w; z])
%% Save results
save('PI_CIM4_ICB_21_1.mat', 'PI')
load(strjoin({cd 'PI_CIM4_ICB_21_1.mat'},'/'),'PI')
N='12'
save(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_x_' N '.mat'},''), strjoin({'x' N},''))
save(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_p_x_' N '.mat'},''), strjoin({'p_x' N},''))
save(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_J' N '.mat'},''), strjoin({'J' N},''))
save(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_n_id' N '.mat'},''), strjoin({'n_id' N},''))
save(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_pCR' N '.mat'},''), strjoin({'pCR' N},''))
save(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_stepSize' N '.mat'},''), strjoin({'stepSize' N},''))

load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_x_6.mat'},''))
load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_p_x_6.mat'},''))
load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_J6.mat'},''))
load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_n_id6.mat'},''))
load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_pCR6.mat'},''))

load(strjoin({cd 'PI_CIM21_ICB_1_DREAM_MCMC_x_1.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))
for i=10:11
    load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_x_' num2str(i) '.mat'},''))
    load(strjoin({cd '/PI_CIM21_ICB_1_DREAM_MCMC_p_x_' num2str(i) '.mat'},''))

end
