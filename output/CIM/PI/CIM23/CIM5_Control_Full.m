%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM23')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_5.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1; 0.1},...
    'variant', {variants(1); variants(2); variants(3)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-8);
set(cs, 'MaximumWallClock', 2.5)
sbioaccelerate(model, cs)
%% Parameter setup
parameters = {'kin_CD8'; 'kin_Treg';'K_IFNg';'KDE_MDSC';'K_MDSC'; 'kin_DC';'S_L'; ...
    'S_R'; 'kill_Treg';'kin_MDSC';'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; 
    'K_DC'; 'ks_PDL1_Tumor'; 'ks_PDL1_Immune'; 'K_CTLA4'; 'K_PDL1';
    'kdeg_PDL1'; 'kdeg_CTLA4'; 'ks_CTLA4_CD8'; 'K_dif'; };
parameters = [parameters; 'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC2_Control',...
    'MOC2_Control_Mean' 'MC38_Control'};
observables={'TV'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD8_E' 'PDL1_T' 'PDL1_I'};
stateVar={'Tumor'  'CD8' 'Treg' 'DC'...
    'GMDSC' 'CD107a_Rel' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI1=getPIData3('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,'output', 'mean','responseGrouping', true, 'kineticGrouping', true);
PI2=getPIData3('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada.mat',...
    stateVar,groups_subset,'output', 'mean','responseGrouping', true, 'kineticGrouping', true);

PI.data = [PI1.data; PI2.data];
PI.n_data = PI1.n_data+PI2.n_data;
PI.tspan = unique([PI1.tspan; PI2.tspan]);
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Relative units []' ...
    'Relative units []' 'Relative units []'};
PI.observablesFields = {'TV'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD8_E' 'PDL1_T' 'PDL1_I'};
PI.normIndx = 6:8;
PI.model = 'CIM Control';
PI.observablesPlot={'TV' 'CD8' 'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};
plotData(PI, PI.observablesPlot, 'responseGrouping', true, 'kineticGrouping', true)
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);
% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1','doseUnits', 'mole');
clearvars doses groups_subset PI1 PI2
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [ ] , 'cell_indx',[ ], 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([.4 .4 .05], [.1 .1 .1], [1 1 1], PI);
lb=([1e-3   1e-3    1e-3    1e-4 1e-3   1e-3   1e-3 1e-3    1e-6    1e-4    1e-3   1e-2     1e-3    1e-3    1e-3    1e-3    1e0     1e0     1e-3    1e-3    1e-3    1e0     ])';
ub=([1e2    1e2     1e2     1e1 1e1    1e2     1e3  1e3     1e3     1e3     1e0    1e3      1e3     1e1     1e5     1e5     1e3     1e6     1e1     1e1     1e4     1e1     ])';

PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones','LB', lb, 'UB', ub);

PI = assignPrior(PI);
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);


paramNames = getParamNames(PI,sim, observables);
PI.paramNames= paramNames;
clearvars beta lb ub sigma_prior SigmaNames 

%% Objective function
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc

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
save('PI_CIM5_Control_Full_0.mat', 'PI')
load(strjoin({cd 'PI_CIM5_Control_Full_0.mat'},'/'),'PI')
N_i='3';
save(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_x_' N_i '.mat'},''), strjoin({'x' N_i},''))
save(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_p_x_' N_i '.mat'},''), strjoin({'p_x' N_i},''))
save(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_J' N_i '.mat'},''), strjoin({'J' N_i},''))
save(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_n_id' N_i '.mat'},''), strjoin({'n_id' N_i},''))
save(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_pCR' N_i '.mat'},''), strjoin({'pCR' N_i},''))
save(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_stepSize' N_i '.mat'},''), strjoin({'stepSize' N_i},''))

load(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_J' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_n_id' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_pCR' N_i '.mat'},''))
load(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_stepSize' N_i '.mat'},''))

for i=1:3
    load(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_x_' num2str(i) '.mat'},''))
    load(strjoin({cd '/PI_CIM22_Control_14_DREAM_MCMC_p_x_' num2str(i) '.mat'},''))
    

end

load(strjoin({cd 'CIM22_Control_14_DREAM_MCMC_x.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

