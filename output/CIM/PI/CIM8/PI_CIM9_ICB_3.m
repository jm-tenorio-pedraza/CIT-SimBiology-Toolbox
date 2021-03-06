%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM8')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_4.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1;},...
    'variant', {variants(5); variants(6)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-10);
set(cs, 'MaximumWallClock', 0.25)
sbioaccelerate(model, cs)
%% Parameter setup
parameters = {'kpro_Tumor';'kill_CD8';'kpro_Tumor_Linear'; 'K_MDSC';'S_L';...
    'S_R'; 'K_CTLA4';'K_PDL1'; 'kill_Treg'};
parameters = [parameters; 'kin_MDSC';'f3';  'T_0'];

% Define outputs% Define outputs
groups_subset = {'MOC1_antiCTLA4', 'MOC1_antiPDL1', 'MOC1_antiCTLA4_antiPDL1' 'MOC2_antiCTLA4',...
    'MOC2_antiPDL1' 'MOC2_antiCTLA4_antiPDL1'  'MC38_antiPD1'};
observables={'TV'};
stateVar={'Tumor'  };
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table

PI=getPIData3('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,'output', 'mean','responseGrouping', true, 'kineticGrouping', true);
% PI2=getPIData3('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada.mat',...
%     stateVar,groups_subset,'output', 'mean','responseGrouping', true, 'kineticGrouping', true);

% PI.data = [PI1.data; PI2.data];
% PI.n_data = PI1.n_data+PI2.n_data;
% PI.tspan = unique([PI1.tspan; PI2.tspan]);
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Relative units []' ...
    'Relative units []' 'Relative units []'};
PI.normIndx = [];
PI.model = 'CIM9_3 ICB';
PI.observablesPlot={'TV'};

% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct, 'parameters', {'kin_MDSC' 'f3'});

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1_optimized','doseUnits', 'mole');

%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-3),PI,'n_sigma', length(observables),...
    'rand_indx', [5] , 'cell_indx',[], 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([1 1 .001], [1, 1 0.001], [1 1 1], PI);

lb=([.01    0.5    0.1    .05  1e-4     1e-4   1e0   1e0    1e-6])';
ub=([.16    7      1.7    .12     1      1     1e4   1e7    1e1])';
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones','LB', lb, 'UB', ub);

prior = {'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U'};

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
% prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));
prior_fun=@(p)getPriorPDF(p,PI, prior);
prior_fun_MCMC=@(p)getPriorPDFMCMC(p,PI, prior);

paramNames = getParamNames(PI,sim, observables);
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
save('PI_CIM9_ICB_3_0.mat', 'PI')
load(strjoin({cd 'PI_CIM9_Control_3_3.mat'},'/'),'PI')

load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

