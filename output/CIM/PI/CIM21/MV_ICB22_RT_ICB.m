%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM21')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_4.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1;0.1},...
    'variant', {variants(5); variants(6);variants(7)});

cs=model.getconfigset;
set(cs,'SolverType', 'sundials');

set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 0.25)
sbioaccelerate(model, cs)
%% Parameter setup
parameters = {'kin_CD8'; 'K_IFNg'; 'KDE_MDSC'; 'K_MDSC';'kin_MDSC'; 'kin_TIC';...
    'kpro_Tumor';'kpro_Tumor_Linear';'kill_CD8'; 'S_L';'S_R'; 'K_CTLA4';...
    'K_PDL1'; 'kill_Treg'};
parameters = [parameters;  'T_0'; 'RT_8Gyx2_on'; 'RT_2Gyx10_on'];

% Define outputs% Define outputs
groups_subset = {'MOC1_8Gyx2' 'MOC1_2Gyx10' 'MOC1_antiPD1_8Gyx2' 'MOC1_antiPD1_2Gyx10'...
    'MC38_antiPD1' 'MC38_8Gyx2' 'MC38_2Gyx10' ...
    'MC38_antiPD1_8Gyx2' 'MC38_antiPD1_2Gyx2'};
observables={'TV'};
stateVar={'Tumor'  };
doses = {'Blood.Dose_antiPDL1'};
%% Obtain data, simulation function and dose table

PI=getPIData3('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada.mat',...
    stateVar,groups_subset,'output', 'mean','responseGrouping', true, 'kineticGrouping', true);
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'CIM21 ICB';
PI.observablesPlot={'TV'};
PI.observablesFields = {'TV'};
plotData(PI,PI.observablesPlot,'responseGrouping',true, 'kineticGrouping',true)
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);
PI.x_0 = [PI.x_0, [0 0 1 1 0 0 0 0 1 1 0 0 1 1 1 1]', [1 1 0 0 0 0 1 1 0 0 1 1 0 0 0 0]'];
% Get simulation function
[sim,~]=initializePI(model, parameters,observables, PI,doses,...
    'MOC1_optimized','doseUnits', 'mole');
PI.u = getDoseSchedule(PI, {'Blood.Dose_antiPD1'},...
    'startTime', [7 7],'freq', [5 5], 'numDoses', [3 3], 'dose', [0.0013 6.6667e-4],...
    'doseUnits', 'mole');
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-3),PI,'n_sigma', length(PI.observablesPlot),...
    'rand_indx', [ 7 8 9] , 'cell_indx',[5 6], 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([1 1 .001], [1, 1 0.001], [1 1 1], PI);

PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones');

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
clearvars ans beta doses groups_subset initialStruct lb  out sigma_prior SigmaNames stateVar ub variants
%% Posterior samples
PI_CIM = load("/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM21/PI_CIM21_Control_14.mat");
PI_ICB = load("/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM21/PI_CIM4_ICB_21_14.mat");

randIndx_CIM = randsample(1:size(PI_CIM.PI.postSamples,1), 1e3, true);
randIndx_ICB = randsample(1:size(PI_ICB.PI.postSamples,1),1e3, true);
varIndx_CIM = [1:6 10:11 13:14 43:44];
varIndx_ICB =   [1:8 repelem(17:18, 1, 5) 9 11 9 11 9 10 ...
    repelem(27:28,1, 5) 19 21 19 21 19 20  ...
    repelem(37:38, 1, 5) 29 31 29 31 29 30 39:42];
postSamples_cim = PI_CIM.PI.postSamples(:,varIndx_CIM);
postSamples_icb = PI_ICB.PI.postSamples(:, varIndx_ICB);
postSamples = [postSamples_cim(randIndx_CIM,1:6), postSamples_icb(randIndx_ICB, 1:8),...
    postSamples_cim(randIndx_CIM, 7:end-2), postSamples_icb(randIndx_ICB, 9:end-4)...
    postSamples_cim(randIndx_CIM, end-1:end), postSamples_icb(randIndx_ICB, end-3:end)];

%% Posterior predictions
simTime = unique([PI.tspan', 1:PI.tspan(end)]);
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,simTime),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H,...
    'output', 'simoutput', 'simTime', simTime);
tic
PI=getPosteriorPredictions2(exp(postSamples(1:end,:)),PI,simFun,PI.observablesFields,...
    'simTime', simTime);
toc
PI=getCredibleIntervals(PI,PI.observablesFields, exp(postSamples(1:end,:)),PI.H,...
    'logit_indx', [],'simTime', simTime);
