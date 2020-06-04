%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/HuSim')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio projects/ThreeComp_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-09);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'ThreeComp (with TMDD)';
variants = get(model, 'variants');
sbioaccelerate(model, cs)
%% Setting up parameters, data and simulations
parameters = {'Tumor'; 'Q23';'ID'};
parameters = [parameters;'Blood'; 'Peripheral'; 'CL_antiPDL1';'Q12';...
    'KD_antiPDL1'; 'koff_antiPDL1'];
% Define outputs
observables={ 'ID_g_Blood_free'};
dataset_file_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Heery_2017.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Herbst_2014.xlsx'
    };

[PI,PI.u]=getDataSets(dataset_file_ext, 'subsetVariables', {'Blood_antiPDL1_free__ID_g_'}, 'species', 'human');

PI.u = PI.u(:,1);
PI.variableUnits={'%ID/g'};
PI.observablesPlot = { 'Blood antiPDL1 Perc'};
PI.observablesFields = {'Blood_antiPDL1_Perc'};
dose = {'Blood.antiPDL1'};
sim=createSimFunction(model,parameters,observables, dose,variants(1),...
    'UseParallel', false);
PI.normIndx = [];
PI.model = 'PK-Three Compartment Model (Human)';
%% Prior parameters
                % Atezo     Avelu
KD_antiPDL1 =   [1.75       0.0467];
koff_antiPDL1 = [1.56e-4    0.753e-4];
CL_antiPDL1 =   [8.33       30.8];
Q21 =           [22.75      31.3];
Blood =         [3280       3420];
Peripheral =    [3630       918];

priorPK = [Blood; Peripheral; CL_antiPDL1; Q21; KD_antiPDL1; koff_antiPDL1];
%% Hierarchical model simulation
PI.H = getHierarchicalStruct(parameters(1:end-7),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');

PI.x_0 = [[PI.data(:).dose]' PI.H.CellIndx*priorPK'];
% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e0  1e0  ];
ub = [1e2  1e2   ];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');

PI = assignPrior(PI);

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
% prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);

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


