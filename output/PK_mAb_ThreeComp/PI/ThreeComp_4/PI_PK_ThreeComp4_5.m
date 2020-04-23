%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_4')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/ThreeComp_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'TwoComp_CE';
variants = get(model,'variants');
%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor';'Peripheral';'CL'; 'Q23';'Q12';'PDL1_Tumor'; 'kdeg_PDL1'; 'ID'};
% Define outputs
observables={'ID_Id_g_Blood' 'ID_Id_g_Tumor'};

dataset_file_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_1.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_2.xlsx'...
    };

[PI,PI.u]=getDataSets(dataset_file_ext, 'subsetVariables', {'Blood_antiPDL1__ID_g_'...
     'Tumor_antiPDL1__ID_g_' ...
      });

% PI.data = [PI.data(8); PI.data(1:7); PI.data(9:10)];
% PI.u = PI.u([8 1:7 9:10],1);
PI.u = PI.u(:,1);
PI.variableUnits={'%ID/g'  '%ID/g' };
PI.observablesPlot = {'Serum_antiPDL1_Indium' ...
    'Tumor_antiPDL1_Indium' };

dose = {'Blood.antiPDL1'};
sim=createSimFunction(model,parameters,observables, dose,variants(5),...
    'UseParallel', false);
PI.normIndx = [];
PI.model = 'PK-Two Compartment Model';
% Get initial values
PI.x_0 =[PI.data(:).dose]';
clear dataset_file_ext dose MODEL 
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');

% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([.1 .1 .001], [.1, .1 0.001], [1 1 1], PI);
lb = [1e-1   1e-3   1e-2      1e-2    1e-3  1e-3    1e0     1e-3];
ub = [1e1    2      1e1       1e1     1e1   1e1     1e6     1e2];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1)-1,beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');

prior = {'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U'};


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

%% Save results
save('PI_PK_TwoComp4_5_TMDD_0.mat', 'PI')
load(strjoin({cd 'PI_PK_TwoComp4_3_TMDD_21.mat'},'/'))

save(strjoin({cd '/PI_PK_TwoComp4_3_TMDD_21_DREAM_MCMC_x.mat'},''), 'x')
save(strjoin({cd '/PI_PK_TwoComp4_3_TMDD_21_DREAM_MCMC_p_x.mat'},''), 'p_x')
load(strjoin({cd '/PI_PK_TwoComp4_2_TMDD_2_red_DREAM_MCMC_p_x.mat'},''))
load(strjoin({cd '/PI_PK_TwoComp4_2_TMDD_2_red_DREAM_MCMC_x.mat'},''))