%% PK general setup

% Search paths
clear all
warning off
%%
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output\PK_mAb_ThreeComp\PI\ThreeComp_7')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\ThreeComp_4.sbproj');

else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_7')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/ThreeComp_4.sbproj');

end

%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock',2.5)
variants = get(model,'variants');
sbioaccelerate(model,cs);
%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor';'CL_antiPDL1'; 'Q23'; ...
     'PDL1_Tumor'; 'kdeg_PDL1'; 'PDL1_Blood';'kint';'ID'};
% Define outputs
observables={'ID_Id_g_Blood' 'ID_g_Blood_free' 'ID_Id_g_Tumor' 'ID_g_Tumor_free'};
if ispc
    dataset_file_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Nedrow_2017_1.xlsx'...
        '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Nedrow_2017_1.xlsx'...
         '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Contreras_2016_2.xlsx'};
else
   dataset_file_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_1.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_2.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Contreras_2016_2.xlsx'...
    };

end


[PI,PI.u]=getDataSets(dataset_file_ext, 'subsetVariables', {'Blood_antiPDL1__ID_g_'...
    'Blood_antiPDL1_free__ID_g_' 'Tumor_antiPDL1__ID_g_' 'Tumor_antiPDL1_free__ID_g_'...
      });

PI.u = PI.u(:,1);
PI.variableUnits={'%ID/g'  '%ID/g' '%ID/g' '%ID/g' };
PI.observablesPlot = {'Blood Serum antiPDL1_{Total}' 'Blood Serum antiPDL1_{Free}'...
    'Tumor antiPDL1_{Total}' 'Tumor antiPDL1_{Free}' };
PI.observablesFields = {'Blood_Serum_Total', 'Blood_Serum_Free', 'Tumor_Total', 'Tumor_Free'};
dose = {'Blood.antiPDL1'};
sim=createSimFunction(model,parameters,observables, dose,variants(3),...
    'UseParallel', false);
PI.normIndx = [];
PI.model ='Three-Compartment model with TMDD';
% Get initial values
PI.x_0 =[PI.data(:).dose]';
clear dataset_file_ext dose MODEL 
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');

% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([.1 .1 .001], [.1 .1 0.001], [1 1 1], PI);
lb = [1e-2   1e-3   1e-4    1e-4   1e0 1e-3 1e0 1e-3];
ub = [1e1    2      1e1     1e1    1e6 1e0  1e6 1e1];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1)-1,beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');
PI = assignPrior(PI);

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);

paramNames = getParamNames(PI,sim, observables);
PI.paramNames = paramNames;
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
save('PI_PK_TwoComp4_7_TMDD_0.mat', 'PI')
load(strjoin({cd 'PI_PK_ThreeComp4_7_TMDD_0.mat'},'/'))

save(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_x2.mat'},''), 'x2')
save(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_p_x2.mat'},''), 'p_x2')
load(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_p_x1.mat'},''))
load(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_x1.mat'},''))