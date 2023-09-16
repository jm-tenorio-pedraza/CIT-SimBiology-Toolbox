%% PK general setup

% Search paths
clear all
warning off
%%
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output\PK_mAb_ThreeComp\PI\ThreeComp_10')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_10_1_PK_2.sbproj');

else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_10')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_10_1_PK_2.sbproj');

end

%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock',0.25)
variants = get(model,'variants');
sbioaccelerate(model,cs);
set(cs, 'Time','hour')
set(cs, 'StopTime',150)

%% Setting up parameters, data and simulations

parameters = {'Blood'; 'T_0';'Peripheral';'CL_antiPDL1';'Q12';...
     'PDL1_Tumor_ss'; 'kdeg_PDL1';'kint'; 'f_L'; 'V_vs'; 'vref_is';  'ID_antiPDL1'};

 % Define outputs
observables={'ID_g_Blood_antiPDL1_Id' 'ID_ml_antiPDL1_Blood_free'...
    'ID_g_Tumor_antiPDL1_Id' 'ID_ml_antiPDL1_Tumor_free'  };
if ispc
    dataset_file_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Nedrow_2017_1.xlsx'...
        '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Nedrow_2017_2.xlsx'...
         '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Contreras_2016_2.xlsx'...
         '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Pang_2018.xlsx'};
else
   dataset_file_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_1.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_2.xlsx'...
        '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Contreras_2016_2.xlsx'...
        '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Pang_2018.xlsx'};
end
%% Load data

[PI,PI.u]=getDataSets(dataset_file_ext, 'subsetVariables', {'Blood_antiPDL1__ID_g_'...
    'Blood_antiPDL1_free__ID_g_'  ...
    'Tumor_antiPDL1__ID_g_' 'Tumor_antiPDL1_free__ID_g_'...
      });

PI.u = PI.u(:,1);
PI.variableUnits={'%ID/g' '%ID/g' '%ID/g' '%ID/g'};
PI.observablesPlot = {'Blood Serum antiPDL1_{Total}' 'Blood Serum antiPDL1_{Free}'...
    'Tumor antiPDL1_{Total}' 'Tumor antiPDL1_{Free}' };
PI.observablesFields = {'Blood_Serum_Total', 'Blood_Serum_Free',...
    'Tumor_Total', 'Tumor_Free'};
dose = {'Blood.antiPDL1'};
PI.normIndx = [];
PI.model ='Four-Compartment model with TMDD';
% Get initial values
PI.x_0 =cell2mat({PI.data(:).dose}');
PI.x_0 = PI.x_0(:,1);

% Add response vector
respV = repelem({'Progressor'}, length(PI.data),1);
[PI.data(1:end).Response] = respV{:,:};
clear dataset_file_ext  MODEL 
%% 
sim=createSimFunction(model,parameters,observables, dose,variants(2),...
    'UseParallel', false);

%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct3(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[1 2 3],'resp_indx', [],'n_indiv', length(PI.u),'CellField', 'Name');

% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([.1 .1 .001 ], [.1 .1 0.001 ], [1 1 1], PI);
lb = [1e-1   1      1e-1    1e-2    1e-2    1e-1    1e-2     1e-3   1e-3    1e-3 1e-3];
ub = [1e1    1e2    2e1     1e2     1e2     1e3     1e2      1e2    0.5     0.2  0.9];
PI.par = getParamStruct3(sim,PI.H,size(PI.data,1)-1,beta,...
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
%%

%% Save results
save('PI_PK_ThreeComp10_1.mat', 'PI')
%%
load(strjoin({cd 'PI_PK_ThreeComp10_1.mat'},'/'))
%%
%% save MCMC results
N_i='7';
model_name = '/PI_PK_ThreeComp4_9_TMDD_8_';
save(strjoin({cd model_name 'x_' N_i '.mat'},''), strjoin({'x' N_i},''))
save(strjoin({cd model_name 'p_x_' N_i '.mat'},''), strjoin({'p_x' N_i},''))
save(strjoin({cd model_name 'J' N_i '.mat'},''), strjoin({'J' N_i},''))
save(strjoin({cd model_name 'n_id' N_i '.mat'},''), strjoin({'n_id' N_i},''))
save(strjoin({cd model_name 'stepSize' N_i '.mat'},''), strjoin({'stepSize' N_i},''))
%% load MCMC results
load(strjoin({cd model_name 'x_' N_i '.mat'},''))
load(strjoin({cd model_name 'p_x_' N_i '.mat'},''))
load(strjoin({cd model_name 'J' N_i '.mat'},''))
load(strjoin({cd model_name 'n_id' N_i '.mat'},''))
load(strjoin({cd model_name 'stepSize' N_i '.mat'},''))


for i=1:6
    load(strjoin({cd model_name 'x_' num2str(i) '.mat'},''))
    load(strjoin({cd model_name 'p_x_' num2str(i) '.mat'},''))
end
