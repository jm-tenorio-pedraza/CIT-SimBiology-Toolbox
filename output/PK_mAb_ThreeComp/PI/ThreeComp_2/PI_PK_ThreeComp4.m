%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_2')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/ThreeComp_2.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'TwoComp_CE';
%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor'; 'Peripheral'; 'CL'; 'Q12'; 'Q23'; 'kint'; 'ID'};
% Define outputs
observables={'ID_Id_g_Blood','ID_Id_g_Tumor',...
    'ID_g_Tumor_free', 'T2B' };

dataset_file_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_1.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Nedrow_2017_2.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Contreras_2016.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Deng_2016.xlsx'};

[PI,PI.u]=getDataSets(dataset_file_ext, 'subsetVariables', { 'Blood_antiPDL1__ID_g_',...
    'Tumor_antiPDL1__ID_g_', ...
    'Tumor_antiPDL1_free__ID_g_', 'Tumor_to_Blood__' });

% Eliminate the data points where ATA was observed in Deng et al (2016):
% deng= PI.data(11:13);
% dataTime_deng = arrayfun(@(x) x.dataTime(1:5), deng, 'UniformOutput', false);
% dataValue_deng = arrayfun(@(x) x.dataValue(1:5,:), deng, 'UniformOutput', false);
% sd_deng =arrayfun(@(x) x.SD(1:5,:), deng, 'UniformOutput', false);
% [PI.data(11:13).dataTime] = dataTime_deng{:,:};
% [PI.data(11:13).dataValue] = dataValue_deng{:,:};
% [PI.data(11:13).SD] = sd_deng{:,:};
% PI.n_data = PI.n_data - 3;
PI.u = PI.u(:,1);
PI.variableUnits={'%ID/g' '%ID/g' '%ID/g' '[]' '%ID/g' '%ID/g' '[]' };
PI.observablesPlot = {'ID_g_Blood' ...
    'ID_g_Tumor' 'ID_g_Tumor_free' 'T2B' };

dose = {'Blood.antiPDL1'};
sim=createSimFunction(model,parameters,observables, dose,[],...
    'UseParallel', false);
PI.normIndx = [];
PI.model = 'PK-Two Compartment Model';
% Get initial values
PI.x_0 =[PI.data(:).dose]';

%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [4],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');

% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-3  1e-3    1e-3    1e-3    1e-3    1e-4    1e-3];
ub = [1e3   1e3     1e3     1e3     1e3     1e3     1e3];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'LB', lb', 'UB', ub');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    (@(x)getCovariance(x,PI.H)),PI.normIndx,'log',true);
prior = {'U' 'U' 'U' 'U' 'U' 'U' 'U'};

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
% prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));
prior_fun=@(p)getPriorPDF(p,PI, prior);

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),...
    exp(finalValues(end-length(observables)+1:end)),PI.normIndx,'log', true);
paramNames = getParamNames(PI,sim, observables);
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc

%% Save results
save('PI_PK_ThreeComp5.mat', 'PI')
load(strjoin({cd 'PI_PK_ThreeComp16.mat'},'/'))

save(strjoin({cd '/PK_red_DREAM_MCMC_x.mat'},''), 'x')
save(strjoin({cd '/PK_red_DREAM_MCMC_p_x.mat'},''), 'p_x')
load(strjoin({cd '/PK_red_DREAM_MCMC_x.mat'},''))
load(strjoin({cd '/PK_red_DREAM_MCMC_p_x.mat'},''))