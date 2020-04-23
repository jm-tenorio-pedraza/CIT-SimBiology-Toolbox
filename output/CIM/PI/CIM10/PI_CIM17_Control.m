%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/CIM/PI/CIM10')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_4.sbproj');
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1; 0.1},...
    'variant', {variants(1); variants(2); variants(3)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-10);
set(cs, 'MaximumWallClock', 0.25)
sbioaccelerate(model, cs)
%% Parameter setup
parameters = {'kin_CD8';'K_IFNg';'KDE_MDSC';'S_L'; 'K_MDSC'; ...
    'kin_MDSC';'kin_TIC';'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; };
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
PI.normIndx = 6:8;
PI.model = 'CIM Control';
PI.observablesPlot={'TV' 'CD8' 'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};

% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1','doseUnits', 'mole');

%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [8 9 10] , 'cell_indx',[6:7], 'n_indiv', length(PI.u));
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues([1 1 .001], [1, 1 0.001], [1 1 1], PI);

lb=([1e-3   1e-3    1e-4 1e-7 1e-3   1e-4    1e-3    1e-2   1e-2     1e-3])';
ub=([1e2     1e2     1e1  1e1 1e1    1e2     1e2     1e1    1e2  1e3])';
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior, 'ref', 'ones','LB', lb, 'UB', ub);

prior = {'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U'};

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
save('PI_CIM17_Control_14.mat', 'PI')
load(strjoin({cd 'PI_CIM17_Control_14.mat'},'/'),'PI')

load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

