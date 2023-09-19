%% General setup
%% Search paths
clear all
warning off
opengl('save', 'software')
%%
if ispc
    wd='\Users\jmten\OneDrive\Dokumente\GitHub\';
    data_ext = '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_CIM10_MOC1_InVitro.mat';
else
    wd='/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/';
    data_ext = '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_CIM10_MOC1_InVitro.mat';
end
    addpath(genpath(strjoin({wd 'CIT-SimBiology-Toolbox'}, '')))
    cd(strjoin({wd 'CIT-SimBiology-Toolbox/output'},''))
    out=sbioloadproject(strjoin({wd 'sbio-projects/CIM_10_1_PKPD.sbproj'},''));
 
%% Load project 
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {0.01; 0.1; 0.1},...
    'variant', {variants(1); variants(3); variants(4)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', .5)
set(cs, 'Time','hour')
set(cs, 'StopTime',72)
    %% Rule setup
% Initial tumor volume
model.Rules(19,1).Active=false;
% Initial tumor cell number
model.Rules(1,1).Active=true;
% Initial PDL1_Tumor
model.Rules(16,1).Active=true;
% Initial CTLA4_CD8
model.Rules(17,1).Active=true;
% Initial CTLA4_Treg
model.Rules(20,1).Active=true;
% Initial PDL1_Immune
model.Rules(30,1).Active=true;
% Carrying capacity
model.Rules(47,1).Active=false;
% Observers
model.Rules(2,1).Active=false;
model.Rules(3,1).Active=false;
model.Rules(4,1).Active=false;
model.Rules(5,1).Active=false;
model.Rules(6,1).Active=false;
model.Rules(7,1).Active=false;
model.Rules(8,1).Active=false;
model.Rules(10,1).Active=false;

%% Reaction setup
BloodoutSet = sbioselect(model.Reactions,'Where', 'Reaction', 'regexp', 'Blood');
PeripheralSet = sbioselect(model.Reactions,'Where', 'Reaction', 'regexp', 'Peripheral');
LNSet = sbioselect(model.Reactions,'Where', 'Reaction', 'regexp', 'LymphNode');

CTLA4Set=sbioselect(model.Reactions, 'Where', 'Reaction', 'regexp', 'CTLA4');
CTLA4RuleSet=sbioselect(model.Rules, 'Where', 'Rule', 'regexp', {'kde_Treg =' 'CTLA4_RO ='});
set(BloodoutSet, 'Active', false);
set(PeripheralSet, 'Active', false);
set(LNSet, 'Active', false);
% set(CTLA4Set, 'Active', false);
set(CTLA4RuleSet, 'Active', false);

%% Parameter setup
parameters = {'kpro_Tumor';'kill_CD8'; ...
    'KDE_MDSC';'IARkill_w_PDL1';'IARkill_w_MDSC'; 'IARpro_w_PDL1';...
    'IARpro_w_MDSC';'K_kill';'T_max'};
parameters = [parameters;'T_0';'CD8_N_0'; 'GMDSC_0'; 'antiPDL1_InVitro'];
% Define outputs% Define outputs
observables={'Impendance'};
stateVar={'Impendance'};
doses = {'Blood.Dose_antiPDL1'};
%% Assemble PI struct
PI=load(data_ext);
PI=PI.SI;
PI.data=PI.data';
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))), {PI.data(:).dataValue}));
PI.tspan=unique([PI.data(:).dataTime]);

PI.variableUnits={'Relative units []' };
PI.observablesFields = {'Impendance'};
PI.normIndx = [];
PI.model = 'CIM10_PKPD_MOC1_InVitro';
PI.observablesPlot={'Impendance change'};
plotData(PI, PI.observablesPlot, 'responseGrouping', false, 'kineticGrouping', false)
%% Set T_0 as a variable parameters
% PI.x_0=PI.x_0(:,2:end);
%% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'InVitro_MOC1','doseUnits', 'mole');
% clearvars   PI1 PI2 variants data_ext data_ext1 doses ans Cell_Field groups_subset groupsResp ImmuneResp index
% close all
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct3(parameters(1:end-4),PI,'n_sigma', length(observables),...
    'rand_indx', [] , 'cell_indx',[1 2 ], 'resp_indx', [],'n_indiv', length(PI.u),...
    'CellField','Group');
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues2([.1 .1 .1 .1], ...
    [.6 .6 .6 .6], [1 1 1 1], PI);
% lb=([1e-2   1e-2    1e-2    1e-2    .01     1e-3    1e-3    1e-4  1e-4       1e-2   1e-3     1e-3    1e-2])';
% ub=([1e4    1e2     1e2     1e5     1       100     1e4     1e6   1e6        1e5    1e3      1e3     1e4])';
paramValues = sim.Parameters.Value(1:end-4);
lb = paramValues.*[1e-1     1e-2        1e-2    1e-2       1e-2      1e-2   1e-2  1e-3 .01]';
ub = paramValues.*[1e1      1e2         1e2     500        500       400    400   1e6  1]';

PI.par = getParamStruct3(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'startSigma',...
    ones(length(PI.H.SigmaParams), 1)*.6, 'ref', 'ones','LB', lb, 'UB', ub);

PI = assignPrior(PI,'sigmaDist', 'JP');
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false,...
    'logTransform',true,'errorModel','additive','constantVar',.01,'indivData',false);
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);


paramNames = getParamNames(PI,sim, observables);
PI.paramNames  =paramNames;
%% Objective function
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Obj function
obj_fun=@(x)(postLogLikelihood(x, prior_fun, likelihood_fun)*(-1));
tic
obj_fun((finalValues))
toc

