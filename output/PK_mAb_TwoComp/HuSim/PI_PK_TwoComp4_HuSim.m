%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_TwoComp/HuSim')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/TwoComp_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'TwoComp_CE';
variants = get(model, 'variants');

%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor'; 'CL'; 'Q23'; 'kint';'PDL1_Tumor'; 'PDL1_Blood'; 'ID'};
% Define outputs
observables={ 'ID_g_Blood_free'};
dataset_file_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Heery_2017.xlsx'...
    '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/Herbst_2014.xlsx'
    };

[PI,PI.u]=getDataSets(dataset_file_ext, 'subsetVariables', {'Blood_antiPDL1_free__ID_g_'}, 'species', 'human');

PI.u = PI.u(:,1);
PI.variableUnits={'%ID/g'};
PI.observablesPlot = { 'Blood_antiPDL1_Perc'};
dose = {'Blood.antiPDL1'};
sim=createSimFunction(model,parameters,observables, dose,variants(1),...
    'UseParallel', false);
PI.normIndx = [];
PI.model = 'PK-Two Compartment Model (Human)';
% Get initial values
PI.x_0 =[PI.data(:).dose]';

%% Hierarchical model simulation
PI.H = getHierarchicalStruct2(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [3],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');

% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-3  1e-3    1e-4    1e-4    1e-4    1e0     1e0];
ub = [1e1   1e1     1e2     1e2     1e1     1e6     1e6];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1)-1,beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
prior = {'U' 'U' 'U' 'U' 'U' 'U' 'U'};


% Log-ikelihood function
likelihood_fun=@(p)likelihood2(exp(p),sim,PI,'censoring',false);
% prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));
prior_fun=@(p)getPriorPDF(p,PI, prior);
prior_fun_MCMC=@(p)getPriorPDFMCMC(p,PI, prior);

paramNames = getParamNames(PI,sim, observables);
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc


%% Posterior samples
PI_Preclinical=load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_TwoComp/PI/TwoComp_4/PI_PK_TwoComp4_4.mat');
postSamples_PopParams = (PI_Preclinical.PI.postSamples(:,[1:7])) ;
omega = (PI_Preclinical.PI.postSamples(:,17));
sigma = PI_Preclinical.PI.postSamples(:,18);
N_pop = 600;
N_indiv = 5;
randIndx = randsample(size(postSamples_PopParams,1),N_pop, true);

Z_indiv = (randn(N_indiv*N_pop, 11)).*repmat(omega(randIndx), N_indiv,11);

scalingExp = [0.8 0.8 0.7 0.7 0 0 0];
scalingFactor = (77/.022).^(scalingExp);
Theta1 = [repmat(postSamples_PopParams(randIndx,:), N_indiv, 1) Z_indiv...
    repmat(omega(randIndx), N_indiv,1) repmat(sigma(randIndx), N_indiv,1) ];

Theta1(:,1:7) = Theta1(:,1:7) + log(scalingFactor);
%  plotBivariateMarginals_2(exp(Theta1(:,1:7)),...
%        'names',parameters([1:7]),'interpreter', 'tex')
% plotBivariateMarginals_2(Theta(:,8:18),...
%       'names',{PI.par([8:18]).name},'interpreter', 'tex')
Delta = repmat(postSamples_PopParams(randIndx,1:4), N_indiv, 1)- mean(repmat(postSamples_PopParams(randIndx,1:4), N_indiv, 1));
Theta2 = [repmat(log([PI.par(1:4).startValue]), N_pop*N_indiv, 1)+Delta Theta1(:,5:end)];
plotBivariateMarginals_2((Theta2(:,1:7)),...
       'names',parameters([1:7]),'interpreter', 'tex')
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi3(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);
tic
PI=getPosteriorPredictions(exp(Theta1),PI,simFun,PI.observablesPlot);
toc
tic
PI2 = getPosteriorPredictions(exp(Theta2), PI,simFun, PI.observablesPlot);
toc
PI=getCredibleIntervals(PI,PI.observablesPlot, (Theta1),PI.H);
plotPosteriorPredictions(PI,PI.observablesPlot,'output','indiv')

PI2 = getCredibleIntervals(PI2, PI2.observablesPlot, Theta2, PI.H);
PI2_subset = PI2;
PI2_subset.data = PI2.data(5:12);
plotPosteriorPredictions(PI2_subset, PI.observablesPlot, 'output', 'indiv')
PI2.postSamples = Theta2;
PI.postSamples = Theta1;
%% Save
save(strjoin({cd 'PI_PK_TwoComp4_HuSim_Theta1'}, '/'), 'PI')
save(strjoin({cd 'PI_PK_TwoComp4_HuSim_Theta2'}, '/'), 'PI2')
