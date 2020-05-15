%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/HuSim')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/ThreeComp_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-15);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-12);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'ThreeComp (with TMDD)';
variants = get(model, 'variants');

%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor'; 'Peripheral'; 'CL_antiPDL1'; 'Q23'; 'Q12';'PDL1_Tumor'; 'kdeg_PDL1'; 'ID'};
parameters = [parameters; 'KD_antiPDL1'; 'koff_antiPDL1'];
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
PI.model = 'PK-Three Compartment Model (Human)';
% Get initial values

%% Hierarchical model simulation
PI.H = getHierarchicalStruct(parameters(1:end-3),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[2 4], 'n_indiv', length(PI.u),'CellField', 'Name');
KD_antiPDL1 =   [1.75 0.0467];
koff_antiPDL1 = [1.56e-4 0.753e-4];

PI.x_0 = [[PI.data(:).dose]' PI.H.CellIndx*KD_antiPDL1' PI.H.CellIndx*koff_antiPDL1'];
% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-3  1e-3    1e-4    1e-4    1e-4    1e0     1e0 1e0];
ub = [1e1   1e1     1e2     1e2     1e1     1e6     1e6 1e6];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1)-1,beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
prior = {'U' 'U' 'U' 'U' 'U' 'U' 'U' 'U'};


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
PI_Preclinical=load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_4/PI_PK_ThreeComp4_4_TMDD_11.mat');
postSamples_PopParams = (PI_Preclinical.PI.postSamples(:,[1:8])) ;
omega = (PI_Preclinical.PI.postSamples(:,[15 16]));
sigma = PI_Preclinical.PI.postSamples(:,17);
N_pop = 500;
N_indiv = 2;
randIndx = randsample(size(postSamples_PopParams,1),N_pop, true);

Z_indiv = (randn(N_indiv*N_pop, 4)).*repmat(omega(randIndx,:), N_indiv,2);
% Z = Z_indiv*PI.H.CellIndx';
scalingExp1 = [0.9 0.9 0.9 0.9 0.9 0.9 0 0];
scalingExp2 = [1 1 1 0.9 0.9 0.9 0 0];
scalingExp3 = [0.8 0.8 0.8 0.8 0.8 0.8 0 0];
scalingExp4 = [0.7 0.7 0.7 0.7 0.7 0.8 0 0];

scalingFactor = (77/.022).^(scalingExp1);
scalingFactor2 = (77/.022).^(scalingExp2);
scalingFactor3 = (77/.022).^(scalingExp3);
scalingFactor4 = (77/.022).^(scalingExp4);

Theta = [repmat(postSamples_PopParams(randIndx,:), N_indiv, 1) Z_indiv...
    repmat(omega(randIndx,:), N_indiv,1) repmat(sigma(randIndx), N_indiv,1) ];

Theta1 = [Theta(:,1:8) + log(scalingFactor) Theta(:,9:end)];
Theta2 = [Theta(:,1:8) + log(scalingFactor2) Theta(:,9:end)];
Theta3 = [Theta(:,1:8) + log(scalingFactor3) Theta(:,9:end)];
Theta4 = [Theta(:,1:8) + log(scalingFactor4) Theta(:,9:end)];

% Alternative parametrization
Delta = repmat(postSamples_PopParams(randIndx,1:6), N_indiv, 1)- mean(repmat(postSamples_PopParams(randIndx,1:6), N_indiv, 1));
Theta5 = [repmat(log([PI.par(1:6).startValue]), N_pop*N_indiv, 1)+Delta Theta1(:,7:end)];

table({PI.par(1:8).name}', exp(mean(Theta1(:,1:8)))',exp(mean(Theta2(:,1:8)))',...
    exp(mean(Theta3(:,1:8)))',exp(mean(Theta5(:,1:8)))')
%% Plot bivariate marginals
plotBivariateMarginals_2(exp(Theta2(:,1:8)),...
       'names',parameters([1:8]),'interpreter', 'tex')
% plotBivariateMarginals_2(Theta(:,8:18),...
%       'names',{PI.par([8:18]).name},'interpreter', 'tex')
%plotBivariateMarginals_2((Theta2(:,1:7)),...
      % 'names',parameters([1:7]),'interpreter', 'tex')
%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H);

PI1 = getPosteriorPredictions(exp(Theta1),PI,simFun,PI.observablesPlot);
PI2 = getPosteriorPredictions(exp(Theta2), PI,simFun, PI.observablesPlot);
PI3 = getPosteriorPredictions(exp(Theta3), PI,simFun, PI.observablesPlot);
PI4 = getPosteriorPredictions(exp(Theta4), PI,simFun, PI.observablesPlot);

PI1=getCredibleIntervals(PI1,PI1.observablesPlot, (Theta1),PI1.H);
PI2=getCredibleIntervals(PI2,PI2.observablesPlot, (Theta2),PI2.H);
PI3=getCredibleIntervals(PI3,PI3.observablesPlot, (Theta3),PI3.H);
PI4=getCredibleIntervals(PI4,PI4.observablesPlot, (Theta3),PI4.H);

plotPosteriorPredictions(PI1,PI1.observablesPlot,'output','indiv','central', 'mean')
plotPosteriorPredictions(PI2,PI2.observablesPlot,'output','indiv', 'central', 'mean')
plotPosteriorPredictions(PI3,PI3.observablesPlot,'output','indiv', 'central', 'mean')
plotPosteriorPredictions(PI4,PI4.observablesPlot,'output','indiv', 'central', 'mean')

PI4.postSamples = Theta4;
PI1.postSamples = Theta1;
PI2.postSamples = Theta2;
PI3.postSamples = Theta3;

%% Save
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_1_2'}, '/'), 'PI1')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_2_2'}, '/'), 'PI2')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_3_2'}, '/'), 'PI3')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_4_2'}, '/'), 'PI4')

load(strjoin({cd 'PI_PK_ThreeComp4_HuSim_3'}, '/'), 'PI3')

