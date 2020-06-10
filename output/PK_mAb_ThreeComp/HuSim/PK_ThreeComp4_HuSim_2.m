%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/HuSim')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/ThreeComp_4.sbproj');
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

parameters = {'Blood'; 'Tumor'; 'Peripheral'; 'CL_antiPDL1'; 'Q23'; 'Q12';...
    'PDL1_Tumor'; 'kdeg_PDL1'; 'ID'};
parameters = [parameters; 'KD_antiPDL1'; 'koff_antiPDL1'];
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
PI.H = getHierarchicalStruct(parameters(1:end-3),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[1 3 4 6], 'n_indiv', length(PI.u),'CellField', 'Name');

PI.x_0 = [[PI.data(:).dose]' PI.H.CellIndx*priorPK(5:6,:)'];
% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-1  1e-1    1e-1    1e-2    1e-2    1e-2    1e0 1e-1];
ub = [1e4   1e4     1e4     1e3     1e3     1e2     1e6 1e1];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
PI = assignPrior(PI);

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
% prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type',{'uniform/normal/inverse gamma/inverse gamma'}));
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);

paramNames = getParamNames(PI,sim, observables);
%% Objective function

% Obj function
obj_fun=@(x)(likelihood_fun(x)*(-1)+prior_fun(x)*(-1));
tic
obj_fun((finalValues))
toc


%% Posterior samples
PI_Preclinical=load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_4/PI_PK_ThreeComp4_4_TMDD_11.mat');
Meta(1).Struct = PI_Preclinical;
N_pop = 500;
N_indiv = 5;
N_cell = 2;
theta = getTheta(Meta, parameters(1:8), N_pop, N_indiv,N_cell);
plotCorrMat(theta, parameters(1:8))
plotBivariateMarginals_2((theta(:,1:8)),...
       'names',parameters([1:8]),'interpreter', 'tex')
%%

randIndx = randsample(size(PI_Preclinical.PI.postSamples,1),N_pop, true);
sigma = repelem(PI_Preclinical.PI.postSamples(randIndx,17),N_indiv*N_cell,1);

scalingExp1 = [1    1       1       0.9     0.9     0.9     0   0]; 
scalingExp2 = [0.9  0.9     0.9     0.9     0.9     0.9     0   0];
scalingExp3 = [0.9  0.9     0.9     0.85    0.85    0.85    0   0];
scalingExp4 = [0.9  0.9     0.9     0.8     0.8     0.8     0   0];
scalingExp5 = [0.9  0.9     0.9     0.75    0.75    0.75    0   0];
scalingExp6 = [0.9  0.9     0.9     0.7     0.7     0.7     0   0];

scalingFactor1 = (70/.02).^(scalingExp1);
scalingFactor2 = (70/.02).^(scalingExp2);
scalingFactor3 = (70/.02).^(scalingExp3);
scalingFactor4 = (70/.02).^(scalingExp4);
scalingFactor5 = (70/.02).^(scalingExp5);
scalingFactor6 = (70/.02).^(scalingExp6);

Theta = [theta repmat(sigma(randIndx), N_cell*N_indiv,1)];

Theta1 = [Theta(:,1:8) + log(scalingFactor1)];
Theta2 = [Theta(:,1:8) + log(scalingFactor2)];
Theta3 = [Theta(:,1:8) + log(scalingFactor3)];
Theta4 = [Theta(:,1:8) + log(scalingFactor4)];
Theta5 = [Theta(:,1:8) + log(scalingFactor5)];
Theta6 = [Theta(:,1:8) + log(scalingFactor6)];

%% Alternative parametrization
Delta = theta(:,1:6) - mean(theta(:,1:6));
priorMu = mean(log(priorPK(1:4,:)),2);
Mu = [priorMu(1) log(PI.par(2).startValue)...
    priorMu(2:3)' log(PI.par(5).startValue)...
    priorMu(4)];
Eta = reshape((log(priorPK(1:4,:)) - priorMu)',1,[]);

Theta7 = [repmat(Mu, N_pop*N_indiv*N_cell, 1)+Delta, theta(:,7:8),...
    repmat(Eta, N_pop*N_indiv*N_cell,1), sigma];
plotBivariateMarginals_2((Theta7(:,1:8)),...
       'names',parameters([1:8]),'interpreter', 'tex')
%% Compare parameter mean estimates
table({PI.par(1:8).name}', exp(mean(Theta1(:,1:8)))',exp(mean(Theta2(:,1:8)))',...
    exp(mean(Theta3(:,1:8)))',exp(mean(Theta4(:,1:8)))',exp(mean(Theta5(:,1:8)))',exp(mean(Theta6(:,1:8)))', exp(mean(Theta7(:,1:8)))')
%% Posterior predictions
simTime = unique([PI.tspan' 1:4:PI.tspan(end)]);
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,simTime),x,...
    [],PI.normIndx, PI.H, 'simTime', simTime);
simFun2=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,simTime),x,...
    @(x)getPhi2(x, PI.H,size(PI.data,1),'initialValue',PI.x_0),...
    PI.normIndx, PI.H, 'simTime', simTime);

PI1 = getPosteriorPredictions(exp(Theta1(1:end,:)),PI,simFun,PI.observablesFields, 'simTime', simTime);
PI2 = getPosteriorPredictions(exp(Theta2), PI,simFun, PI.observablesFields, 'simTime', simTime);
PI3 = getPosteriorPredictions(exp(Theta3), PI,simFun, PI.observablesFields, 'simTime', simTime);
PI4 = getPosteriorPredictions(exp(Theta4), PI,simFun, PI.observablesFields, 'simTime', simTime);
PI5 = getPosteriorPredictions(exp(Theta5), PI,simFun, PI.observablesFields, 'simTime', simTime);
PI6 = getPosteriorPredictions(exp(Theta6), PI,simFun, PI.observablesFields, 'simTime', simTime);
PI7 = getPosteriorPredictions(exp(Theta7), PI,simFun2, PI.observablesFields, 'simTime', simTime);

PI1=getCredibleIntervals(PI1,PI1.observablesFields, Theta1,PI1.H, 'sigma', sigma);
PI2=getCredibleIntervals(PI2,PI2.observablesFields, Theta2,PI2.H);
PI3=getCredibleIntervals(PI3,PI3.observablesFields, Theta3,PI3.H);
PI4=getCredibleIntervals(PI4,PI4.observablesFields, Theta4,PI4.H);
PI5=getCredibleIntervals(PI5,PI5.observablesFields, Theta5,PI5.H);
PI6=getCredibleIntervals(PI6,PI6.observablesFields, Theta6,PI6.H);
PI7=getCredibleIntervals(PI7,PI7.observablesFields, Theta7,PI7.H);

plotPosteriorPredictions(PI1,1,'output','indiv','central', 'median','color', 'cell','simTime', simTime)
plotPosteriorPredictions(PI2,1,'output','indiv','central', 'median','color', 'cell')
plotPosteriorPredictions(PI3,1,'output','indiv','central', 'median','color', 'cell')
plotPosteriorPredictions(PI4,1,'output','indiv','central', 'median','color', 'dataset')
plotPosteriorPredictions(PI5,1,'output','indiv','central', 'median','color', 'dataset')
plotPosteriorPredictions(PI6,1,'output','indiv','central', 'median','color', 'dataset')
plotPosteriorPredictions(PI7,1,'output','indiv','central', 'median','color', 'cell')

PI4.postSamples = Theta4;
PI1.postSamples = Theta1;
PI2.postSamples = Theta2;
PI3.postSamples = Theta3;
PI5.postSamples = Theta5;
PI6.postSamples = Theta6;
PI7.postSamples = Theta7;

%% Calculate MSRE
table([mean(PI1.MSE); mean(PI2.MSE); mean(PI3.MSE); mean(PI4.MSE);mean(PI5.MSE);
    mean(PI6.MSE); mean(PI7.MSE)], 'RowNames',{'scaling1' 'scaling2' 'scaling3'...
    'scaling4' 'scaling5' 'scaling6' 'scaling7'}, 'VariableNames', {'Mean_Squared_Error'})

%% Save
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_1_2'}, '/'), 'PI1')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_2_2'}, '/'), 'PI2')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_3_2'}, '/'), 'PI3')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_4_2'}, '/'), 'PI4')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_5_2'}, '/'), 'PI5')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_6_2'}, '/'), 'PI6')
save(strjoin({cd 'PI_PK_ThreeComp4_HuSim_7_2'}, '/'), 'PI7')

load(strjoin({cd 'PI_PK_ThreeComp4_HuSim_3'}, '/'), 'PI3')

