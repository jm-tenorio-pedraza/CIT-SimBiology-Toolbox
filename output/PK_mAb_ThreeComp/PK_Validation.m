%% PK general setup

% Search paths
clear all
warning off
%%
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output\PK_mAb_ThreeComp\PI\ThreeComp_9')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_10.sbproj');

else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_9')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_10.sbproj');

end

%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock',2.5)
set(cs, 'TimeUnits','hour')
set(cs, 'Time','hour')
set(cs, 'StopTime',500)

variants = get(model,'variants');
sbioaccelerate(model,cs);
%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Peripheral';'CL_antiPDL1'; 'Q12';...
      'kdeg_PDL1'; 'PDL1_Blood';'PDL1_Peripheral';'kint';'ID_antiPDL1'; 'ID_antiCTLA4'};
observables={'ID_ml_antiPDL1_Blood_free' 'ID_ml_antiCTLA4_Blood_free'};
stateVar={ 'antiPDL1_ID_g_Blood_free' 'antiCTLA4_ID_g_Blood_free'};
doses = {'Blood.antiPDL1' 'Blood.antiCTLA4'};

 %%
% Define outputs
if ispc
    dataset_file_ext = '\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Schofield.mat';
else
   dataset_file_ext = '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Schofield.mat';

end

PI=load(dataset_file_ext);
PI=PI.PI;
dataValue = arrayfun(@(x)x.dataValue(:, [ 2 4]), PI.data,'UniformOutput', false);
sd = arrayfun(@(x)x.SD(:, [ 2 4]), PI.data,'UniformOutput', false);
[PI.data(1:end).dataValue]=dataValue{:,:};
[PI.data(1:end).SD] = sd{:,:};
PI.data = PI.data';
PI.n_data = sum(arrayfun(@(x)sum(sum(~isnan(x.dataValue))), PI.data));
PI.variableUnits={'%ID/g'  '%ID/g'};
PI.observablesPlot = {'Blood Serum anti-PD-L1_{Free}' 'Blood Serum anti-CTLA-4_{Free}' };
PI.observablesFields = {'Blood_Serum_antiPDL1_ID_ml','Blood_Serum_antiCTLA4_ID_ml'};

dose = {'Blood.antiPDL1' 'Blood.antiCTLA4'};
sim=createSimFunction(model,parameters,observables, dose,variants(8),...
    'UseParallel', false);
PI.normIndx = [];
PI.model ='Three-Compartment model with TMDD';
% Get initial values
PI.x_0 =cell2mat({PI.data(:).dose}');
clear dataset_file_ext  MODEL 
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-2),PI,'n_sigma', length(observables),...
    'rand_indx', [],'cell_indx',[], 'n_indiv', length(PI.u),'CellField', 'Name');

% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([.1 .1 .001], [.1 .1 0.001], [1 1 1], PI);
lb = [1e-2   1e-3    1e-4    1e-4     1e-2 1e0 1e0 1e-3];
ub = [1e1    1e1     1e1     1e1      1e1  1e6 1e6 1e1];
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

%% Load posterior samples
PI_fit = load(strjoin({cd 'PI_PK_ThreeComp4_9_TMDD_8.mat'},'/'));
PI_fit =PI_fit.PI;

[varIndx, refIndx]= ismember(parameters(1:8), {PI_fit.par(:).name}');
[sigmaIndx, sigmaRefIndx] = ismember({PI_fit.par(:).name}', 'sigma_ID_ml_antiPDL1_Blood_free');

% Simulate IEV
psiIndx = ismember({PI_fit.par(:).name}', 'psi_CL_antiPDL1');
psi_CL_antiPDL1 = PI_fit.postSamples(:,psiIndx);

z_CL_antiPDL1 = randn(size(psi_CL_antiPDL1)).*exp(psi_CL_antiPDL1);
etaIndx = (ismember(parameters(1:8), 'CL_antiPDL1'));

% Parameter set with population params and sigma params:
postSamples = [PI_fit.postSamples(:, refIndx) PI_fit.postSamples(:, sigmaIndx) PI_fit.postSamples(:, sigmaIndx)];
postSamples(:,etaIndx) = postSamples(:,etaIndx)+ z_CL_antiPDL1;

% Adjusting CL and Blood
bloodIndx = ismember(parameters(1:8), 'Blood');
CL_antiPDL1Indx = ismember(parameters(1:8), 'CL_antiPDL1');

bloodSample = postSamples(:, bloodIndx);
bloodDev = bloodSample - mean(bloodSample);

clSample = postSamples(:,CL_antiPDL1Indx);
clDev = clSample - mean(clSample);

postSamples(:, bloodIndx) = bloodDev + log(0.0117);
postSamples(:,CL_antiPDL1Indx) = clDev + log(0.0378); 

% Adjusting Q12
% q12Indx = ismember(parameters(1:8), 'Q12');
% peripheralIndx = ismember(parameters(1:8), 'Peripheral');
% postSamples(:, q12Indx) = postSamples(:,q12Indx) - log(10);
%%
simTime = unique([PI.tspan', 1:PI.tspan(end)]);
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,simTime),x,...
    @(p)getPhi2(p,PI.H,length(PI.u),'initialValue',PI.x_0),PI.normIndx, PI.H,...
    'output', 'simoutput', 'simTime', simTime);
tic
PI=getPosteriorPredictions2(exp(postSamples(1:end,:)),PI,simFun,PI.observablesFields,...
    'simTime', simTime);
toc
PI=getCredibleIntervals(PI,PI.observablesFields, exp(postSamples(1:end,:)),PI.H,...
    'logit_indx', [],'simTime', simTime);
%% Plotting individual predictions

for i =1:length(PI.observablesPlot)
 plotPosteriorPredictions(PI,i,'outputs','indiv', ...
        'newFig', true, 'TimeUnit', 'hours','color', 'cell','simTime', ...
        simTime, 'YScale', 'log')

end
%% Save results
save('PK_Validation_ThreeComp4_9_TMDD_8.mat', 'PI')
load(strjoin({cd 'PI_PK_ThreeComp4_7_TMDD_12.mat'},'/'))

save(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_x2.mat'},''), 'x2')
save(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_p_x2.mat'},''), 'p_x2')
load(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_p_x1.mat'},''))
load(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_11_DREAM_MCMC_x1.mat'},''))