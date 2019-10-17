%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreComp/PI')
sensitivity = false;
PI=1;
%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/PK_mAb_ThreeComp.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)
[name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
    parameters=name(value>0);
    if PI ==1
exclude_parameters = {'k12' 'k21' 'k23' 'k32' 'ID_g_Tumor'...
        'T2B'  'ID_g_Blood' 'ka_Tumor' 'ka_IP', 'ID','ke_Central','ke_Tumor'};
    else
        exclude_parameters = {'k12' 'k21' 'k23' 'k32' 'ID_g_Tumor'...
        'T2B'  'ID_g_Blood' 'ka_Tumor' 'ka_IP', 'ID','ke_Central'};

    end
parameters = setdiff(parameters, [exclude_parameters]);
parameters = ['Blood'; 'Tumor';'ke_Central';parameters; 'IP'; 'ID'];
% Define outputs
observables={'ID_g_Blood' 'Blood.antiPDL1' 'ID_g_Tumor' 'Tumor.antiPDL1' 'T2B'};
stateVar={'Tumor.antiPDL1' 'Tumor_to_Blood' 'Blood.antiPDL1'};

dataset_file_ext = {'/Users/migueltenorio/Documents/Data/Nedrow_2017_1.xlsx'...
    '/Users/migueltenorio/Documents/Data/Nedrow_2017_2.xlsx'...
    '/Users/migueltenorio/Documents/Data/Contreras_2016.xlsx'};

[PI,u]=getDataSets(dataset_file_ext);
PI.variableUnits={'%ID/g' 'mg/l' '%ID/g' 'mg/l' '[]'};
variants=getvariant(model);
MOC1=variants(strcmp(get(variants,'Name'), 'MOC1'));
dose = {'Blood.Dose_antiPDL1'};
sim=createSimFunction(model,parameters,observables, dose,MOC1,...
    'UseParallel', false);
normIndx = [];
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1});
% Get initial values
x_0 =[PI.data(:).dose]';

%% Optimization setup
% Hierarchical structure
H = getHierarchicalStruct(parameters(1:end-1),'n_sigma', length(observables), 'n_rand', 2, 'n_indiv', length(u));
try
    sigmaNames=arrayfun(@(x)strjoin({'Omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false)';
    sigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'b', x}, '_'),observables,'UniformOutput', false);
catch
    sigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables','UniformOutput', false);
end
% Generating PI
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
     repelem(1, length([H.IndividualParams(:).Index]),1);...
    repelem(0.1, length(H.SigmaParams),1)];
PI.par = getParamStruct2(sim,H,size(PI.data,1),repelem(0.5,length(H.SigmaParams),1),...
    sigmaNames,'Sigma', sigma_prior);
 finalValues =log([PI.par(:).startValue]);

% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,144,u,1:.1:144),PI,...
    @(x)getPhi2(x,H,length(u),'initialValue',x_0),(@(x)getCovariance(x,H)),normIndx);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,H,'type','uniform'));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,144,u,1:.1:144),PI,...
    @(x)getPhi2(x,H,length(u),'initialValue',x_0),exp(finalValues(end-length(observables)+1:end)),normIndx);

%% Save results
save('PI_PK.mat', 'PI')
load(strjoin({cd 'PI_PK.mat'},'/'))
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kuznetsov/parameters_hat.mat','parameters_hat')
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))