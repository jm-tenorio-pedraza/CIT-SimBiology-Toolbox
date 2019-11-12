%% Kuznetsov general setup

%% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kuznetsov/PI')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kuznetsov.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)
%% Load parameters to estimate and data
if sensitivity
    [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
    parameters=name(value>0);
    % Define parameters to exclude from SA
    exclude_parameters = {'k12' 'k21' 'K_D_antiPDL1' 'k23' 'k32' 'vol_Tumor'...
        'T_0'  'CD8' 'CD8_E' 'Tumor_P' 'TV' 'kde_Tumor' 'CD8_logit' 'Q12' 'Q23' ...
        'ka_Central' 'ke_Central_antiCTLA4' 'ke_Central_antiPDL1' 'E_0'};
    parameters = setdiff(parameters, [exclude_parameters]);
    parameters = [ parameters; 'T_0'];
else
    % Define parameters to estimate
    parameters = load(strjoin({cd 'parameters_hat_2.mat'},'/'));
    parameters = parameters.parameters_hat;
    parameters = setdiff(parameters, 'kel_Effector','stable');
    parameters = [parameters; 'T_0'];
%     parameters = [{'kpro_Tumor'}; setdiff(parameters,'kpro_Tumor','stable')];
end
% Define outputs
observables={'TV' 'CD8_logit'};
stateVar={'Tumor' 'CD8_logit'};
groups_subset = {'MOC1_Control' 'MOC1_antiPDL1' 'MOC1_Control_Mean'...
    'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1' 'MOC2_Control' 'MOC2_Control_Mean'...
    'MOC2_antiPDL1' 'MOC2_antiCTLA4'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction','input','mergePhenotypes',true);
PI.variableUnits={'Volume [mL]' 'Logit []'};

[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1');
normIndx = [];
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1});
% Get initial values
x_0 = getInitialValues([PI.data(:).Group], initialStruct);
% close all
%% Optimization setup
% Hierarchical structure
H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx',1 , 'n_indiv', length(u),'cell_indx', 3);
SigmaNames=[];
try
    cellSigmaNames=arrayfun(@(x)strjoin({'lambda', x.name}, '_'),H.CellParams,'UniformOutput',false)';
    SigmaNames = [SigmaNames;cellSigmaNames];
    indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false)';
    SigmaNames = [SigmaNames;indivSigmaNames];
catch
    try
        indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),H.IndividualParams,'UniformOutput',false)';
        SigmaNames = [SigmaNames;indivSigmaNames];
    catch
    end
    
end
SigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'sigma', x}, '_'),...
        observables,'UniformOutput', false);
% Generating PI
alpha = [repelem(0.01, length([H.CellParams(:).OmegaIndex]),1);...
    repelem(0.01, length([H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(H.SigmaParams, [H.CellParams(:).OmegaIndex ...
    H.IndividualParams(:).OmegaIndex])),1)];
beta = [repelem(0.1, length([H.CellParams(:).OmegaIndex]),1);...
    repelem(0.1, length([H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(H.SigmaParams, [H.CellParams(:).OmegaIndex ...
    H.IndividualParams(:).OmegaIndex])),1)];
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
    repelem(1, length([H.CellParams(:).Index]),1);
     repelem(1, length([H.IndividualParams(:).Index]),1);...
    alpha];
PI.par = getParamStruct2(sim,H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior);
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,PI.tspan(end),u,PI.tspan),PI,...
    @(x)getPhi2(x,H,size(PI.data,1),'initialValue',x_0),(@(x)getCovariance(x,H)),normIndx);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,H,'type','uniform'));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),u,PI.tspan),PI,...
    @(x)getPhi2(x,H,length(u),'initialValue',x_0),exp(finalValues(end-length(observables)+1:end)),normIndx);
paramNames = ['\eta_{kin_{max}}' 'kill_{max}'...
    '\eta_{kpro_{Tumor}}' 'f2'...
    {PI.par([H.CellParams(:).Index H.IndividualParams(:).Index]).name}, '\lambda_{kin_{max}}',...
    '\omega_{kpro_{Tumor}}', '\sigma_{TV}'  '\sigma_{CD8}'];
%% Save results
save('PI_Kuznetsov_4.mat', 'PI')
load(strjoin({cd 'PI_Kuznetsov_4.mat'},'/'))
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kuznetsov/parameters_hat.mat','parameters_hat')
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))