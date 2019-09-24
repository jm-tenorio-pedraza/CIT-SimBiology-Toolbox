%% Kuznetsov general setup

%% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/MATLAB/SimBiology/Kuznetsov/output/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kuznetsov.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.25)
if sensitivity
    [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
    parameters=name(value>0);
    % Define parameters to exclude from SA
    exclude_parameters = {'k12' 'k21' 'K_D_antiPDL1' 'k23' 'k32' 'vol_Tumor'...
        'T_0'  'CD8' 'CD8_E' 'Tumor_P' 'TV' 'kde_Tumor' 'CD8_logit'};
    parameters = setdiff(parameters, [exclude_parameters]);

else
    % Define parameters to estimate
    parameters = load('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kuznetsov/parameters_hat.mat');
    parameters = parameters.parameters_hat;
end
% Define outputs
observables={'TV' 'CD8_logit'};
stateVar={'Tumor' 'CD8_logit'};
groups_subset = {'MOC1_Control' 'MOC1_antiPDL1' 'MOC1_Control_Mean' 'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1'};
doses = {'Dose_antiPDL1' 'Dose_antiCTLA4'};
PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',stateVar,groups_subset,observables);
PI.variableUnits={'Volume [mL]' 'Percentage [%]'};

[sim,u]=initializePI(model,parameters,observables,PI,doses, 'MOC1');

%% Optimization setup
% Hierarchical structure
H.PopulationParams=1:length(parameters);
H.IndividualParams=struct('name', parameters(1:2),...
    'Index', {length(parameters)+1:length(parameters)+length(u);...
    length(parameters)+length(u)+1:length(parameters)+length(u)*2},...
    'EtaIndex', {1;2},...
    'OmegaIndex', {length(parameters)+length(u)*2+1;length(parameters)+length(u)*2+2});
H.SigmaParams=length(parameters)+length(u)*2+1:length(parameters)+length(observables)+length(u)*2+2;
% Generating PI
sigmaNames={'Omega_kpro_max_effector';'Omega_kin_max';'b_TV'; 'b_CD8'};
sigma_prior= [ repelem(1,length(H.PopulationParams), 1);...
     repelem(1, length([H.IndividualParams(:).Index]),1);...
    repelem(1, length(H.SigmaParams),1)];
PI.par = getParamStruct2(sim,H,size(PI.data,1),repelem(0.5,length(H.SigmaParams),1),...
    sigmaNames,'Sigma', sigma_prior);

% Residuals function
residuals_fun=@(p)getNormResiduals(p,@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,size(PI.data,1)),(@(x)getCovariance(x,H)),[]);

% Log-ikelihood function
likelihood_fun=@(p)sum(residuals_fun(exp(p))*(-1));
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,H,'type','uniform')-norm(p));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,100,u,1:1:100),PI,...
    @(x)getPhi2(x,H,length(u)),finalValues(end-length(observables)+1:end),[]);

%% Save results
save('PI_Kuznetsov.mat', 'PI')
save('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kuznetsov/parameters_hat.mat','parameters_hat')
