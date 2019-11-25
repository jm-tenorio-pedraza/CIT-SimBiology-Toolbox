%% General setup for Kosinsky et al with k_pro as the parameter to vary for each individual
%% Search paths
warning off
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/Kosinsky/PI')

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_modified.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', 0.2)
sensitivity = false;
%% load data and previous results
stateVar={'Tumor' 'CD8_logit' 'CD107a_logit' 'DC_Rel' 'GMDSC_Rel'...
    'Tumor_PDL1_Rel'};

% Create function handle for simulations
% Define parameters to estimate
if sensitivity
    [name,I] = sort(get(model.Parameters, 'Name'));
    value = cell2mat(get(model.Parameters, 'Value'));
    value = value(I);
    % Define parameters to estimate
    parameters=name(value>0);
    exclude_parameters = {'CD107a' 'CD8' 'K_D_antiPDL1' 'K_D_antiCTLA4'...
        'Total_Cell_Count' 'vol_Tcell' 'vol_Tumor', 'T_0' 'CD107a_logit'...
        'CD8_logit' 'kin_Naive' 'kin_Effector' 'ka_Central'};
    parameters = setdiff(parameters, exclude_parameters);
    parameters = [ parameters; 'T_0'];
else
    parameters = {'kin_max'; 'kdif'; 'K_Naive'; 'f3';'kpro_Tumor';'kill_max'};
    parameters = [parameters; 'T_0'];

end
% Define outputs% Define outputs
observables={'TV' 'CD8_logit' 'CD107a_logit' 'DCm' 'ISC' 'PDL1'};

% Create PI with data
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean', 'MOC1_antiPDL1', ...
    'MOC1_antiCTLA4' 'MOC1_antiCTLA4_antiPDL1', 'MOC2_Control', 'MOC2_Control_Mean',...
    'MOC2_antiPDL1', 'MOC2_antiCTLA4'};

doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table
PI=getPIData('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables, 'zeroAction', 'omit','mergePhenotypes', true);

PI.variableUnits={'Volume [mL]' 'Logit []' 'Logit []' 'Relative units []'...
    'Relative units []' 'Relative units []'};

variants = getvariant(model);

PI.normIndx = 4:6;
initialStruct = struct('name', {'MOC1';'MOC2'}, 'initialValue', {5; 0.1},'variant', {variants(2); variants(3)});
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group], initialStruct);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'MOC1','doseUnits', 'mole');


%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', 1, 'cell_indx',5, 'n_indiv', length(PI.u));
try
    cellSigmaNames=arrayfun(@(x)strjoin({'lambda', x.name}, '_'),PI.H.CellParams,'UniformOutput',false)';
    indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),PI.H.IndividualParams,'UniformOutput',false)';
    SigmaNames = [cellSigmaNames; indivSigmaNames];
    SigmaNames(end+1:end+length(observables),1) =  cellfun(@(x) strjoin({'sigma', x}, '_'),...
        observables,'UniformOutput', false);
catch
    try
    indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),PI.H.IndividualParams,'UniformOutput',false)';
    
    SigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables','UniformOutput', false);
    SigmaNames=[indivSigmaNames; SigmaNames];
    catch
            SigmaNames= cellfun(@(x) strjoin({'b', x}, '_'),observables','UniformOutput', false);
    end
end
% Generating PI
alpha = [repelem(0.01, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.01, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
beta = [repelem(0.1, length([PI.H.CellParams(:).OmegaIndex]),1);...
    repelem(0.1, length([PI.H.IndividualParams(:).OmegaIndex]),1);...
    repelem(0.01, length(setdiff(PI.H.SigmaParams, [PI.H.CellParams(:).OmegaIndex ...
    PI.H.IndividualParams(:).OmegaIndex])),1)];
sigma_prior= [ repelem(1,length(PI.H.PopulationParams), 1);...
    repelem(1, length([PI.H.CellParams(:).Index]),1);
     repelem(1, length([PI.H.IndividualParams(:).Index]),1);...
    alpha];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior);
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end

% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun=@(p)(createPriorDistribution3(exp(p),PI,PI.H,'type','uniform'));

% Residuals 
residuals_fn = @(x) getResiduals(exp(x),@(x)sim(x,PI.tspan(end),PI.u,PI.tspan),PI,...
    @(x)getPhi2(x,PI.H,length(PI.u),'initialValue',PI.x_0),exp(finalValues(end-length(observables)+1:end)),PI.normIndx);

%% Save results
save('PI_Kosinsky_4.mat', 'PI')
load(strjoin({cd 'PI_Kosinsky_4.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_p.mat'},'/'))
load(strjoin({cd 'DREAM_MCMC_logP.mat'},'/'))

sbiosaveproject '/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/Kosinsky_2.sbproj' 'model'