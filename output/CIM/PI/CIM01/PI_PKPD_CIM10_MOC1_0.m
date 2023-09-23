%% General setup
%% Search paths
clear all
warning off
opengl('save', 'software')
%%
if ispc
    wd='\Users\jmten\OneDrive\Dokumente\GitHub\';
else
    wd='/Users/migueltenorio/Documents/GitHub/';
end
    addpath(genpath(strjoin({wd 'CIT-SimBiology-Toolbox'}, '')))
    cd(strjoin({wd 'CIT-SimBiology-Toolbox/output'},''))
    out=sbioloadproject(strjoin({wd 'QSP-models/CIM_10_1_PKPD.sbproj'},''));
        data_ext = strjoin({ wd 'CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'},'/');
    data_ext1 = strjoin({wd 'CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'},'/');

%% Load project 
% Extract model
model=out.m1;
variants = getvariant(model);
initialStruct = struct('name', {'MOC1';'MOC2';'MC38'}, 'initialValue', {5; 0.1; 0.1},...
    'variant', {variants(1); variants(3); variants(4)});

cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-9);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-6);
set(cs, 'MaximumWallClock', .5)
set(cs, 'Time','day')
set(cs, 'StopTime',100)
    %% Rule setup
% Initial tumor volume
model.Rules(21,1).Active=true;
% Initial tumor cell number
model.Rules(1,1).Active=true;
% Initial PDL1_Tumor
model.Rules(18,1).Active=true;
% Initial CTLA4_CD8
model.Rules(19,1).Active=true;
% Initial CTLA4_Treg
model.Rules(22,1).Active=true;
% Initial PDL1_Immune
model.Rules(32,1).Active=true;
% Carrying capacity
model.Rules(50,1).Active=false;
% Observers
model.Rules(2,1).Active=true;
model.Rules(3,1).Active=true;
model.Rules(4,1).Active=true;
model.Rules(5,1).Active=true;
model.Rules(6,1).Active=true;
model.Rules(9,1).Active=true;
model.Rules(10,1).Active=true;
model.Rules(12,1).Active=true;
% saturation_antiCTLA4
model.Rules(51,1).Active = false;
%% Parameter setup
parameters = {'kin_CD8'; 'kin_Treg'; 'kin_DC';'kin_MDSC';...
    'kpro_Tumor';'kill_CD8'; ...
    'PDL1_Tumor_ss';'PDL1_Immune_ss';
    'K_MDSC';'kel_DC'; 'kel_Effector';'kel_MDSC'; 'kel_Treg';};
parameters = [parameters;'T_0'];
% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_Control_Mean',  'MOC1_Control_Immune' ...
    'MOC1_antiCTLA4' 'MOC1_antiPDL1' 'MOC1_antiCTLA4_antiPDL1'};
observables={'Tumor'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};
stateVar={'Tumor'  'CD8' 'Treg' 'DC'...
    'GMDSC' 'CD107a' 'Tumor_PDL1_Rel' 'Myeloid_PDL1_Rel'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};
%% Obtain data, simulation function and dose table
if ispc
    data_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'};
    data_ext1 = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'};
else
    data_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo_2.mat'};
    data_ext1={'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada_3.mat'};
end
PIClavijo=getPIData4(data_ext, stateVar,groups_subset,'output', 'mean',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'replace','logTransform',false);
PIMorisada=getPIData4(data_ext1, stateVar,groups_subset,'output', 'mean',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'replace','logTransform',false);

PIIndiv_Clavijo=getPIData4(data_ext, stateVar,groups_subset,'output', 'individual',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'noAction','logTransform',false);
PIIndiv_Morisada=getPIData4(data_ext1, stateVar,groups_subset,'output', 'individual',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'replace','logTransform',false);
% Remove extreme points
PIMorisada.data(1).dataTime = PIMorisada.data(1).dataTime(1:end-3);
PIMorisada.data(1).dataValue = PIMorisada.data(1).dataValue(1:end-3,:);
PIMorisada.data(1).SD = PIMorisada.data(1).SD(1:end-3,:);
PIMorisada.data(1).zero_indx = PIMorisada.data(1).zero_indx(1:end-3);
PIMorisada.n_data = sum(cellfun(@(x)sum(sum(~isnan(x))), {PIMorisada.data(:).dataValue}));

Morisada_Names=({PIIndiv_Morisada.data(:).Name}');
morisada_unique = {PIMorisada.data(:).Name};
morisadaData = PIIndiv_Morisada.data;
structIndx=1;
for i=1:length(morisada_unique)
    index=ismember(Morisada_Names,morisada_unique(i));
    morisadaData(structIndx:structIndx+sum(index)-1) = PIIndiv_Morisada.data(index);
    morisada_counts(i) = sum(index);
        structIndx = structIndx+sum(index);

end

Clavijo_Names =categorical({PIIndiv_Clavijo.data(:).Name}');
clavijo_unique = {PIClavijo.data(:).Name};
% Exclude responders
clavijoData = PIIndiv_Clavijo.data(ismember(Clavijo_Names,clavijo_unique));
structIndx=1;
% For each unique group count how many mice are in it
for i=1:length(clavijo_unique)
    index = ismember(Clavijo_Names,clavijo_unique(i));
    clavijoData(structIndx:structIndx+sum(index)-1)= PIIndiv_Clavijo.data(index);
    clavijo_counts(i)=sum(index);
    structIndx = structIndx+sum(index);
end
PIClavijo.n_data = sum(cellfun(@(x)sum(sum(~isnan(x))), {clavijoData(:).dataValue}));

%% Assemble PI structure
PI=[];
PI.data=[PIClavijo.data; PIMorisada.data];
PI.n_data=sum(cellfun(@(x)sum(sum(~isnan(x))), {PI.data(:).dataValue}));
PI.tspan=unique([PIClavijo.tspan; PIMorisada.tspan]);
%Rename Response vector for control groups
ImmuneResp={'Control_Clavijo' 'Control_Clavijo' 'Control_Morisada' 'Control_Morisada'};
[PI.data([1:2 6:7]).Response] = ImmuneResp{:,:};
% Add cell field
Cell_Field = cellfun(@(x)x(1:find(x=='_',1)-1),[PI.data(:).Group],...
    'UniformOutput',false);
[PI.data(1:end).Cell] = Cell_Field{:,:};
% Add counts of individual data corresponding to each group
counts=[clavijo_counts';morisada_counts'];
counts=mat2cell(counts,ones(length(counts),1));
[PI.data(1:end).Count] = counts{:,:};
PI.IndivData=[clavijoData; morisadaData];
PI.indivN_data =sum(cellfun(@(x)sum(sum(~isnan(x))), {PI.IndivData(:).dataValue}));
clearvars PIClavijo PIIndiv PIIndiv_Clavijo PIIndiv_Morisada PIMorisada Morisada_Names Clavijo_Names morisada_counts morisada_unique morisadaData clavijo_counts Clavijo_Names clavijo_unique clavijoData counts
%%
PI.variableUnits={'Volume [mL]' 'Percentage [%]' 'Percentage [%]'  'Percentage [%]' ...
     'Percentage [%]'   'Percentage [%]'  ...
    'Relative units []' 'Relative units []'};
PI.observablesFields = {'TV'  'CD8'  'Treg' 'DCm'...
    'MDSC' 'CD8_E' 'PDL1_T' 'PDL1_I'};
PI.normIndx = 7:8;
PI.model = 'CIM10_PKPD_MOC1_kill_CD8';
PI.observablesPlot={'TV' 'CD8' 'Treg' 'DCm'...
    'MDSC' 'CD107a' 'PDL1_T' 'PDL1_I'};
plotData(PI, PI.observablesPlot, 'responseGrouping', true, 'kineticGrouping', false)
% Get initial values
[PI.x_0, PI.variants] = getInitialValues([PI.data(:).Group],...
    initialStruct);
PI.varDist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
%% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'PKPD_Fit_MOC1','doseUnits', 'mole');
clearvars   PI1 PI2 variants data_ext data_ext1 doses ans Cell_Field groups_subset groupsResp ImmuneResp index
% close all
%% Optimization setup
% Hierarchical structure
PI.H = getHierarchicalStruct3(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [] , 'cell_indx',[], 'resp_indx', [],'n_indiv', length(PI.u),...
    'CellField','Cell');
SigmaNames = getVarNames(PI, stateVar);
[beta, sigma_prior] = getVarValues2([.1 .1 .1 .1], ...
    [.6 .6 .6 .6], [1 1 1 1], PI);
% lb=([1e-2   1e-2    1e-2    1e-2    .01     1e-3    1e-3    1e-4  1e-4       1e-2   1e-3     1e-3    1e-2])';
% ub=([1e4    1e2     1e2     1e5     1       100     1e4     1e6   1e6        1e5    1e3      1e3     1e4])';
paramValues = sim.Parameters.Value(1:end-1);
lb = paramValues.*[1e-6  1e-6    1e-6   1e-6    0.5     1e-6    1e-6    1e-6    1e-10  1e-6    1e-6    1e-6    1e-6]';
ub = paramValues.*[1e6    1e6      1e6     1e6      2       1e6      1e6      1e6      1e6    1e6      1e6     1e6     1e6]';

initialSigma=ones(length(PI.H.SigmaParams),1);
initialSigmaVars=ones(length(PI.varDist),1);
initialSigmaVars(ismember(PI.varDist,'Beta')) =.1;
initialSigma(end-(length(PI.varDist)-1):end)=initialSigmaVars;
PI.par = getParamStruct3(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'startSigma',...
    initialSigma, 'ref', 'ones','LB', lb, 'UB', ub);

PI = assignPrior(PI,'sigmaDist', 'JP');
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false,...
    'logTransform',true,'errorModel','additive','constantVar',.01,'indivData',true);
prior_fun=@(p)getPriorPDFMCMC2(exp(p),PI);


paramNames = getParamNames(PI,sim, observables);
PI.paramNames  =paramNames;
clearvars beta lb ub sigma_prior SigmaNames 
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

 %% Save results
% save(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_MCMC_kin_CD8_2.mat'},'/'), 'PI')
% %% Load results
% load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_MCMC_kin_CD8_2.mat'},'/'),'PI')
