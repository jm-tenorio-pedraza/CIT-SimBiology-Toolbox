%% Initial values
% Search paths
clear all
warning off
if ispc
    addpath(genpath('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox'))
    cd('\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\output')
    out = sbioloadproject('\Users\jmten\OneDrive\Dokumente\GitHub\sbio-projects\CIM_10_1_PKPD.sbproj');
    data_ext = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'};
    data_ext1 = {'\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'};
else
    addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
    cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output')
    out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/sbio-projects/CIM_10_1_PK.sbproj');
    data_ext = {'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo_2.mat'};
    data_ext1={'/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Morisada_3.mat'};
end
%% Load project 
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 2.5)
MODEL = 'CIM 10';
variants = get(model, 'variants');
set(cs, 'Time','day')
set(cs, 'StopTime',365*3)

model.Rules(21,1).Active=true;
model.Rules(1,1).Active=true;
model.Rules(18,1).Active=true;
model.Rules(19,1).Active=true;
model.Rules(23,1).Active=true;
model.Rules(33,1).Active=true;

%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Peripheral';'CL_antiPDL1';'Q12';...
    'f_L'; 'V_vs';'kin_CD8'; 'kin_Treg'; 'kin_DC';'kin_MDSC';'kin_CD45';...
    'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; ...
    'K_MDSC';'K_IFNg';'K_CTLA4'; 'K_PDL1'; 'KDE_MDSC';'KDE_Treg';
    'V_intra';'PDL1_Tumor_ss';'PDL1_Immune_ss';'kel_CD45';'kel_DC';'kel_Effector';'kel_MDSC';...
    'kel_Naive';'kel_Treg';'T_0'};
% Define outputs
groups_subset = {'MOC1_Control'};
observables={'Tumor'   'TumorInter.Cancer'  'TumorInter.CD8_E' 'TumorInter.CD8_N'...
    'TumorInter.DC' 'TumorInter.GMDSC'...
    'TumorInter.CD4Foxp3' 'TumorInter.Ag'  'TumorInter.CD45'...
    'TumorInter.CTLA4_CD8' 'TumorInter.CTLA4_Treg'...
    'TumorInter.PDL1_Tumor' 'TumorInter.PDL1_Immune'...
    };
stateVar={'Tumor'};
doses = {'Blood.antiPDL1' 'Blood.antiCTLA4'};
%%
PI=getPIData4(data_ext, stateVar,groups_subset,'output', 'mean',...
    'responseGrouping', true, 'kineticGrouping', false, 'zeroHandling',...
    'replace','logTransform',false);
cellgroup=repelem({'MOC1'},1);
[PI.data(1:end).Cell]=cellgroup{:,:};
PI.data(1).Group='MOC1_Control';
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'Human simulations';
PI.observablesPlot={'Tumor' 'Cancer' 'CD8_E' 'CD8_N' 'DC' 'GMDSC'...
    'CD4Foxp3' 'Ag' 'CD45' 'CTLA4_CD8' 'CTLA4_Treg' 'PDL1_Tumor' 'PDL1_Immune'};
PI.tspan = 1:14:(365*3);

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'Human_InitialValues','doseUnits', 'mole');
PI = getDoseRegimen(PI);

%% Hierarchical model simulation
PI.H = getHierarchicalStruct3(parameters(1:end),PI,'n_sigma', length(observables),...
    'rand_indx', [] , 'cell_indx',[], 'resp_indx', [],'n_indiv', length(PI.u),...
    'CellField','Cell');% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues2([.1 .1 .1 .1], ...
    [.6 .6 .6 .6], [1 1 1 1], PI);
PI.par = getParamStruct3(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'startSigma',...
    ones(length(PI.H.SigmaParams), 1)*.6, 'ref', 'ones');
try
    finalValues =log([PI.par(:).finalValue]);
catch
    finalValues =log([PI.par(:).startValue]);

end
PI= assignPrior(PI);
% Log-ikelihood function
likelihood_fun=@(p)likelihood(exp(p),sim,PI,'censoring',false);
prior_fun_MCMC=@(p)getPriorPDFMCMC2(exp(p),PI);
paramNames = getParamNames(PI,sim, observables);
%% Load PIs and posterior parameter samples
PI_PK = load(strjoin({cd 'PK_mAb_ThreeComp/PI/CIM10/PI_CIM10_PK_1.mat'},'/'));
PI_MOC1 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC1_PKPD_1_8_kinCD8.mat'},'/'));
PI_MOC2 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC2_PKPD_1_8.mat'},'/'));
PI_MC38 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MC38_PKPD_1_8.mat'},'/'));

Meta=[];
Meta(1).Struct = PI_MC38;
Meta(2).Struct = PI_MOC1;
Meta(3).Struct = PI_MOC2;
values=repelem({sim.Parameters.Value},length(Meta),1);
[Meta(1:end).Values]=values{:,:};
%% Select which dimensions to sample from 
N_pop =1e3;
[theta1,~] = getTheta2(Meta, parameters, N_pop, 'covStructure', 'dependent','fixedParams','fixed');
[theta2,~] = getTheta2(Meta, parameters, N_pop, 'covStructure', 'independent');
% Indexes of parameter uncertainties to propagate
PK_params={'Blood' 'Peripheral' 'CL_antiPDL1'...
    'Q12' 'f_L' 'V_vs'};
PD_params = {'kpro_Tumor'};
indx_pk = find(ismember(parameters,PK_params));
indx_pd = find(ismember(parameters, PD_params));
priorMu_pk = log([PI.par(indx_pk).startValue]);
priorMu_pd = log([PI.par(indx_pd).startValue]);

Theta = theta1;
histogram(Theta(:,ismember(parameters,'kin_CD8')))
% PK parameters
[pk_indx, pk_order]=ismember({PI_PK.PI.par(:).name},PK_params);
pk_order=pk_order(~pk_order==0);
theta_pk = PI_PK.PI.postSamples(:,ismember({PI_PK.PI.par(:).name},PK_params));
sampleIndx_pk=randsample(size(theta_pk,1),N_pop*3, true);
delta_pk = theta_pk(sampleIndx_pk,:)-mean(theta_pk);
Theta(:,indx_pk) =priorMu_pk+delta_pk(:,pk_order);

% PD parameters
% % [pd_indx, pd_order]=ismember(parameters,PD_params);
% % pd_order=pd_order(~pd_order==0);
% % theta_pd = Theta(:,pd_indx);
% % delta_pd = theta_pd-mean(theta_pd);
% % Theta(:,indx_pd) =priorMu_pd + delta_pd(:,pd_order);
histogram(exp(Theta(:,indx_pd)))
% Sample initial tumor cell values
T0_range= [0.1 5]; %  0.5 to 5 cm
T_0 = log((rand(size(Theta,1),1)*(T0_range(2)-T0_range(1))... % Initial T_0 (1e6 cells)
    +T0_range(1)));
Theta(:,ismember(parameters,'T_0')) = T_0;
histogram(exp(T_0))
% Sample baseline tumor values
TV_range = [.5 2.5].^3*(4/3)*pi;
TV_0 = log((rand(size(Theta,1),1)*(TV_range(2)-TV_range(1))... % Initial T_0 (1e6 cells)
    +TV_range(1)));
histogram(exp(TV_0))
table(parameters,exp(mean(Theta))')
%% Set linear tumor growth to 0
% kproTumorLinear_indx=ismember(parameters,'kpro_Tumor_Linear');
% Theta(:,kproTumorLinear_indx)=-log(1000);
%% Scaling parameters
a=1;    b=.9;  c=.8;   d=.75;  e=.65; f=.25; g=-.25;

scalingExp1 =   [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % none
scalingExp2 =   [0 0 0 0 0 0 d 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % CD8 T cell infiltration at .75
scalingExp3 =   [0 0 0 0 0 0 d 0 0 0  0  0  0  d  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % CD8 T-Cell infiltration at .75 and killing rate at .75
scalingExp4 =   [0 0 0 0 0 0 d 0 0 0  0  0  0  d  0  0  0  0  0  0  0  0  0  0  0  g  0  0  0  0]; % CD8 T-Cell infiltration at .75, killing rate at .75, naive T-cell death rate at -.25
scalingExp5 =   [0 0 0 0 0 0 a 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % CD8 T cell infiltration at 1
scalingExp6 =   [0 0 0 0 0 0 a 0 0 0  0  0  0  a  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % CD8 T cell infiltration at 1 and killing rate at 1
scalingExp7 =   [0 0 0 0 0 0 a 0 0 0  0  0  0  a  0  0  0  0  0  0  0  0  0  0  0  g  0  0  0  0]; % CD8 T-Cell infiltration at 1, killing rate at 1, naive T-cell death rate at -.25
scalingExp8 =   [0 0 0 0 0 0 d d 0 0  0  0  0  d  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % T cell infiltration at .75 and killing rate at .75
scalingExp9 =   [0 0 0 0 0 0 a a 0 0  0  0  0  a  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % T cell infiltration at 1 and killing rate at 1
scalingExp10=   [0 0 0 0 0 0 d d 0 0  0  0  0  d  0  0  0  0  0  0  0  0  0  0  0  g  0  g  g  0]; % T cell infiltration at .75 and killing rate at .75 and T-cell death rate at -.25
scalingExp11 =  [0 0 0 0 0 0 a a 0 0  0  0  0  a  0  0  0  0  0  0  0  0  0  0  0  g  0  g  g  0]; % T cell infiltration at .1 and killing rate at 1 and T-cell death rate at -.25
scalingExp12 =  [0 0 0 0 0 0 d d d d  d  0  0  d  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % Cell infiltration at .75 and killing rate at .75
scalingExp13 =  [0 0 0 0 0 0 a a a a  a  0  0  a  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]; % Cell infiltration at 1 and killing rate at 1
scalingExp14 =  [0 0 0 0 0 0 d d d d  d  0  0  d  0  0  0  0  0  0  0  0  0  g  g  g  g  g  g  0]; % Cell infiltration at .75 and killing rate at .75 and death rate at -.25
scalingExp15 =  [0 0 0 0 0 0 a a a a  a  0  0  a  0  0  0  0  0  0  0  0  0  g  g  g  g  g  g  0]; % Cell infiltration at 1 and killing rate at 1 and death rate at -.25
scalingExp16 =  [0 0 0 0 0 0 d 0 0 0  0  0  0  d  a  a  0  0  a  a  0  0  0  0  0  0  0  g  0  0]; % CD8 T cell infiltration at .75, cell amounts at 1, killing rate at .75 and naive T-cell death rate at -.25

scalingFactor1 = (77/.022).^(scalingExp1);
scalingFactor2 = (77/.022).^(scalingExp2);
scalingFactor3 = (77/.022).^(scalingExp3);
scalingFactor4 = (77/.022).^(scalingExp4);
scalingFactor5 = (77/.022).^(scalingExp5);
scalingFactor6 = (77/.022).^(scalingExp6);
scalingFactor7 = (77/.022).^(scalingExp7);
scalingFactor8 = (77/.022).^(scalingExp8);
scalingFactor9 = (77/.022).^(scalingExp9);
scalingFactor10 = (77/.022).^(scalingExp10);
scalingFactor11 = (77/.022).^(scalingExp11);
scalingFactor12 = (77/.022).^(scalingExp12);
scalingFactor13 = (77/.022).^(scalingExp13);
scalingFactor14 = (77/.022).^(scalingExp14);
scalingFactor15 = (77/.022).^(scalingExp15);
scalingFactor16 = (77/.022).^(scalingExp16);

Theta1 = Theta + log(scalingFactor1);
Theta2 = Theta + log(scalingFactor2);
Theta3 = Theta + log(scalingFactor3);
Theta4 = Theta + log(scalingFactor4);
Theta5 = Theta + log(scalingFactor5);
Theta6 = Theta + log(scalingFactor6);
Theta7 = Theta + log(scalingFactor7);
Theta8 = Theta + log(scalingFactor8);
Theta9 = Theta + log(scalingFactor9);
Theta10= Theta + log(scalingFactor10);
Theta11= Theta + log(scalingFactor11);
Theta12= Theta + log(scalingFactor12);
Theta13= Theta + log(scalingFactor13);
Theta14= Theta + log(scalingFactor14);
Theta15= Theta + log(scalingFactor15);
Theta16= Theta + log(scalingFactor16);
%% Plot bivariate marginals
plotBivariateMarginals_2(exp(Theta1(:,:)),...
       'names',parameters,'interpreter', 'tex')
plotCorrMat(Theta1, parameters)
%Plot 1st cell line samples
plotBivariateMarginals_2((theta1(1:N_pop,:)),...
       'names',parameters,'interpreter', 'tex')

%% Simulate
simFun=@(x)getOutput2(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    {},PI.normIndx,'Output', 'data');
dataOutput = simFun(exp(mean(Theta1(N_pop+100,:),1)));
[PI.data(1:end).simValue] = dataOutput{:,:};
plotSimValue(PI,'scale','linear')

PI1=getPosteriorPredictions2(exp(Theta1(:,:)),PI,simFun,PI.observablesPlot);
PI2=getPosteriorPredictions2(exp(Theta2(:,:)),PI,simFun,PI.observablesPlot);
PI3=getPosteriorPredictions2(exp(Theta3(:,:)),PI,simFun,PI.observablesPlot);
PI4=getPosteriorPredictions2(exp(Theta4(:,:)),PI,simFun,PI.observablesPlot);
PI5=getPosteriorPredictions2(exp(Theta5(:,:)),PI,simFun,PI.observablesPlot);
PI6=getPosteriorPredictions2(exp(Theta6(:,:)),PI,simFun,PI.observablesPlot);
PI7=getPosteriorPredictions2(exp(Theta7(:,:)),PI,simFun,PI.observablesPlot);
PI8=getPosteriorPredictions2(exp(Theta8(:,:)),PI,simFun,PI.observablesPlot);
PI9=getPosteriorPredictions2(exp(Theta9(:,:)),PI,simFun,PI.observablesPlot);
PI10=getPosteriorPredictions2(exp(Theta10(:,:)),PI,simFun,PI.observablesPlot);
PI11=getPosteriorPredictions2(exp(Theta11(:,:)),PI,simFun,PI.observablesPlot);
PI12=getPosteriorPredictions2(exp(Theta12(:,:)),PI,simFun,PI.observablesPlot);
PI13=getPosteriorPredictions2(exp(Theta13(:,:)),PI,simFun,PI.observablesPlot);
PI14=getPosteriorPredictions2(exp(Theta14(:,:)),PI,simFun,PI.observablesPlot);
PI15=getPosteriorPredictions2(exp(Theta15(:,:)),PI,simFun,PI.observablesPlot);
PI16=getPosteriorPredictions2(exp(Theta16(:,:)),PI,simFun,PI.observablesPlot);

%% plot Each species
figure
nrow=ceil(sqrt(length(PI.observablesPlot)));
ncol=ceil(length(PI.observablesPlot)/nrow);
colors=linspecer(length(PI.data));
for i=1:length(PI.observablesPlot)
    output_i=PI.observablesPlot{i};
    subplot(nrow,ncol,i)
    hold on
for j=1:length(PI.data)
    plot(PI14.tspan,mean(PI14.output(j).(output_i),1,'omitnan'),'Color','black','LineStyle','--','LineWidth',3);
    plot(PI14.tspan,(PI14.output(j).(output_i)),'Color',colors(j,:));
end
set(gca,'YScale','log')
title(output_i)
end
legend([PI.data(:).Group])
%% Identify initial values 
PI1=initialValues(PI1,TV_0,Theta1);
PI2=initialValues(PI2,TV_0,Theta2);
PI3=initialValues(PI3,TV_0,Theta3);
PI4=initialValues(PI4,TV_0,Theta4);
PI5=initialValues(PI5,TV_0,Theta5);
PI6=initialValues(PI6,TV_0,Theta6);
PI7=initialValues(PI7,TV_0,Theta7);
PI8=initialValues(PI8,TV_0,Theta8);
PI9=initialValues(PI9,TV_0,Theta9);
PI10=initialValues(PI10,TV_0,Theta10);
PI11=initialValues(PI11,TV_0,Theta11);
PI12=initialValues(PI12,TV_0,Theta12);
PI13=initialValues(PI13,TV_0,Theta13);
PI14=initialValues(PI14,TV_0,Theta14);
PI15=initialValues(PI15,TV_0,Theta15);
PI16=initialValues(PI16,TV_0,Theta16);

%% Save output
for i=1:16
    num_i = num2str(i);
    save(strjoin({cd '/CIM/HuSim/CIM10/PI' num_i '.mat'},''),strjoin({'PI', num_i}, ''))
end
