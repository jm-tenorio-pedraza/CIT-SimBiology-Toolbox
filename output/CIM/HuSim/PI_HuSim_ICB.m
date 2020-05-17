%% PK general setup

% Search paths
clear all
warning off
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output')
sensitivity = false;

%% Load project 
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_4.sbproj');
% Extract model
model=out.m1;
cs=model.getconfigset;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-11);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-9);
set(cs, 'MaximumWallClock', 0.25)
MODEL = 'CIM 4';
variants = get(model, 'variants');

%% Setting up parameters, data and simulations

parameters = {'Blood'; 'Tumor'; 'Peripheral'; 'CL_antiPDL1'; 'Q23'; 'Q12';...
    'kdeg_PDL1';'kin_CD8';'K_IFNg';'KDE_MDSC';'K_MDSC'; ...
    'kin_MDSC';'kin_TIC';'kpro_Tumor'; 'kpro_Tumor_Linear';'kill_CD8'; 'S_L';...
    'S_R'; 'K_CTLA4';'K_PDL1'; 'kill_Treg'};
parameters = [parameters; 'T_0'];
% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1', 'MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
observables={'TV' 'CD8' 'Tumor.antiPDL1' 'Tumor.antiCTLA4'};
stateVar={'Tumor'};
doses = {'Blood.Dose_antiPDL1' 'Blood.Dose_antiCTLA4'};

PI=getPIData2('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/data/PI_Clavijo.mat',...
    stateVar,groups_subset,observables,'zeroAction', 'omit','mergePhenotypes',...
    false,'output', 'mean','maxIIV', true);
PI.data = PI.data([1 4 2 3]);
PI.variableUnits={'Volume [mL]'};
PI.normIndx = [];
PI.model = 'Human simulations';
PI.observablesPlot={'TV' 'CD8' 'Tumor_antiPDL1' 'Tumor_antiCTLA4'};
PI.tspan = 1:3:365;

% Get simulation function
[sim,PI.u]=initializePI(model,parameters,observables,PI,doses, 'Human','doseUnits', 'mole');
% Get initial values
%% Doses (Durvalumab + Tremelimumab)
n_doses_antiPDL1_mono = 27;
n_doses_antiCTLA4_mono = 9;
dosing_times_antiPDL1_mono = 0:14:365;
dosing_times_antiCTLA4_mono = [0:28:28*6 (28*6+12*7):12*7:365];
dosing_times_antiPDL1_combi = [0:28:28*3 (28*3+7*2):14:365];
dosing_times_antiCTLA4_combi = 0:28:28*3;

maxDosingFreq = dosing_times_antiPDL1_mono; % Maximum dosing frequency
nDoses_max = length(maxDosingFreq);

doseScalingFactor = 77*1e-3/1.5e5*1e6; % kg*mg/kg*g/mg*mole/g*µmol/mol

control = table(maxDosingFreq', repelem(0,nDoses_max,1),...
    repelem(0/60, nDoses_max,1),'VariableNames',{'Time' 'Amount' 'Rate'});

antiPDL1_mono = table(dosing_times_antiPDL1_mono', repelem(10*doseScalingFactor,n_doses_antiPDL1_mono,1),...
    repelem(10*doseScalingFactor/60, n_doses_antiPDL1_mono,1),'VariableNames',{'Time' 'Amount' 'Rate'});

antiCTLA4_mono = table(dosing_times_antiCTLA4_mono',...
    repelem(10*doseScalingFactor,n_doses_antiCTLA4_mono,1),...
    repelem(10*doseScalingFactor/60, n_doses_antiCTLA4_mono,1),...
    'VariableNames',{'Time' 'Amount' 'Rate'});

antiPDL1_combi = table(dosing_times_antiPDL1_combi',...
    [repelem(20*doseScalingFactor, 4,1); repelem(10*doseScalingFactor, length(dosing_times_antiPDL1_combi)-4,1)],...
    [repelem(20*doseScalingFactor/60, 4,1); repelem(10*doseScalingFactor/60, length(dosing_times_antiPDL1_combi)-4,1)],...
    'VariableNames',{'Time' 'Amount' 'Rate'});

antiCTLA4_combi = table(dosing_times_antiCTLA4_combi',...
    repelem(10*doseScalingFactor,4,1), repelem(10*doseScalingFactor/60,4,1),...
    'VariableNames',{'Time' 'Amount' 'Rate'});

control.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiPDL1_mono.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiCTLA4_mono.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiPDL1_combi.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiCTLA4_combi.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};

PI.u = {};
PI.u(1,1:2) = {control control};
PI.u(2,1:2) = {antiPDL1_mono control};
PI.u(3, 1:2) = {control antiCTLA4_mono};
PI.u(4, 1:2) = {antiPDL1_combi antiCTLA4_combi};
clearvars n_doses_antiPDL1_mono n_doses_antiCTLA4_mono dosing_times_antiPDL1_mono...
    dosing_times_antiCTLA4_mono dosing_times_antiPDL1_combi maxDosingFreq doseScalingFactor ...
    control antiPDL1_mono antiCTLA4_mono antiPDL1_combi antiPDL1_combi
%% Hierarchical model simulation
PI.H = getHierarchicalStruct(parameters(1:end-1),PI,'n_sigma', length(observables),...
    'rand_indx', [14:16],'cell_indx',[2 4 12 13], 'n_indiv', length(PI.u),'CellField', 'Name');
% Generating PI
SigmaNames = getVarNames(PI, observables);
[beta, sigma_prior] = getVarValues([0.1 0.1 .001], [0.1, 0.1 0.001], [1 1 1], PI);
lb = [1e-3  1e-3    1e-4    1e-4    1e-4    1e0     1e0    1e-3   1e-3    1e-4 1e-3   1e-4    1e-3    1e-2   1e-2     1e-3 1e-7     1e-7   1e0   1e3    1e-4];
ub = [1e1   1e1     1e2     1e2     1e1     1e6    1e6   1e2     1e2     1e1 1e1    1e2     1e2     1e1    1e2  1e3 1        1      1e4   1e6     1e3];
PI.par = getParamStruct2(sim,PI.H,size(PI.data,1),beta,...
    SigmaNames,'Sigma', sigma_prior,'LB', lb', 'UB', ub');
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
PI_PK = load(strjoin({cd 'PK_mAb_ThreeComp/PI/ThreeComp_4/PI_PK_ThreeComp4_4_TMDD_11.mat'},'/'));
PI_CIM = load(strjoin({cd 'CIM/PI/CIM10/PI_CIM21_Control_14.mat'},'/'));
PI_ICB = load(strjoin({cd 'CIM/PI/CIM8/PI_CIM4_ICB_21_1.mat'},'/'));
Meta=[];
Meta(1).Struct = PI_PK;
Meta(2).Struct = PI_CIM;
Meta(3).Struct = PI_ICB;
%% Select which dimensions to sample from 
N_pop = 500;
N_indiv = 10;
[theta,H] = getTheta(Meta, parameters, N_pop, N_indiv);
Theta1 = theta;

kpro_Tumor_human = log(0.0072);
delta = theta(:,14)-mean(theta(:,14));
Theta1(:,14) = kpro_Tumor_human+delta;

TV_range = [pi*4/3*.5^3 pi*4/3*2.5^3];
TV = log((rand(size(theta,1),1)*(TV_range(2)-TV_range(1))+TV_range(1))/0.00153);
Theta1(:,end) = TV;

scalingExp1 = [0.8 0.8 0.8 0.8 0.8 0.8 0 0 0 0 0 0 0 0 0 -1/4 0 0 0 0 0 0];
scalingExp2 = [0.8 0.8 0.8 0.8 0.8 0.8 0 0 0 0 0 0 0 0 0 -1/4 0 0 0 0 0 0];

scalingFactor1 = (77/.022).^(scalingExp1);

Theta1 = Theta1+log(scalingFactor1);
Theta1 = exp(Theta1);
table(parameters, exp(mean(log(Theta1))'))
sigma_indx = randsample(size(PI_CIM.PI.postSamples,1),N_pop*N_indiv,true);
sigma = nan(N_pop*N_indiv,2);
sigma(:,2) = PI_CIM.PI.postSamples(sigma_indx,end-6);
sigma_indx = randsample(size(PI_ICB.PI.postSamples,1),N_pop*N_indiv,true);
sigma(:,1) = PI_ICB.PI.postSamples(sigma_indx,end);
sigma = exp(sigma);
%% Plot bivariate marginals
plotBivariateMarginals_2((Theta1(:,:)),...
       'names',parameters,'interpreter', 'tex')
plotCorrMat(Theta1, parameters)

%% Posterior predictions
simFun=@(x)getOutput2(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    {},PI.normIndx,'Output', 'data');
dataOutput = simFun((Theta1(1,:)));
[PI.data(1:end).simValue] = dataOutput{:,:};
figure
hold on
arrayfun(@(x)plot(PI.tspan, x.simValue(:,1), 'Color', x.colors), PI.data, 'UniformOutput', false)

tic
PI1=getPosteriorPredictions2((Theta1(:,:)),PI,simFun,PI.observablesPlot);
toc

%% Plot posterior predictions (Killing rate scaling exponent: -1/4]
colors = linspecer(4);
colors_i = zeros(size(Theta1,1),3);
groups = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'};
response2 = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1),zeros(4,1),zeros(4,1),zeros(4,1),...
    'VariableNames', {'Response' 'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'});
lineStyle = {'-' '-.' '--' ':'};
response2_pred = response2;
% Calculate SLD and response %
for i=1:4
    PI1.output(i).SLD = ((PI1.output(i).TV./(Theta1(:,end)*0.00153)).^(1/3))*100-100;
    PI1.output(i).TV_sigma = exp(log(PI1.output(i).TV)+randn(N_indiv*N_pop,...
        size(PI1.output(i).TV,2)).*sigma(:,1));
    PI1.output(i).CD8_sigma = exp(log(PI1.output(i).CD8)+randn(N_indiv*N_pop,...
        size(PI1.output(i).CD8,2)).*sigma(:,2));
    PI1.output(i).SLD_sigma = ((PI1.output(i).TV_sigma./(Theta1(:,end)*0.00153)).^(1/3))*100-100;
    
    pd = and(any(PI1.output(i).SLD>20,2),...
        any(PI1.output(i).TV/(4/3*pi).^(1/3)...
        -(PI1.output(i).TV(:,1)/(4/3*pi)).^(1/3)>0.5,2));
    
    pr = and((and(any(PI1.output(i).SLD<=-30,2), ~any(PI1.output(i).SLD<-99,2))),~pd);
    cr = and((any(PI1.output(i).SLD<-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response2{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];

    pd_sigma = (any(PI1.output(i).SLD_sigma>20,2));
    pr_sigma = and((and(any(PI1.output(i).SLD_sigma<=-30,2), ~any(PI1.output(i).SLD_sigma<-99,2))),~pd);
    cr_sigma = and((any(PI1.output(i).SLD_sigma<-99,2)),~pd);
    sd_sigma = and(and(~pd_sigma,~pr_sigma), ~cr_sigma);
    response2_pred{1:4,i+1} = [mean(pd_sigma); mean(sd_sigma); mean(pr_sigma); mean(cr_sigma)];
    
    colors_i(pd,:) = repmat(colors(1,:),sum(pd),1);
    colors_i(sd,:) = repmat(colors(2,:),sum(sd),1);
    colors_i(pr,:) = repmat(colors(3,:),sum(pr),1);
    colors_i(cr,:) = repmat(colors(4,:),sum(cr),1);
    PI1.output(i).Response = [pd sd pr cr];
    PI1.output(i).Response_sigma = [pd_sigma sd_sigma pr_sigma cr_sigma];

    PI1.output(i).colors = colors_i;
end
% Plot SLD
indx = zeros(4,1);
figure
for i=1:4
subplot(2,2,i)
h=plot(PI1.tspan/7, PI1.output(i).SLD);

for j=1:length(h)
    h(j).Color = PI1.output(i).colors(j,:);
end
hold on
plot(PI1.tspan/7, repelem(20,1,length(PI1.tspan)), '--k')
plot(PI1.tspan/7, repelem(-30,1,length(PI1.tspan)), '--k')
for k=1:4
indx(k) = find(PI1.output(i).Response(:,k),1);
end
legend(h(indx),{'PD'; 'SD'; 'PR'; 'CR'})
set(gca, 'YSCale', 'linear','YLim', [-100 100])
xlabel('Time [weeks]')
ylabel('Change in tumor diameter wrt baseline')
title (strjoin({'Tumor response in virtual patients (' groups{i} ')'},''))
end
% Plot mean SLD
figure('Renderer', 'painters', 'Position', [10 10 800 700])
for i =1:4
   subplot(2,2,i)
   hold on
   for j=1:4
       mean_ij=mean(PI1.output(i).SLD(PI1.output(i).Response(:,j),:),'omitnan');
       ci_ub = quantile(PI1.output(i).SLD(PI1.output(i).Response(:,j),:),.975);
       ci_lb = quantile(PI1.output(i).SLD(PI1.output(i).Response(:,j),:),.025);
       
       pi_ub = quantile(PI1.output(i).SLD_sigma(PI1.output(i).Response(:,j),:),.975);
       pi_lb = quantile(PI1.output(i).SLD_sigma(PI1.output(i).Response(:,j),:),.025);
       
    h=plot(PI1.tspan/7, mean_ij);
    ci_plot = patch('XData', [PI1.tspan/7 PI1.tspan(end:-1:1)/7], 'YData', [ci_ub ci_lb(end:-1:1)]);
    pi_plot = patch('XData', [PI1.tspan/7 PI1.tspan(end:-1:1)/7], 'YData', [pi_ub pi_lb(end:-1:1)]);
    h.Color = colors(j,:);
    h.LineWidth = 2;
    ci_plot.FaceColor = colors(j,:);
    ci_plot.LineStyle = lineStyle{j};
    ci_plot.EdgeColor = colors(j,:);
    ci_plot.FaceAlpha = 0.2;
    
    pi_plot.LineStyle = 'none';
    pi_plot.FaceColor = colors(j,:);
    pi_plot.FaceAlpha = 0.01;
   end
   ax = gca;
   legend(ax.Children(end:-3:1),{'PD'; 'SD'; 'PR'; 'CR'})

 set(gca, 'YSCale', 'linear','YLim', [-100 100])
 xlabel('Time [weeks]')
ylabel('Change in tumor diameter wrt baseline')
title (strjoin({'Mean tumor response in virtual patients (' groups{i} ')'},''),'interpreter', 'none')

end

% Mean CD8
figure('Renderer', 'painters', 'Position', [10 10 800 700])
for i=1:4
subplot(2,2,i)
hold on
 for j=1:4
        mean_ij=mean(PI1.output(i).CD8(PI1.output(i).Response(:,j),:),'omitnan');
       ci_ub = quantile(PI1.output(i).CD8(PI1.output(i).Response(:,j),:),.975);
       ci_lb = quantile(PI1.output(i).CD8(PI1.output(i).Response(:,j),:),.025);
       
       pi_ub = quantile(PI1.output(i).CD8_sigma(PI1.output(i).Response(:,j),:),.975);
       pi_lb = quantile(PI1.output(i).CD8_sigma(PI1.output(i).Response(:,j),:),.025);
       
    h=plot(PI1.tspan/7, mean_ij);
    ci_plot = patch('XData', [PI1.tspan/7 PI1.tspan(end:-1:1)/7], 'YData', [ci_ub ci_lb(end:-1:1)]);
    pi_plot = patch('XData', [PI1.tspan/7 PI1.tspan(end:-1:1)/7], 'YData', [pi_ub pi_lb(end:-1:1)]);
    h.Color = colors(j,:);
    h.LineWidth = 2;
    ci_plot.FaceColor = colors(j,:);
    ci_plot.LineStyle = lineStyle{j};
    ci_plot.EdgeColor = colors(j,:);
    ci_plot.FaceAlpha = 0.2;
    
    pi_plot.LineStyle = 'none';
    pi_plot.FaceColor = colors(j,:);
    pi_plot.FaceAlpha = 0.01;
end
     ax = gca;
   legend(ax.Children(end:-3:1),{'PD'; 'SD'; 'PR'; 'CR'},'FontSize', 10)
xlabel('Time [weeks]')
ylabel('Percentage [%]')
title (strjoin({'CD8+ T-cell response in virtual patients (' groups{i} ')'},''))
ylim([0 5])
end

% Mean antiPDL1
figure
for i=1:4
subplot(2,2,i)
hold on
 for j=1:4
       mean_ij=mean(PI1.output(i).Tumor_antiPDL1(PI1.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI1.output(i).Tumor_antiPDL1(PI1.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI1.tspan/7, mean_ij,std_ij);
    h.Color = colors(j,:);
 end
   legend({'PD'; 'SD'; 'PR'; 'CR'})

xlabel('Time [weeks]')
ylabel('Concentration [mg/L]')
% ylim([0 0.5])
title (strjoin({'anti-PD-L1 concentrations in tumor tissue of virtual patients (' groups{i} ')'},''))

end

% Mean antiCTLA4
figure
for i=1:4
subplot(2,2,i)
hold on
 for j=1:4
       mean_ij=mean(PI1.output(i).Tumor_antiCTLA4(PI1.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI1.output(i).Tumor_antiCTLA4(PI1.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI1.tspan/7, mean_ij,std_ij);
    h.Color = colors(j,:);
 end
   legend({'PD'; 'SD'; 'PR'; 'CR'})

xlabel('Time [weeks]')
ylabel('Concentration [mg/L]')
ylim([0 0.2])
title (strjoin({'anti-PD-L1 concentrations in tumor tissue of virtual patients (' groups{i} ')'},''))

end

%% PFS in theta1
% Calculate SLD and response %
[PI1, response1] = getPFS(PI1, groups, 181);
[PI1,T, censor] = getSurvivalTime(PI1, 181, groups, Theta1);
plotSurvival(T,censor,Theta1,groups)
% subplot(1,2,1)
plotSurvivalFunction(PI1,181,groups)
title('Progresson-free survival at 6 months (kill_{CD8} \propto M^{-0.75})')

%% Analysing parameter-output relations
corrInOut = plotInputToOutput(Theta1, {'SLD'}, PI1, groups, parameters);
corrInOut2 = plotInputToOutput(PI1.Theta, {'SLD'}, PI1, groups, parameters,'plotOutput', true);

