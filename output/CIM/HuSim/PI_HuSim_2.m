%% File paths
warning on
clear all
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output')

%% Load PIs and posterior parameter samples
PI_PK = load(strjoin({cd 'PK_mAb_TwoComp/PI/PI_PK_CE.mat'},'/'));
PI_CIM = load(strjoin({cd 'CIM/PI/CIM3/PI_CIM_Control_3_red.mat'},'/'));
PI_ICB = load(strjoin({cd 'CIM/PI/CIM3/PI_CIM_ICB_1.mat'},'/'));

%% Select which dimensions to sample from 
N_pop = 500;
N_indiv = 10;
pk_col_indx = [PI_PK.PI.H.PopulationParams PI_PK.PI.H.IndividualParams.OmegaIndex];
cim_col_indx = [3 5 6];
icb_col_indx = [PI_ICB.PI.H.PopulationParams PI_ICB.PI.H.CellParams.OmegaIndex ...
    PI_ICB.PI.H.IndividualParams.OmegaIndex PI_ICB.PI.H.SigmaParams(end)];

par_PK = {PI_PK.PI.par(pk_col_indx).name};
par_CIM = {PI_CIM.PI.par(cim_col_indx).name};
par_ICB = {PI_ICB.PI.par(icb_col_indx).name};

pk_row_indx = randsample(size(PI_PK.PI.postSamples,1),N_pop, true);
cim_row_indx = randsample(size(PI_CIM.PI.postSamples,1), N_pop, true);
icb_row_indx = randsample(size(PI_ICB.PI.postSamples,1), N_pop, true);

par = [par_PK par_CIM par_ICB];
%% Generate sample (Theta)
pk_postSamples = PI_PK.PI.postSamples(pk_row_indx,pk_col_indx);
cim_postSamples = PI_CIM.PI.postSamples(cim_row_indx,cim_col_indx);
icb_postSamples = PI_ICB.PI.postSamples(icb_row_indx,icb_col_indx);

sigmaSamples = [PI_ICB.PI.postSamples(cim_row_indx, PI_ICB.PI.H.SigmaParams(4)) PI_CIM.PI.postSamples(cim_row_indx, PI_CIM.PI.H.SigmaParams(4))];
postSamples = [pk_postSamples cim_postSamples icb_postSamples];
pop_indx = [3 4 7:9 13:16];
var_indx = [5 6 17 18 19];
eta_indx = [1 2 10 11 12];

theta = repmat(postSamples(:,[pop_indx eta_indx]), N_indiv,1);
sigma = exp(repmat(sigmaSamples,N_indiv,1));

delta = theta(:,13)-mean(theta(:,13));
kpro_Tumor_human = log(0.0072);
theta(:,13) = kpro_Tumor_human+delta;
z = randn(N_indiv*N_pop,length(eta_indx)).*exp(postSamples(var_indx));
theta(:,eta_indx)= theta(:,eta_indx)+z;
% plotBivariateMarginals_2(exp(theta),...
%      'names',par([pop_indx eta_indx]),'interpreter', 'tex')
%% Adjust for human differences
Theta1 = exp(theta);
Theta2 = exp(theta);
allo_factor1 = [1 1 0.9 0.9 0];
allo_factor2 = [1 1 0.9 0.9 -1/4];

Theta1(:,[1 2 10 11 12])= Theta1(:,[1 2 10 11 12]).*((77/.022).^allo_factor1);
Theta2(:,[1 2 10 11 12])= Theta2(:,[1 2 10 11 12]).*((77/.022).^allo_factor2);

TV_range = [pi*4/3*.5^3 pi*4/3*2.5^3];
TV = (rand(size(Theta1,1),1)*(TV_range(2)-TV_range(1))+TV_range(1))/0.00153;
Theta1 = [Theta1 TV];
Theta2 = [Theta2 TV];
%% Load project
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_3.sbproj');
% Extract model
model=out.m1;
variants = get(model,'variants');
human = variants(2);
% Solver configuration
cs=model.getconfigset;
cs.StopTime = 365;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-12);
set(cs, 'MaximumWallClock', 0.25)

%% Parameter setup
parameters = [par([pop_indx eta_indx]) 'T_0'];
parameters(1) = {'CL_antiPDL1'};
% plotBivariateMarginals_2((Theta1),...
%      'names',parameters,'interpreter', 'none')

% Define outputs% Define outputs
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
PI.H.PopulationParams = 1:15;
PI.H.SigmaParams = 16;
PI.H.IndividualParams.OmegaIndex = [];
PI.H.CellParams.OmegaIndex = [];
% Get simulation function
sim = createSimFunction(model,parameters,observables, doses,human);

%% Doses (Atezolizumab + Ipilimumab)
n_doses_antiPDL1 = 18;
n_doses_antiCTLA4 = 4;
dosing_times_antiPDL1 = 0:21:21*(n_doses_antiPDL1-1);
dosing_times_antiCTLA4 = 0:21:21*(n_doses_antiCTLA4-1);

antiPDL1_dose = 77*10*1e-3/1.5e5*1e6; % kg*mg/kg*g/mg*mole/g*µmol/mol
antiCTLA4_dose = 77*3*1e-3/1.5e5*1e6;   % kg*mg/kg*g/mg*mole/g*µmol/mol

control = table(dosing_times_antiPDL1', repelem(0,n_doses_antiPDL1,1),...
    repelem(0/60, n_doses_antiPDL1,1),'VariableNames',{'Time' 'Amount' 'Rate'});
antiPDL1 = table(dosing_times_antiPDL1', repelem(antiPDL1_dose,n_doses_antiPDL1,1),...
    repelem(antiPDL1_dose/60, n_doses_antiPDL1,1),'VariableNames',{'Time' 'Amount' 'Rate'});
antiCTLA4 = table(dosing_times_antiPDL1', [repelem(antiCTLA4_dose,n_doses_antiCTLA4,1);
    repelem(0,n_doses_antiPDL1-n_doses_antiCTLA4,1)],...
    [repelem(antiCTLA4_dose/60, n_doses_antiCTLA4,1); repelem(0,n_doses_antiPDL1-n_doses_antiCTLA4,1)]...
    ,'VariableNames',{'Time' 'Amount' 'Rate'});

control.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiPDL1.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};
antiCTLA4.Properties.VariableUnits = {'day' 'micromole' 'micromole/minute'};

PI.u = {};
PI.u(1,1:2) = {control control};
PI.u(2,1:2) = {antiPDL1 control};
PI.u(3, 1:2) = {control antiCTLA4};
PI.u(4, 1:2) = {antiPDL1 antiCTLA4};

%% Posterior predictions
simFun=@(x)getOutput(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    {},PI.normIndx, PI.H);
PI = simFun(Theta1(5,:));
tic
PI1=getPosteriorPredictions(Theta1,PI,simFun,PI.observablesPlot);
toc

tic
PI2=getPosteriorPredictions(Theta2,PI,simFun,PI.observablesPlot);
toc

%% Plot posterior predictions (Tumor growth scaling exponent: 0]
colors = linspecer(4);
colors_i = zeros(size(Theta1,1),3);
groups = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'};

response1 = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1),zeros(4,1),zeros(4,1),zeros(4,1),...
    'VariableNames', {'Response' 'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'});

% Calculate SLD and response %
for i=1:4
    PI1.output(i).SLD = ((PI1.output(i).TV./(Theta1(:,end)*0.00153)).^(1/3))*100-100;
    pd = (any(PI1.output(i).SLD>20,2));
    pr = and((and(any(PI1.output(i).SLD<=-30,2), ~any(PI1.output(i).SLD<-99,2))),~pd);
    cr = and((any(PI1.output(i).SLD<-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response1{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
    colors_i(pd,:) = repmat(colors(1,:),sum(pd),1);
    colors_i(sd,:) = repmat(colors(2,:),sum(sd),1);
    colors_i(pr,:) = repmat(colors(3,:),sum(pr),1);
    colors_i(cr,:) = repmat(colors(4,:),sum(cr),1);
    PI1.output(i).Response = [pd sd pr cr];
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
figure
for i =1:4
   subplot(2,2,i)
   hold on
   for j=1:4
       mean_ij=mean(PI1.output(i).SLD(PI1.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI1.output(i).SLD(PI1.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI1.tspan/7, mean_ij,std_ij);
    h.Color = colors(j,:);
   end
   legend({'PD'; 'SD'; 'PR'; 'CR'})

 %set(gca, 'YSCale', 'linear','YLim', [-100 100])
 xlabel('Time [weeks]')
ylabel('Change in tumor diameter wrt baseline')
title (strjoin({'Mean tumor response in virtual patients (' groups{i} ')'},''))

end

% Mean CD8
figure
for i=1:4
subplot(2,2,i)
hold on
 for j=1:4
       mean_ij=mean(PI1.output(i).CD8(PI1.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI1.output(i).CD8(PI1.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI1.tspan/7, mean_ij,std_ij);
    h.Color = colors(j,:);
 end
   legend({'PD'; 'SD'; 'PR'; 'CR'})

xlabel('Time [weeks]')
ylabel('Percentage [%]')
title (strjoin({'CD8+ T-cell response in virtual patients (' groups{i} ')'},''))

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
ylim([0 0.5])
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


%% Plot posterior predictions (Killing rate scaling exponent: -1/4]
response2 = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1),zeros(4,1),zeros(4,1),zeros(4,1),...
    'VariableNames', {'Response' 'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'});
lineStyle = {'-' '-.' '--' ':'};
response2_pred = response2;
% Calculate SLD and response %
for i=1:4
    PI2.output(i).SLD = ((PI2.output(i).TV./(Theta1(:,end)*0.00153)).^(1/3))*100-100;
    PI2.output(i).TV_sigma = exp(log(PI2.output(i).TV)+randn(N_indiv*N_pop,...
        size(PI2.output(i).TV,2)).*sigma(:,1));
    PI2.output(i).CD8_sigma = exp(log(PI2.output(i).CD8)+randn(N_indiv*N_pop,...
        size(PI2.output(i).CD8,2)).*sigma(:,2));
    PI2.output(i).SLD_sigma = ((PI2.output(i).TV_sigma./(Theta1(:,end)*0.00153)).^(1/3))*100-100;
    
    pd = and(any(PI2.output(i).SLD>20,2),...
        any(PI2.output(i).TV/(4/3*pi).^(1/3)...
        -(PI2.output(i).TV(:,1)/(4/3*pi)).^(1/3)>0.5,2));
    
    pr = and((and(any(PI2.output(i).SLD<=-30,2), ~any(PI2.output(i).SLD<-99,2))),~pd);
    cr = and((any(PI2.output(i).SLD<-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response2{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];

    pd_sigma = (any(PI2.output(i).SLD_sigma>20,2));
    pr_sigma = and((and(any(PI2.output(i).SLD_sigma<=-30,2), ~any(PI2.output(i).SLD_sigma<-99,2))),~pd);
    cr_sigma = and((any(PI2.output(i).SLD_sigma<-99,2)),~pd);
    sd_sigma = and(and(~pd_sigma,~pr_sigma), ~cr_sigma);
    response2_pred{1:4,i+1} = [mean(pd_sigma); mean(sd_sigma); mean(pr_sigma); mean(cr_sigma)];
    
    colors_i(pd,:) = repmat(colors(1,:),sum(pd),1);
    colors_i(sd,:) = repmat(colors(2,:),sum(sd),1);
    colors_i(pr,:) = repmat(colors(3,:),sum(pr),1);
    colors_i(cr,:) = repmat(colors(4,:),sum(cr),1);
    PI2.output(i).Response = [pd sd pr cr];
    PI2.output(i).Response_sigma = [pd_sigma sd_sigma pr_sigma cr_sigma];

    PI2.output(i).colors = colors_i;
end
% Plot SLD
indx = zeros(4,1);
figure
for i=1:4
subplot(2,2,i)
h=plot(PI2.tspan/7, PI2.output(i).SLD);

for j=1:length(h)
    h(j).Color = PI2.output(i).colors(j,:);
end
hold on
plot(PI2.tspan/7, repelem(20,1,length(PI2.tspan)), '--k')
plot(PI2.tspan/7, repelem(-30,1,length(PI2.tspan)), '--k')
for k=1:4
indx(k) = find(PI2.output(i).Response(:,k),1);
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
       mean_ij=mean(PI2.output(i).SLD(PI2.output(i).Response(:,j),:),'omitnan');
       ci_ub = quantile(PI2.output(i).SLD(PI2.output(i).Response(:,j),:),.975);
       ci_lb = quantile(PI2.output(i).SLD(PI2.output(i).Response(:,j),:),.025);
       
       pi_ub = quantile(PI2.output(i).SLD_sigma(PI2.output(i).Response(:,j),:),.975);
       pi_lb = quantile(PI2.output(i).SLD_sigma(PI2.output(i).Response(:,j),:),.025);
       
    h=plot(PI2.tspan/7, mean_ij);
    ci_plot = patch('XData', [PI2.tspan/7 PI2.tspan(end:-1:1)/7], 'YData', [ci_ub ci_lb(end:-1:1)]);
    pi_plot = patch('XData', [PI2.tspan/7 PI2.tspan(end:-1:1)/7], 'YData', [pi_ub pi_lb(end:-1:1)]);
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
        mean_ij=mean(PI2.output(i).CD8(PI2.output(i).Response(:,j),:),'omitnan');
       ci_ub = quantile(PI2.output(i).CD8(PI2.output(i).Response(:,j),:),.975);
       ci_lb = quantile(PI2.output(i).CD8(PI2.output(i).Response(:,j),:),.025);
       
       pi_ub = quantile(PI2.output(i).CD8_sigma(PI2.output(i).Response(:,j),:),.975);
       pi_lb = quantile(PI2.output(i).CD8_sigma(PI2.output(i).Response(:,j),:),.025);
       
    h=plot(PI2.tspan/7, mean_ij);
    ci_plot = patch('XData', [PI2.tspan/7 PI2.tspan(end:-1:1)/7], 'YData', [ci_ub ci_lb(end:-1:1)]);
    pi_plot = patch('XData', [PI2.tspan/7 PI2.tspan(end:-1:1)/7], 'YData', [pi_ub pi_lb(end:-1:1)]);
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
       mean_ij=mean(PI2.output(i).Tumor_antiPDL1(PI2.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI2.output(i).Tumor_antiPDL1(PI2.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI2.tspan/7, mean_ij,std_ij);
    h.Color = colors(j,:);
 end
   legend({'PD'; 'SD'; 'PR'; 'CR'})

xlabel('Time [weeks]')
ylabel('Concentration [mg/L]')
ylim([0 0.5])
title (strjoin({'anti-PD-L1 concentrations in tumor tissue of virtual patients (' groups{i} ')'},''))

end

% Mean antiCTLA4
figure
for i=1:4
subplot(2,2,i)
hold on
 for j=1:4
       mean_ij=mean(PI2.output(i).Tumor_antiCTLA4(PI2.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI2.output(i).Tumor_antiCTLA4(PI2.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI2.tspan/7, mean_ij,std_ij);
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
subplot(1,2,1)
plotSurvivalFunction(PI1,181,groups)
title('Progresson-free survival at 6 months (kill_{CD8} \propto M^{0})')

[PI2, response2] = getPFS(PI2, groups, 181);
[PI2,T2, censor2] = getSurvivalTime(PI2, 181, groups, Theta2);
plotSurvival(T2,censor2,Theta2,groups)
subplot(1,2,2)
plotSurvivalFunction(PI2,181,groups)
title('Progresson-free survival at 6 months (kill_{CD8} \propto M^{-0.25})', 'FontSize', 14)

%% Analysing parameter-output relations
corrInOut = plotInputToOutput(Theta1, {'SLD'}, PI1, groups, parameters);
corrInOut2 = plotInputToOutput(PI2.Theta, {'SLD'}, PI2, groups, parameters,'plotOutput', true);

%%
PI1.Theta = Theta1;
PI1.parameters = parameters;
PI2.Theta = Theta2;
PI2.parameters = parameters;


save(strjoin({cd '/CIM/HuSim/PI1.mat'},''),'PI1')
save(strjoin({cd '/CIM/HuSim/PI2.mat'},''),'PI2')

load(strjoin({cd '/CIM/HuSim/PI2.mat'},''))
