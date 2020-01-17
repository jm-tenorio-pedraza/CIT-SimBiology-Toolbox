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

postSamples = [pk_postSamples cim_postSamples icb_postSamples];
pop_indx = [3 4 7:9 13:16];
var_indx = [5 6 17 18 19];
eta_indx = [1 2 10 11 12];

theta = repmat(postSamples(:,[pop_indx eta_indx]), N_indiv,1);

delta = theta(:,13)-mean(theta(:,13));
kpro_Tumor_human = log(0.0072);
theta(:,13) = kpro_Tumor_human+delta;
z = randn(N_indiv*N_pop,length(eta_indx)).*exp(postSamples(var_indx));
theta(:,eta_indx)= theta(:,eta_indx)+z;
plotBivariateMarginals_2(exp(theta),...
     'names',par([pop_indx eta_indx]),'interpreter', 'tex')
%% Adjust for human differences
Theta1 = exp(theta);
Theta2 = exp(theta);
allo_factor1 = [1 1 0.9 0.9 0];
allo_factor2 = [1 1 0.9 0.9 -1/4];

Theta1(:,[1 2 10 11 13])= Theta1(:,[1 2 10 11 13]).*((77/.022).^allo_factor1);
Theta2(:,[1 2 10 11 13])= Theta2(:,[1 2 10 11 13]).*((77/.022).^allo_factor2);

TV_range = [pi*4/3*.5^3 pi*4/3*2.5^3];
TV = (rand(size(Theta1,1),1)*(TV_range(2)-TV_range(1))+TV_range(1))/0.00153;
Theta1 = [Theta1 TV];
Theta2 = [Theta2 TV];
%% Load project
out = sbioloadproject('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/sbio projects/CIM_3.sbproj');
% Extract model
model=out.m1;
variants = get(model,'variants');
human = variants(5);
% Solver configuration
cs=model.getconfigset;
cs.StopTime = 365;
set(cs.SolverOptions, 'AbsoluteTolerance', 1.0e-12);
set(cs.SolverOptions, 'RelativeTolerance', 1.0e-12);
set(cs, 'MaximumWallClock', 0.25)

%% Parameter setup
parameters = [par([pop_indx eta_indx]) 'T_0'];
parameters(1) = {'CL_antiPDL1'};
% Define outputs% Define outputs
groups_subset = {'MOC1_Control', 'MOC1_antiPDL1','MOC1_antiCTLA4', 'MOC1_antiCTLA4_antiPDL1'};
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

%% Doses
dosing_times = 0:21:365;

antiPDL1_dose = 770*1e-3/1.5e5*1e6;
antiCTLA4_dose = 385*1e-3/1.5e5*1e6;

control = table(dosing_times', repelem(0,length(dosing_times),1),...
    repelem(0/60, length(dosing_times),1),'VariableNames',{'Time' 'Amount' 'Rate'});
antiPDL1 = table(dosing_times', repelem(antiPDL1_dose,length(dosing_times),1),...
    repelem(antiPDL1_dose/60, length(dosing_times),1),'VariableNames',{'Time' 'Amount' 'Rate'});
antiCTLA4 = table(dosing_times', repelem(antiCTLA4_dose,length(dosing_times),1),...
    repelem(antiCTLA4_dose/60, length(dosing_times),1),'VariableNames',{'Time' 'Amount' 'Rate'});

control.Properties.VariableUnits={'day' 'micromole' 'micromole/minute'};
antiPDL1.Properties.VariableUnits={'day' 'micromole' 'micromole/minute'};
antiCTLA4.Properties.VariableUnits={'day' 'micromole' 'micromole/minute'};

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

colors = linspecer(4);
colors_i = zeros(size(Theta1,1),3);
groups = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1+antiCTLA4'};
%% Plot posterior predictions (Tumor growth scaling exponent: 0]
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


%% Plot posterior predictions (Tumor growth scaling exponent: -1/4]
response2 = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1),zeros(4,1),zeros(4,1),zeros(4,1),...
    'VariableNames', {'Response' 'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'});

% Calculate SLD and response %
for i=1:4
    PI2.output(i).SLD = ((PI2.output(i).TV./(Theta1(:,end)*0.00153)).^(1/3))*100-100;
    pd = (any(PI2.output(i).SLD>20,2));
    pr = and((and(any(PI2.output(i).SLD<=-30,2), ~any(PI2.output(i).SLD<-99,2))),~pd);
    cr = and((any(PI2.output(i).SLD<-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response2{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
    colors_i(pd,:) = repmat(colors(1,:),sum(pd),1);
    colors_i(sd,:) = repmat(colors(2,:),sum(sd),1);
    colors_i(pr,:) = repmat(colors(3,:),sum(pr),1);
    colors_i(cr,:) = repmat(colors(4,:),sum(cr),1);
    PI2.output(i).Response = [pd sd pr cr];
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
figure
for i =1:4
   subplot(2,2,i)
   hold on
   for j=1:4
       mean_ij=mean(PI2.output(i).SLD(PI2.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI2.output(i).SLD(PI2.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI2.tspan/7, mean_ij,std_ij);
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
       mean_ij=mean(PI2.output(i).CD8(PI2.output(i).Response(:,j),:),'omitnan');
       std_ij = std(PI2.output(i).CD8(PI2.output(i).Response(:,j),:),'omitnan');
    h=errorbar(PI2.tspan/7, mean_ij,std_ij);
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

%% PFS
% Calculate SLD and response %
cutoff_indx = find(PI.tspan==181);
for i=1:4
    pd = any(PI1.output(i).SLD(:,1:cutoff_indx)>20,2);
    pr = and((and((PI1.output(i).SLD(:,cutoff_indx)<=-30), ~(PI1.output(i).SLD(:,cutoff_indx)<-99))),~pd);
    cr = and(((PI1.output(i).SLD(:,cutoff_indx)<-99)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response1{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
    PI1.output(i).PFS = [pd sd pr cr];
end


cutoff_indx = find(PI.tspan==181);
T = zeros(size(Theta1,1)*4,1);
for i=1:4
    rowindx = (i-1)*size(Theta1,1);
    pd = rowfun(@(x)detectPD(x,PI.tspan(1:cutoff_indx)),table(PI1.output(i).SLD(:,1:cutoff_indx)));
    T(rowindx+1:rowindx+size(Theta1,1),1)=pd{:,1};
end
zeroV = zeros(size(Theta1,1),1);
oneV = ones(size(Theta1,1),1);
X = [[oneV; zeroV; zeroV; zeroV] [zeroV; oneV; zeroV; oneV], [zeroV; zeroV; oneV; oneV]]...
    ;
[b,logl,H,stats] = coxphfit(X,T);


figure()
ax1 = gca;
censor = T==181;

[f,x] = ecdf(T(1:5e3,1),'Censoring',censor(1:5e3),'function','survivor');
stairs(x,f,'--k','LineWidth', 2)
hold on
[f,x] = ecdf(T(5e3+1:10e3,1),'Censoring',censor(5e3+1:10e3),'function','survivor');
stairs(x,f,'--r','LineWidth', 2)
[f,x] = ecdf(T(10e3+1:15e3,1),'Censoring',censor(10e3+1:15e3),'function','survivor');
stairs(x,f,'--g','LineWidth', 2)
[f,x] = ecdf(T(15e3+1:20e3,1),'Censoring',censor(15e3+1:20e3),'function','survivor');
stairs(x,f,'--b','LineWidth', 2)

legend('Control','antiPDL1','antiCTLA4', 'antiCTLA4+antiPDL1')
title('Progression-free survival (6 months)')
ylabel('S(t)')
xlabel('Time [days]')
PI1.Theta = Theta1;
PI1.parameters = parameters;
PI2.Theta = Theta2;
PI2.parameters = parameters;

%%
save(strjoin({cd '/CIM/HuSim/PI1.mat'},''),'PI1')
save(strjoin({cd '/CIM/HuSim/PI2.mat'},''),'PI2')
