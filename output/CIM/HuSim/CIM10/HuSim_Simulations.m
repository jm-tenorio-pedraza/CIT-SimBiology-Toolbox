%% Load PIs and posterior parameter samples
PI_MOC1 = load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_MOC1_PKPD_1_8_kinCD8.mat'},'/'));

Meta=[];
Meta(1).Struct = PI_MOC1;
values=repelem({sim.Parameters.Value},length(Meta),1);
[Meta(1:end).Values]=values{:,:};

%% Extendend theta
HumanPI=[];
HumanPI(1).Theta=Theta1;
HumanPI(2).Theta=Theta2;
HumanPI(3).Theta=Theta3;
HumanPI(4).Theta=Theta4;
HumanPI(5).Theta=Theta5;
HumanPI(6).Theta=Theta6;
HumanPI(7).Theta=Theta7;
HumanPI(8).Theta=Theta8;
% % HumanPI(9).Theta=Theta9;
% % HumanPI(10).Theta=Theta10;
% HumanPI(11).Theta=Theta11;
% HumanPI(12).Theta=Theta12;
% HumanPI(13).Theta=Theta13;
% HumanPI(14).Theta=Theta14;
% HumanPI(15).Theta=Theta15;
% HumanPI(16).Theta=Theta16;
%% Posterior predictions
simFun=@(x)getOutput2(PI,@(p)sim(p,PI.tspan(end),PI.u,PI.tspan),x,...
    {},PI.normIndx,'Output', 'data');
dataOutput = simFun(exp((HumanPI(1).Theta(5,:))));
[PI.data(1:end).simValue] = dataOutput{:,:};
plotSimValue(PI,'scale','linear')    
for i =3:length(HumanPI)
    HumanPI(i).PI=getPosteriorPredictions2(exp(HumanPI(i).Theta(1:end,:)),PI,simFun,PI.observablesPlot);
end
%% Check simulations for NaNs
checkPI = arrayfun(@(x) checkSimulations(x.PI,x.Theta),HumanPI,'UniformOutput',false);
[HumanPI(1:end).PI]=checkPI{:,:};
theta=cellfun(@(x)x.Theta,checkPI,'UniformOutput',false);
[HumanPI(1:end).Theta]=theta{:,:};

%% Sigma
sigma_indx= ismember({Meta(1).Struct.PI.par(:).name},'sigma_Tumor');
sigma = arrayfun(@(x)exp(randsample(Meta(1).Struct.PI.postSamples(:,sigma_indx), ...
    size(x.Theta,1),true,ones(size(Meta(1).Struct.PI.postSamples,1),1))),HumanPI,'UniformOutput',false);
[HumanPI(1:end).Sigma]=sigma{:,:};
%%
for i=1:length(HumanPI)
   plotSimulations(HumanPI(i).PI,'YScale','log','YLabel', PI.variableUnits)
end
%% Obtain SLD change and PFS
treatments = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'};
cutoff=24*30;
% Calculate SLD and response %
for i =1:length(HumanPI)
    HumanPI(i).PI=getSLD(HumanPI(i).PI, HumanPI(i).Sigma, 'TV_0',...
        HumanPI(i).PI.output(1).Tumor(1:end,1));
    [HumanPI(i).PI, HumanPI(i).ORR] = getORR(HumanPI(i).PI, treatments,'TV_0',...
        HumanPI(i).PI.output(1).Tumor(1:end,1),'N',size( HumanPI(i).Theta,1),'cutoff_value',12*30);
    [HumanPI(i).PI, HumanPI(i).Response] = getPFS(HumanPI(i).PI, treatments,'cutoff_value',cutoff);
    HumanPI(i).PI = getSurvivalTime(HumanPI(i).PI,...
        treatments,'N',size(HumanPI(i).Theta,1),'cutoff_value',cutoff);
    HumanPI(i).PI = getSurvivalTime(HumanPI(i).PI,...
        treatments,'N',size(HumanPI(i).Theta,1),'tumorField','Tumor','survivalType','OS','CC',...
        75);
end
%% Plot SLD and PFS
for i = 1:length(HumanPI)
    plotORR(HumanPI(i).PI, treatments,'output', {'SLD'})
end

for i = 1:length(HumanPI)
%     figure
    plotSurvivalFunction(HumanPI(i).PI,30*12,treatments)
end
%% Analysing parameter-output relations
corrInOut = arrayfun(@(x)plotInputToOutput(x.Theta, {'SLD'}, x.PI, treatments, parameters,'plotOutput',true),HumanPI,'UniformOutput',false);
corrInOut2 = plotInputToOutput(HumanPI(1).Theta, {'SLD'},HumanPI(1).PI, treatments, parameters,'plotOutput', true);
figure
hold on
arrayfun(@(x)plot((x.CD8_E(:,end)),(x.Tumor(:,end)), '+'), HumanPI(1).PI.output)
%% CD8 vs tumor
HumPI = HumanPI(1);
scatter3(exp(HumPI.Theta(:,ismember(parameters,'kill_CD8'))),...
    exp(HumPI.Theta(:,ismember(parameters,'kin_CD8'))),...
    (HumPI.PI.output(2).Tumor(:,end)),20,HumPI.PI.output(2).colors...
   )
xlabel('kill_{CD8} [1/day]') ;ylabel('kin_{CD8}[1e6 cells/day]');zlabel('Final tumor volume  [mL]');
title('Association between T-cell parameters and final tumor volumes')
grid on

scatter3(exp(HumPI.Theta(:,ismember(parameters,'T_0'))),...
    exp(HumPI.Theta(:,ismember(parameters,'kin_CD8'))),...
    (HumPI.PI.output(2).Tumor(:,end)),20,HumPI.PI.output(2).colors...
   )
xlabel('T_0 [10^6 cells]') ;ylabel('kin_{CD8}[10^6 cells/day]');zlabel('Final tumor volume  [mL]');
title('Association between T-cell parameters and final tumor volumes')
grid on
%% Comparing ORR
ORR = {HumanPI(1:end).Response};
parametrizations=num2cell(1:length(HumanPI));
parametrizationNames =   cellfun(@(x)strjoin({'Par' num2str(x)},''),parametrizations,'UniformOutput',false);
controlORR = array2table(cell2mat(cellfun(@(x) x{:,2}, ORR,...
    'UniformOutput', false)), 'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames);
antiPDL1ORR = array2table(cell2mat(cellfun(@(x) x{:,3}, ORR, ...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames );
antiCTLA4ORR = array2table(cell2mat(cellfun(@(x) x{:,4}, ORR,...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames) ;
antiPDL1_antiCTLA4ORR = array2table(cell2mat(cellfun(@(x) x{:,5}, ORR, ...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames);
antiPDL1ORR{:,:}-controlORR{:,:}
antiCTLA4ORR{:,:}-controlORR{:,:}
antiPDL1_antiCTLA4ORR{:,:}-controlORR{:,:}
%% 
pfs = arrayfun(@(x)getMedianPFS(x.PI,treatments),HumanPI,'UniformOutput',false);
[HumanPI(1:end).PFS]=pfs{:,:};
os = arrayfun(@(x)getMedianPFS(x.PI,treatments,'survivalType','OS'),HumanPI,'UniformOutput',false);
[HumanPI(1:end).OS]=os{:,:};
%% Plot predictions of ORR against realized values
RECIST_Siu = table({'PD' 'SD' 'PR' 'CR'}', 'VariableNames', {'Response'});
RECIST_Siu(:, end+1) = {nan(1,1)};
RECIST_Siu(:, end+1) = {64.6  6.2     9.2 0}';
RECIST_Siu(:, end+1) = {69.8  0       1.6 0}';
RECIST_Siu(:, end+1) = {64.3  5.4     7.8 0}';
RECIST_Siu.Properties.VariableNames(2:end) = treatments;

ORR_Siu = table(treatments','VariableNames', {'Treatment'});
ORR_Siu{:,end+1} = [nan 9.2      1.6     7.8]';
ORR_Siu{:,end+1} = [nan 3.46     0.04    3.78]';
ORR_Siu{:,end+1} = [nan 19.02    8.53    13.79]';
ORR_Siu.Properties.VariableNames = {'Treatments' 'ORR' 'CI_LB','CI_UB'};

ORR_Ferris = table(treatments','VariableNames', {'Treatment'});
ORR_Ferris{:,end+1} = [17.3 17.9    nan    18.2]';
ORR_Ferris{:,end+1} = [12.8 13.3    nan    13.6]';
ORR_Ferris{:,end+1} = [22.5 23.4    nan    23.6]';
ORR_Ferris.Properties.VariableNames = {'Treatments' 'ORR' 'CI_LB','CI_UB'};

colors=linspecer(8);
ORR_Preclinical = HumanPI(1).ORR;
markers = {'d'};

figure
hold on
for i=1:length(treatments)
    indx = (i-1)*2+2;
    CI = cellfun(@(x)str2double(x), ORR_Preclinical{:, indx+1});
    
    h = errorbar(ORR_Siu{i,2}, ORR_Preclinical{3, indx}*100,...
        (ORR_Preclinical{3, indx}-CI(3,1))*100,(CI(3,2)-ORR_Preclinical{3, indx})*100,...
        ORR_Siu{i,2}-ORR_Siu{i,3},ORR_Siu{i,4}-ORR_Siu{i,2});
    h1 = errorbar(ORR_Ferris{i,2}, ORR_Preclinical{3, indx}*100,...
        (ORR_Preclinical{3, indx}-CI(3,1))*100,(CI(3,2)-ORR_Preclinical{3, indx})*100,...
        ORR_Ferris{i,2}-ORR_Ferris{i,3},ORR_Ferris{i,4}-ORR_Ferris{i,2});
    try
    h.Color = colors(i,:);
    h.LineStyle = 'none';
    h.LineWidth = 2;
    h.Marker = markers{1};
    h.MarkerEdgeColor = colors(i,:);
    h.MarkerFaceColor = colors(i,:);
    h.MarkerSize = 15;
    catch
    end
    try
    h1.Color = colors(i+4,:);
    h1.LineStyle = 'none';
    h1.LineWidth = 2;
    h1.Marker = markers{1};
    h1.MarkerEdgeColor = colors(i+4,:);
    h1.MarkerFaceColor = colors(i+4,:);
    h1.MarkerSize = 15;
    
    catch
    end
end
ax = gca;
plot(-5:1:30, -5:1:30, '--k','LineWidth',2)
xlim([0 25])
ylim([0 25])
grid on
legendNames = [cellfun(@(x)strjoin({x,'(Ferris 2020)'},' '),treatments(end:-1:1), 'UniformOutput',false)
    cellfun(@(x)strjoin({x,'(Siu 2020)'},' '),treatments(end:-1:1), 'UniformOutput',false)];

legend(ax.Children(2:1:end),legendNames,'interpreter', 'none')
title('RECIST 1.1 ORR comparison between simulations and realized clinical values')
xlabel('Clinical estimates [% of patients]')
ylabel('Clinical simulation [% of simulations]')
%% Median PFS
PFS_Siu = table(treatments');
PFS_Siu(:,end+1) = {nan 1.9 2 1.9}';
PFS_Siu(:,end+1) = {nan 1.8 1.8 1.9}';
PFS_Siu(:,end+1) = {nan 2.8 2 2.1}';
PFS_Preclinical = HumanPI(1).PFS;

PFS_Ferris = table(treatments');
PFS_Ferris(:,end+1) = {3.7 2.1 nan 2}';
PFS_Ferris(:,end+1) = {3.1 1.9 nan 1.9}';
PFS_Ferris(:,end+1) = {3.8 3   nan 2.3}';
PFS_Ferris.Properties.VariableNames = {'Treatments' 'PFS' 'CI_LB','CI_UB'};

figure
hold on
for i=1:4

h=errorbar(PFS_Siu{i,2}, PFS_Preclinical{i,2},PFS_Preclinical{i,2}-PFS_Preclinical{i,3},...
    PFS_Preclinical{i,4}-PFS_Preclinical{i,2}, PFS_Siu{i,2}-PFS_Siu{i,3},...
    PFS_Siu{i,4}-PFS_Siu{i,2});
h.MarkerFaceColor = colors(i,:);
h1=errorbar(PFS_Ferris{i,2}, PFS_Preclinical{i,2},  PFS_Preclinical{i,2}-PFS_Preclinical{i,3},...
    PFS_Preclinical{i,4}-PFS_Preclinical{i,2},PFS_Ferris{i,2}-PFS_Ferris{i,3},...
    PFS_Ferris{i,4}-PFS_Ferris{i,2});
try
    h.Color = colors(i,:);
    h.LineStyle = 'none';
    h.LineWidth = 2;
    h.Marker = markers{1};
    h.MarkerEdgeColor = colors(i,:);
    h.MarkerFaceColor = colors(i,:);
    h.MarkerSize = 15;
catch
end
try
    h1.Color = colors(i+4,:);
    h1.LineStyle = 'none';
    h1.LineWidth = 2;
    h1.Marker = markers{1};
    h1.MarkerEdgeColor = colors(i+4,:);
    h1.MarkerFaceColor = colors(i+4,:);
    h1.MarkerSize = 15;
    
catch
end
end
ax = gca;
plot(0:1:10, 0:1:10, '--k','LineWidth',2)
xlim([0 10])
ylim([0 10])
grid on
legendNames = [cellfun(@(x)strjoin({x,'(Ferris 2020)'},' '),treatments(end:-1:1), 'UniformOutput',false)
    cellfun(@(x)strjoin({x,'(Siu 2020)'},' '),treatments(end:-1:1), 'UniformOutput',false)];

legend(ax.Children(2:1:end),legendNames,'interpreter', 'none')
title('PFS comparison between simulations and realized clinical values')
xlabel('Clinical estimates [months]')
ylabel('Clinical simulation [months]')
%% PFS curves
Zandberg_antiPDL1 = readtable("C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Zandberg_2018_antiPDL1.csv");
Siu_antiPDL1 = readtable("C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Siu_2019_antiPDL1.csv");
Siu_antiCTLA4_antiPDL1 = readtable("C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Siu_2019_antiCTLA4_antiPDL1.csv");
Siu_antiCTLA4 = readtable("C:\Users\jmten\OneDrive\Dokumente\GitHub\CIT-SimBiology-Toolbox\data\Siu_2019_antiCTLA4.csv");

% Siu_antiCTLA4_antiPDL1 = [ 133 39 16 13 5 1 0];
% Siu_antiPDL1 = [67 24 12 8 5 0 0];
% Siu_antiCTLA4 = [67 12 1 1 0 0 0];
% Siu_time = [0 3 6 9 12 15 18]*30;
% 
% S_t_Siu_antiCTLA4_antiPDL1 = [1 -diff(Siu_antiCTLA4_antiPDL1)/Siu_antiCTLA4_antiPDL1(1)] ;
% S_t_Siu_antiPDL1 = [1 -diff(Siu_antiPDL1)/Siu_antiPDL1(1)] ;
% S_t_Siu_antiCTLA4 = [1 -diff(Siu_antiCTLA4)/Siu_antiCTLA4(1)] ;
% 
% Zandberg_antiPDL1 = [112 51 27 18 10 3 2 0 0];
% Zandberg_time = [0 3 6 9 12 15 18 21 24]*30;
% S_t_Zandberg_antiPDL1 = [1 -diff(Zandberg_antiPDL1)/Zandberg_antiPDL1(1)] ;
% 
for i=1:length(HumanPI)
    plotSurvivalFunction(HumanPI(i).PI,30*23,treatments)
    ax = gca;
    colors = [ax.Children(6).Color; ax.Children(4).Color; ax.Children(2).Color];
    hold on
    % plot(Zandberg_antiPDL1{:,1}*30, Zandberg_antiPDL1{:,2}, 'Color', colors(1,:),'LineWidth',2,'LineStyle','-.')
    plot(Siu_antiPDL1{:,1}*30, Siu_antiPDL1{:,2},'Color', colors(1,:),'LineWidth',2,'LineStyle','--')
    plot(Siu_antiCTLA4_antiPDL1{:,1}*30, Siu_antiCTLA4_antiPDL1{:,2},'Color', colors(3,:),'LineWidth',2,'LineStyle','--')
    plot(Siu_antiCTLA4{:,1}*30, Siu_antiCTLA4{:,2},'Color', colors(2,:),'LineWidth',2,'LineStyle','--')
    ICB_PFS_legends= { 'anti-PD-L1 (Siu et al., 2019)'...
        'anti-CTLA-4 (Siu et al., 2019)' 'anti-CTLA-4 + anti-PD-L1(Siu et al., 2019)'};
    ax.Legend.String = [ax.Legend.String(1:4) ICB_PFS_legends];
end

%% Calculate PFS MSE

antiPDL1_time=round(Siu_antiPDL1{:,1}*30);
antiCTLA4_time=round(Siu_antiCTLA4{:,1}*30);
antiCTLA4_antiPDL1_time=round(Siu_antiCTLA4_antiPDL1{:,1}*30);

dataTime = unique([antiPDL1_time; antiCTLA4_time; antiCTLA4_antiPDL1_time]);
simTime=0:1:(365*2);

timeIndx_antiPDL1 = ismember(simTime,antiPDL1_time);
timeIndx_antiCTLA4 = ismember(simTime,antiCTLA4_time);
timeIndx_antiCTLA4_antiPDL1 = ismember(simTime,antiCTLA4_antiPDL1_time);

MSE = table({'antiPDL1'; 'antiCTLA4';'antiCTLA4_antiPDL1'});


for i=1:length(HumanPI)
    S_t_antiPDL1 = HumanPI(i).PI.output(2).S_t;
    S_t_antiCTLA4 =  HumanPI(i).PI.output(3).S_t;
    S_t_antiCTLA4_antiPDL1 =  HumanPI(i).PI.output(4).S_t;
    
    MSE{1,i+1} = mean((S_t_antiPDL1(timeIndx_antiPDL1) - Siu_antiPDL1{:,2}').^2);
    MSE{2,i+1} = mean((S_t_antiCTLA4(timeIndx_antiCTLA4) - Siu_antiCTLA4{:,2}').^2);
    MSE{3,i+1} = mean((S_t_antiCTLA4( timeIndx_antiCTLA4_antiPDL1) - Siu_antiCTLA4_antiPDL1{:,2}').^2,'omitnan');

end
mean(MSE{:,2:end},1,'omitnan')
%% Calculate median PFS MSE
PFS_MSE = table({'antiPDL1'; 'antiCTLA4';'antiCTLA4_antiPDL1'},'VariableNames',{'Treatments'});
for i=1:length(HumanPI)
    PFS_MSE{1,i+1} = (HumanPI(i).PFS{2,2}-PFS_Siu{2,2}).^2;
    PFS_MSE{2,i+1} = (HumanPI(i).PFS{3,2}-PFS_Siu{3,2}).^2;
    PFS_MSE{3,i+1} = (HumanPI(i).PFS{4,2}-PFS_Siu{4,2}).^2;

end
mean(PFS_MSE{:,2:end},1)
PFS_MSE_Ferris = table({'antiPDL1'; 'antiCTLA4_antiPDL1'},'VariableNames',{'Treatments'});
for i=1:length(HumanPI)
    PFS_MSE_Ferris{1,i+1} = (HumanPI(i).PFS{2,2}-PFS_Ferris{2,2}).^2;
    PFS_MSE_Ferris{2,i+1} = (HumanPI(i).PFS{4,2}-PFS_Ferris{4,2}).^2;
end
mean(PFS_MSE_Ferris{:,2:end},1)
%% Calculate median ORR MSE
ORR_MSE = table({'antiPDL1'; 'antiCTLA4';'antiCTLA4_antiPDL1'},'VariableNames',{'Treatments'});
for i=1:length(HumanPI)
    ORR_sim = HumanPI(i).ORR{:,treatments(2:end)};
    ORR_MSE{1,i+1} = (ORR_sim(3,1)+ORR_sim(4,1)-ORR_Siu{2,2}).^2;
    ORR_MSE{2,i+1} = (ORR_sim(3,2)+ORR_sim(4,2)-ORR_Siu{3,2}).^2;
    ORR_MSE{3,i+1} = (ORR_sim(3,3)+ORR_sim(4,3)-ORR_Siu{4,2}).^2;

end
mean(ORR_MSE{:,2:end},1)


%% Median OS
OS_Siu = table(treatments');
OS_Siu(:,end+1) = {nan 6       5.5 7.6}';
OS_Siu(:,end+1) = {nan 4       3.9 4.9}';
OS_Siu(:,end+1) = {nan 11.3    7   10.6}';
OS_Preclinical = HumanPI(1).OS;

OS_Ferris = table(treatments');
OS_Ferris(:,end+1) = {8.3 7.6  nan 6.5}';
OS_Ferris(:,end+1) = {7.3 6.1  nan 5.5}';
OS_Ferris(:,end+1) = {9.2 9.8  nan 8.2}';
OS_Ferris.Properties.VariableNames = {'Treatments' 'PFS' 'CI_LB','CI_UB'};

figure
hold on
for i=1:4

h=errorbar(OS_Siu{i,2}, OS_Preclinical{i,2},OS_Preclinical{i,2}-OS_Preclinical{i,3},...
    OS_Preclinical{i,4}-OS_Preclinical{i,2}, OS_Siu{i,2}-OS_Siu{i,3},...
    OS_Siu{i,4}-OS_Siu{i,2});
h.MarkerFaceColor = colors(i,:);
h1=errorbar(OS_Ferris{i,2}, OS_Preclinical{i,2},  OS_Preclinical{i,2}-OS_Preclinical{i,3},...
    OS_Preclinical{i,4}-OS_Preclinical{i,2},OS_Ferris{i,2}-OS_Ferris{i,3},...
    OS_Ferris{i,4}-OS_Ferris{i,2});
try
    h.Color = colors(i,:);
    h.LineStyle = 'none';
    h.LineWidth = 2;
    h.Marker = markers{1};
    h.MarkerEdgeColor = colors(i,:);
    h.MarkerFaceColor = colors(i,:);
    h.MarkerSize = 15;
catch
end
try
    h1.Color = colors(i+4,:);
    h1.LineStyle = 'none';
    h1.LineWidth = 2;
    h1.Marker = markers{1};
    h1.MarkerEdgeColor = colors(i+4,:);
    h1.MarkerFaceColor = colors(i+4,:);
    h1.MarkerSize = 15;
    
catch
end
end
ax = gca;
plot(0:1:30, 0:1:30, '--k','LineWidth',2)
xlim([0 30])
ylim([0 30])
grid on
legendNames = [cellfun(@(x)strjoin({x,'(Ferris 2020)'},' '),treatments(end:-1:1), 'UniformOutput',false)
    cellfun(@(x)strjoin({x,'(Siu 2020)'},' '),treatments(end:-1:1), 'UniformOutput',false)];

legend(ax.Children(2:1:end),legendNames,'interpreter', 'none')
title('OS comparison between simulations and realized clinical values')
xlabel('Clinical estimates [months]')
ylabel('Clinical simulation [months]')


