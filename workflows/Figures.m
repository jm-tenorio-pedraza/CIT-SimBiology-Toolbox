%% Figures

%% Figure 6-1
PI_PK=[];
PI_PK(1).struct =  load(strjoin({cd 'PI_PK_TwoComp4_6_wo_TMDD_0.mat'},'/'));
PI_PK(2).struct = load(strjoin({cd 'PI_PK_TwoComp4_6_TMDD_0.mat'},'/'));
PI_PK(3).struct =load(strjoin({cd 'PI_PK_ThreeComp4_6_wo_TMDD_0.mat'},'/'));
PI_PK(4).struct = load(strjoin({cd 'PI_PK_ThreeComp4_6_TMDD_0.mat'},'/'));
figure('Position', [10 10 1.5e3 1e3])
for i=1:4
subplot(2,2,i)
PI = PI_PK(i).struct.PI;
plotFit(PI,'newFig', false)
ax=gca;
ax.FontSize = 14;
xlim([1e-2 1e2])
ylim([1e-2 1e2])
end

%% Figure 6-2
PI_full = load(strjoin({cd 'PI_PK_ThreeComp4_6_TMDD_0.mat'},'/'));
figure('Position', [10 10 1.5e3 1e3])
subplot(3,2,[1 4])
pc = plotPSS(PI_full.PI.pcs,4,PI_full.PI.paramNames(1:10),'threshold',-1,'newFig', false);
plotPRCC(PI_full.PI,PI_full.PI.paramNames(1:10),PI_full.PI.observablesPlot,'time',...
    1:1:PI.tspan(end),'output', 'mean','kpi', 'max', 'nrow', 2, 'ncol', 3,...
    'plotIndx', [2 3 5 6],'newFig', false,'timeUnit', 'hours')

%% Figure 6-3
PI_red = load(strjoin({cd 'PI_PK_ThreeComp4_6_TMDD_reduced_2_0.mat'},'/'));
figure('Position', [10 10 1.5e3 1e3])
for i=1:4
 subplot(2,2,i)
 plotError2(PI_red.PI,i,'all', false, 'indiv', false, 'addErrorVar', false, 'newFig', ...
     false, 'group', 'Cell', 'TimeUnit', 'hours')
end

%% Figure 6-4
PI_red = load(strjoin({cd 'PI_PK_ThreeComp4_6_TMDD_reduced_2_2.mat'},'/'));
figure('Position', [10 10 1.5e3 1e3])
paramNames = getParamNames(PI_red,sim, observables);

for i=1:4
 subplot(2,2,i)
 plotSimOutput(PI_red.PI,i,'all', false, 'indiv', false, 'addErrorVar', false,...
     'newFig', false, 'TimeUnit', 'hours')
%  set(gca, 'YScale','log')
end

%% Figure 6-5


%% Figure 6-12
PI_CIM = load(strjoin({cd 'PI_CIM22_Control_Reduced_1_0.mat'},'/'),'PI');
ncol = 3;
nrow=4;
figure('Position', [10 10 1.5e3 1e3])
subplot(nrow,nrow,[1:3])
pc = plotPSS(PI_CIM.PI.pcs,4,PI_CIM.PI.paramNames(PI_CIM.PI.H.PopulationParams),'threshold',-1,'newFig', false);
ax=gca;
ax.FontSize = 14;
ylabel('Relative PC weight')
plotPRCC(PI_CIM.PI,PI_CIM.PI.paramNames(PI_CIM.PI.H.PopulationParams),PI_CIM.PI.observablesPlot,'time',...
    1:1:PI_CIM.PI.tspan(end),'output', 'mean','kpi', 'max', 'nrow', nrow, 'ncol',ncol,...
    'plotIndx', [4:12],'newFig', false,'timeUnit', 'days')
