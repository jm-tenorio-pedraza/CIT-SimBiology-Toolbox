
%% Figure 3 PK
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_4')
load(strjoin({cd 'PI_PK_ThreeComp4_4_TMDD_11.mat'},'/'))
finalValues = [PI.par(:).finalValue];
subplot(1,2,1)
plotFit(PI,'sigma', exp(finalValues(setdiff(PI.H.SigmaParams,...
    [PI.H.IndividualParams.OmegaIndex PI.H.CellParams.OmegaIndex]))),'newFig', false)
for i=1:length(ax.Children)
    ax.Children(i).MarkerSize =10 ;
    ax.Children(i).MarkerFaceColor = 'none';
    ax.Children(i).LineWidth = 1.5;

end

subplot(1,2,2)
plotCorrMat(PI.postSamples(:,[PI.H.PopulationParams PI.H.SigmaParams]), ...
    PI.paramNames([PI.H.PopulationParams PI.H.SigmaParams]), 'interpreter','tex');
%
subplot(2,2,[1,2])
plotPSS(PI.pcs,5,PI.paramNames(PI.H.PopulationParams),'threshold',-1,'newFig', false)

plotIIVParams(PI.postSamples, PI,'name', PI.paramNames,'newFig', false,...
    'n_row', 2,'n_col',2,'figIndx', 3:4,'panel', true,'dim',true)

