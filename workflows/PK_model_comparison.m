%% Compare PK models wrt number of compartments and TMDD
addpath(genpath('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox'))
cd('/Users/migueltenorio/Documents/GitHub/CIT-SimBiology-Toolbox/output/PK_mAb_ThreeComp/PI/ThreeComp_4')

PI_TwoComp_wo_TMDD = load(strjoin({cd '/PI_PK_TwoComp4_4_wo_TMDD_0.mat'},''));
PI_TwoComp_w_TMDD = load(strjoin({cd '/PI_PK_TwoComp4_4_TMDD_0.mat'},''));
PI_ThreeComp_wo_TMDD = load(strjoin({cd '/PI_PK_ThreeComp4_4_wo_TMDD_0.mat'},''));
PI_ThreeComp_w_TMDD = load(strjoin({cd '/PI_PK_ThreeComp4_4_TMDD_0.mat'},''));
MetaPI = [PI_TwoComp_wo_TMDD PI_TwoComp_w_TMDD PI_ThreeComp_wo_TMDD PI_ThreeComp_w_TMDD]

for i=1:4
subplot(2,2,i)
PI = MetaPI(i).PI;
plotFit(PI,'sigma', [PI.par((setdiff( PI.H.SigmaParams,...
    [ PI.H.IndividualParams.OmegaIndex  PI.H.CellParams.OmegaIndex]))).finalValue],...
    'newFig', false)

end