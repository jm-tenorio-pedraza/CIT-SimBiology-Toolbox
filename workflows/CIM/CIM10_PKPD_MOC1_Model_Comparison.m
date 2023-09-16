%%
if ispc
    wd='\Users\jmten\OneDrive\Dokumente\GitHub\';
else
    wd='/Users/migueltenorio/Documents/GitHub/';
end
addpath(genpath(strjoin({wd 'CIT-SimBiology-Toolbox'}, '')))
cd(strjoin({wd 'CIT-SimBiology-Toolbox/output'},''))
out=sbioloadproject(strjoin({wd 'sbio-projects/CIM_10_1_PKPD.sbproj'},''));
data_ext = strjoin({ wd 'CIT-SimBiology-Toolbox\data\PI_Clavijo_2.mat'},'/');
data_ext1 = strjoin({wd 'CIT-SimBiology-Toolbox\data\PI_Morisada_3.mat'},'/');

PI_kinCD8=load(strjoin({cd '\CIM\PI\CIM01\MCMC\PI_CIM10_PKPD_MOC1_MCMC_kin_CD8.mat'},'/'),'PI');
PI_kproTumor=load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_kpro_Tumor.mat'},'/'),'PI');
PI_killCD8=load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_kill_CD8.mat'},'/'),'PI');
PI_KDEMDSC=load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_KDE_MDSC.mat'},'/'),'PI');
PI_kin_Treg=load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_kin_Treg.mat'},'/'),'PI');
PI_KDETreg=load(strjoin({cd 'CIM/PI/CIM33/PI_CIM10_PKPD_MOC1_KDE_Treg.mat'},'/'),'PI');

AIC=[PI_kinCD8.PI.AIC PI_kproTumor.PI.AIC PI_killCD8.PI.AIC  PI_KDEMDSC.PI.AIC  PI_kin_Treg.PI.AIC PI_KDETreg.PI.AIC];