function PI = getDoseRegimen(PI,varargin)
inputs = inputParser;
inputs.addParameter('n_doses_antiPDL1_mono', 27);
inputs.addParameter('n_doses_antiCTLA4_mono', 9);
inputs.addParameter('dosing_times_antiPDL1_mono', 0:14:365)
inputs.addParameter('dosing_times_antiCTLA4_mono', [0:28:28*6 (28*6+12*7):12*7:365])
inputs.addParameter('dosing_times_antiPDL1_combi',[0:28:28*3 (28*3+7*2):14:365])
inputs.addParameter('dosing_times_antiCTLA4_combi', 0:28:28*3)
inputs.addParameter('BW',77)
inputs.addParameter('MW', 1.5e5)
inputs.parse(varargin{:})
inputs = inputs.Results;



%% Doses (Durvalumab + Tremelimumab)
n_doses_antiPDL1_mono       = inputs.n_doses_antiPDL1_mono;
n_doses_antiCTLA4_mono      = inputs.n_doses_antiCTLA4_mono;
dosing_times_antiPDL1_mono  = inputs.dosing_times_antiPDL1_mono;
dosing_times_antiCTLA4_mono = inputs.dosing_times_antiCTLA4_mono;
dosing_times_antiPDL1_combi = inputs.dosing_times_antiPDL1_combi;
dosing_times_antiCTLA4_combi = inputs.dosing_times_antiCTLA4_combi;

maxDosingFreq = dosing_times_antiPDL1_mono; % Maximum dosing frequency
nDoses_max = length(maxDosingFreq);

doseScalingFactor = inputs.BW*1e-3*1/inputs.MW*1e6; % kg*g/mg*mole/g*µmol/mol

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
    repelem(1*doseScalingFactor,4,1), repelem(1*doseScalingFactor/60,4,1),...
    'VariableNames',{'Time' 'Amount' 'Rate'});

control.Properties.VariableUnits            = {'day' 'micromole' 'micromole/minute'};
antiPDL1_mono.Properties.VariableUnits      = {'day' 'micromole' 'micromole/minute'};
antiCTLA4_mono.Properties.VariableUnits     = {'day' 'micromole' 'micromole/minute'};
antiPDL1_combi.Properties.VariableUnits     = {'day' 'micromole' 'micromole/minute'};
antiCTLA4_combi.Properties.VariableUnits    = {'day' 'micromole' 'micromole/minute'};

PI.u = {};
PI.u(1,1:2) = {control control};
PI.u(2,1:2) = {antiPDL1_mono control};
PI.u(3, 1:2) = {control antiCTLA4_mono};
PI.u(4, 1:2) = {antiPDL1_combi antiCTLA4_combi};
return
