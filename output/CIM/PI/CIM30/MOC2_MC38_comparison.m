PI_MOC2 = load(strjoin({cd 'PI_CIM30_MOC2_red.mat'},'\'),'PI');
PI_MC38 = load(strjoin({cd 'PI_CIM30_MC38_red.mat'},'\'),'PI');
PI_MOC2 = PI_MOC2.PI;
PI_MC38 = PI_MC38.PI;
popIndx = PI.H.PopulationParams;
popIndx_MC38 = PI_MC38.H.PopulationParams;
popIndx_MOC2 = PI_MOC2.H.PopulationParams;
mc38_indx = ismember({PI.par(PI.H.PopulationParams).name}, {PI_MC38.par(popIndx_MC38).name});
moc2_indx = ismember({PI.par(PI.H.PopulationParams).name}, {PI_MOC2.par(popIndx_MOC2).name});

phat_mc38 = [PI_MC38.par(:).finalValue];
phat_moc2 = [PI_MOC2.par(:).finalValue];
phat = finalValues(popIndx);

phat1 = exp(phat);
phat1(mc38_indx) = phat_mc38(popIndx_MC38);

phat2 = exp(phat);
phat2(moc2_indx) = phat_moc2(popIndx_MOC2);

table({PI.par(popIndx).name}', (phat1(popIndx)' - phat2(popIndx)'),...
    phat1(popIndx)'./exp(phat(popIndx)'),...
    phat2(popIndx)'./exp(phat(popIndx))','VariableNames',...
    {'Name' 'Difference' 'Ratio of MC38 wrt finalValues' 'Ration of MOC2 wrt finalValues'})

cellParams = {PI.par(PI.H.PopulationParams).name};
cellParams = cellParams(and(phat1'./exp(phat')~=1, phat2'./exp(phat')~=1));
cellIndx = ismember({PI.par(popIndx).name}, cellParams);

w = zeros(length(cellParams)*2,1);
w(1:2:end) =  log(phat1(cellIndx)) - phat(cellIndx);
w(2:2:end) = log(phat2(cellIndx)) -  phat(cellIndx);

psi = ones(length(cellParams),1);
for i=1:length(cellParams)
    startIndx = (i-1)*2 +1;
    psi(i) = std(w(startIndx:startIndx+1));
end
sigma= phat(PI.H.SigmaParams(end-7:end));

finalValues(popIndx) = phat;
finalValues([PI.H.CellParams.Index]) = w;
finalValues(PI.H.SigmaParams) = [psi' sigma];

