function PI=getOutput(PI,simFun,p,getPhi,normIndx,time)
nVar=size(PI.data(1).dataValue,2);

% Generate parameter structure
phi=getPhi(p);

% Simulate model with parameter structure
simdata=simFun(phi);

% Obtain simulation output at pre-designated time points
simdata=resample(simdata,time);
[T,Y,~]=getdata(simdata);

% Incorporate simulations into data structure array
try
[PI.data(1:length(T)).('simTime')]=T{:,:};
[PI.data(1:length(T)).('simValue')]=Y{:,:};
catch
    T= {T};
    Y = {Y};
    [PI.data(1:length(T)).('simTime')]=T{:,:};
[PI.data(1:length(T)).('simValue')]=Y{:,:};

end
% Normalizing
try
    simOutput=arrayfun(@(x)x.simValue./repmat([ones(1,nVar-length(normIndx)) x.simValue(x.simTime==40,normIndx)],...
        size(x.simValue(:,1),1),1),PI.data,...
        'UniformOutput',false);
catch
    simOutput=arrayfun(@(x)x.simValue./repmat([ones(1,nVar-length(normIndx)) x.simValue(x.simTime==40,normIndx)],...
    size(x.simValue(:,1),1),1),PI.data,...
    'UniformOutput',false);
end

% Input into data array
[PI.data(1:length(simOutput)).('simOutput')]=simOutput{:,:};
end
