function dataOutput=getOutput2(PI,simFun,p,getPhi,normIndx,varargin)
par = inputParser;
par.addParameter('prob', 0.95)
par.addParameter('n_samples',1e3)
par.addParameter('Output', 'data')
par.parse(varargin{:})
par=par.Results;
nVar=size(PI.data(1).dataValue,2);

% Generate parameter structure
try
    phi=getPhi(p);
catch
    phi = repmat(p, length(PI.data),1);
end
% Simulate model with parameter structure
simdata=simFun(phi);

% Obtain simulation output at pre-designated time points
[T,Y,~]=getdata(simdata);
% Obtain simulation output at pre-designated time points
simdata=resample(simdata,PI.tspan);
[~,Y_data,~]=getdata(simdata);
% Incorporate simulations into data structure array
try
    [PI.data(1:length(T)).('simTime')]=T{:,:};
    [PI.data(1:length(T)).('simValue')]=Y{:,:};
    [PI.data(1:length(T)).('y_hat')]=Y_data{:,:};

catch
    T= {T};
    Y = {Y};
    Y_data = {Y_data};

    [PI.data(1:length(T)).('simTime')]=T{:,:};
    [PI.data(1:length(T)).('simValue')]=Y{:,:};
    [PI.data(1:length(T)).('y_hat')]=Y_data{:,:};

end
% Normalizing

if strcmp(par.Output, 'data')
dataOutput=arrayfun(@(x)x.y_hat./repmat([ones(1,nVar-length(normIndx)) x.y_hat(end,normIndx)],...
        size(x.y_hat(:,1),1),1),PI.data,...
        'UniformOutput',false);


elseif strcmp(par.Output, 'sim')
    dataOutput=arrayfun(@(x)x.simValue./repmat([ones(1,nVar-length(normIndx)) x.simValue(x.simTime==x.dataTime(end),normIndx)],...
        size(x.simValue(:,1),1),1),PI.data,...
        'UniformOutput',false);
end
return
