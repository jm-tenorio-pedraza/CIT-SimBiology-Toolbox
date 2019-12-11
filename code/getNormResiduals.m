function [residuals,PI]=getNormResiduals(p,simFun,PI,getPhi,getCovariance,normIndx,varargin)
% Calculates normalized residuals with the paramaters (double) input mapping inputFun (handle)
% simFun (handle), function and the data (table) provided
par=inputParser;
par.addParameter('log',true);
par.parse(varargin{:})
par=par.Results;
nVar=size(PI.data(1).dataValue,2);
if size(p,1)>size(p,2)
    p=p';
end
% Generate parameter structure
phi=getPhi(p);

% Generate Sigma structure
sigma=getCovariance(p); % 1xp vector

% Simulate model with parameter structure
try
simdata=simFun(phi);
catch ME
    if strcmp(ME.identifier,'SimBiology:SimFunction:SomeSimulationsFailed')
        residuals = 1e9;
        return  
    end
end
 
% Obtain simulation output at pre-designated time points
simdata=resample(simdata,PI.tspan);
[T,Y,~]=getdata(simdata);

% Incorporate simulations into data structure array
[PI.data(1:length(T)).('simTime')]=T{:,:};
[PI.data(1:length(T)).('simValue')]=Y{:,:};

% Interpolate simulations to match observed time points
try
    simOutput=arrayfun(@(x)x.simValue(ismember(x.simTime,x.dataTime),:),PI.data,...
    'UniformOutput',false);
% simOutput=arrayfun(@(x)interp1(x.simTime,x.simValue,x.dataTime),PI.data,...
%     'UniformOutput',false);
catch
    residuals=1e9;
    return
end
% Normalizing to final value for DCm, ISC and PDL1

simOutput=cellfun(@(x)x./([ones(1,nVar-length(normIndx)) x(end,normIndx)]),simOutput,'UniformOutput',false);

% Input into data array
[PI.data(1:length(simOutput)).('y_hat')]=simOutput{:,:};

% Errors of log-transformed data
residuals=getErrors(PI,sigma,'log',par.log);

if (length(residuals)~= PI.n_data)
    residuals=1e9;
elseif any([isinf(residuals), ~isreal(residuals)])
    residuals=1e9;
end
return