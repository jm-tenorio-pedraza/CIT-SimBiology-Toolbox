function [residuals,PI]=getNormResiduals(p,simFun,PI,getPhi,normIndx,varargin)
% Calculates normalized residuals with the paramaters (double) input mapping inputFun (handle)
% simFun (handle), function and the data (table) provided
par=inputParser;
par.addParameter('log',true);
par.addParameter('errorModel','additive');
par.addParameter('constantVar',0.1);
par.addParameter('missingRes','NA');
par.addParameter('indivData',false);
par.parse(varargin{:})
par=par.Results;
nVar=size(PI.data(1).dataValue,2);
if size(p,1)>size(p,2)
    p=p';
end
% Generate parameter structure
phi=getPhi(p);
% Generate Sigma structure
sigma=p(setdiff(PI.H.SigmaParams, [[PI.H.CellParams.OmegaIndex],...
    [PI.H.IndividualParams.OmegaIndex], [PI.H.RespParams.OmegaIndex]])); % 1xp vector
% Simulate model with parameter structure
try
simdata=simFun(phi);
catch ME
    if strcmp(ME.identifier,'SimBiology:SimFunction:SomeSimulationsFailed')
        residuals = inf;
        return  
    end
end
 
% Obtain simulation output at pre-designated time points
try
simdata=resample(simdata,PI.tspan);
[T,Y,~]=getdata(simdata);
catch
    residuals = inf;
    return
end
% Incorporate simulations into data structure array
try
[PI.data(1:length(T)).('simTime')]=T{:,:};
[PI.data(1:length(T)).('simValue')]=Y{:,:};
catch
    PI.data.simTime = T;
    PI.data.simValue=Y;
end
% Add simulations to individual data
if par.indivData
    indivData=repelem({PI.data(:).simValue},[PI.data(:).Count]);
    indivTime=repelem({PI.data(:).simTime},[PI.data(:).Count]);
    [PI.IndivData(1:end).simTime]=indivTime{:,:};
    [PI.IndivData(1:end).simValue]=indivData{:,:};
     % Match observed time points
     try
     simOutput=arrayfun(@(x)x.simValue(ismember(x.simTime,x.dataTime),:),PI.IndivData,...
         'UniformOutput',false);
     catch
         residuals = inf;
         return
     end
     simOutput=cellfun(@(x)x./([ones(1,nVar-length(normIndx)) x(end,normIndx)]),...
         simOutput,'UniformOutput',false);
    sigma = arrayfun(@(x)repmat(sigma, size(x.dataValue,1),1), PI.IndivData, 'UniformOutput', false);
    [PI.IndivData(1:end).y_hat]= simOutput{:,:};
    [PI.IndivData(1:end).sigma]= sigma{:,:};
else
    sigma = arrayfun(@(x)repmat(sigma, size(x.dataValue,1),1), PI.data, 'UniformOutput', false);

    % Interpolate simulations to match observed time points
    try
        simOutput=arrayfun(@(x)x.simValue(ismember(x.simTime,x.dataTime),:),PI.data,...
            'UniformOutput',false);
    catch
        try
            simOutput= PI.data.simValue(ismember(PI.data.simTime, PI.data.dataTime),:);
        catch
            residuals=inf;
            return
        end
    end
    % Normalizing to final value for DCm, ISC and PDL1
    simOutput=cellfun(@(x)x./([ones(1,nVar-length(normIndx)) x(end,normIndx)]),...
        simOutput,'UniformOutput',false);
    % Input into data array
    [PI.data(1:length(simOutput)).('y_hat')]=simOutput{:,:};
    [PI.data(1:length(simOutput)).('sigma')] = sigma{:,:};
    if strcmp(par.errorModel, 'additive')
    else
        sigma = arrayfun(@(x)par.constantVar+x.y_hat.*x.sigma,PI.data,'UniformOutput', false);
        [PI.data(1:length(simOutput)).('sigma')] = sigma{:,:};
    end
end

% Errors of log-transformed data
residuals=getErrors(PI,'log',par.log,'indivData',par.indivData);
if par.indivData
    if any([isinf(residuals), ~isreal(residuals)])
        residuals = inf;
    elseif (length(residuals)< PI.indivN_data)
        if strcmp(par.missingRes, 'supplant')
            rescompliment= repelem(max(abs(residuals)),1,PI.n_data-length(residuals));
            residuals = [residuals rescompliment];
        else
            residuals = inf;
        end
    end
else
    if any([isinf(residuals), ~isreal(residuals)])
        residuals = inf;
    elseif (length(residuals)< PI.n_data)
        if strcmp(par.missingRes, 'supplant')
            rescompliment= repelem(max(abs(residuals)),1,PI.n_data-length(residuals));
            residuals = [residuals rescompliment];
        else
            residuals = inf;
        end
    end
end
return