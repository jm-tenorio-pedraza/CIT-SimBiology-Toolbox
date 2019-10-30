function residuals = getResiduals(p,simFun,PI,getPhi,sigma,normIndx,varargin)
% Calculates normalized residuals with the paramaters (double) input mapping inputFun (handle)
% simFun (handle), function and the data (table) provided


par=inputParser;
par.addParameter('addSigma',false)
par.addParameter('output','residuals')
par.parse(varargin{:});
par=par.Results;

nVar=size(PI.data(1).dataValue,2);
if size(p,1)>size(p,2)
    p=p';
end
% Generate parameter structure
phi=getPhi(p);


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

simOutput=cellfun(@(x)x./repmat([ones(1,nVar-length(normIndx)) x(end,normIndx)],size(x(:,1),1),1),simOutput,'UniformOutput',false);

% Input into data array
[PI.data(1:length(simOutput)).('y_hat')]=simOutput{:,:};

% Errors of log-transformed data
if par.addSigma
    error=arrayfun(@(x)reshape((log(x.y_hat)-log(x.dataValue)).^2./...% squared residuals
    ((2*sigma.^2))+...% normalized by their variance
    (log(sigma*sqrt(pi*2))),1,[]),PI.data,'UniformOutput', false);% normalized by their variance

else
error=arrayfun(@(x)reshape((log(x.y_hat)-log(x.dataValue))./...% squared residuals
    (sqrt(2)*(sigma)),1,[]),PI.data,'UniformOutput', false);% normalized by their variance
end

error=cellfun(@(x)x(~isnan(x)),error,'UniformOutput',false);
[PI.data(1:end).residuals] = error{:,:};
residuals=[error{:,:}];


if (length(residuals)~= PI.n_data)
    residuals=repelem(1e5, 1,PI.n_data);
elseif any([isinf(residuals), ~isreal(residuals)])
    residuals=repelem(1e5, 1,PI.n_data);
end

if strcmp(par.output, 'residuals')
else 
    residuals = PI;
end
return