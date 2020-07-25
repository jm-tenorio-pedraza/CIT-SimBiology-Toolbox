function output=getOutput(PI,simFun,p,getPhi,normIndx,H,varargin)
par = inputParser;
par.addParameter('prob', 0.95)
par.addParameter('n_samples',1e3)
par.addParameter('output', PI)
par.addParameter('simTime', PI.tspan)
par.addParameter('logTransform',true)

par.parse(varargin{:})
par=par.Results;
nVar=size(PI.data(1).dataValue,2);

% Generate parameter structure
try
    phi=getPhi(p);
catch
    phi = [repmat(p, length(PI.data),1) PI.x_0];
end
% Simulate model with parameter structure
simdata=simFun(phi);

% Obtain simulation output at pre-designated time points
[T,Y,~]=getdata(simdata);
% Obtain simulation output at pre-designated time points
simdata=resample(simdata,par.simTime);
[T_data,Y_data,~]=getdata(simdata);
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

simOutput=arrayfun(@(x)x.simValue./repmat([ones(1,nVar-length(normIndx)) x.simValue(x.simTime==x.dataTime(end),normIndx)],...
        size(x.simValue(:,1),1),1),PI.data,...
        'UniformOutput',false);
dataOutput=arrayfun(@(x)x.y_hat./repmat([ones(1,nVar-length(normIndx)) x.y_hat(end,normIndx)],...
        size(x.y_hat(:,1),1),1),PI.data,...
        'UniformOutput',false);

% Input into data array
[PI.data(1:length(simOutput)).('simOutput')]=simOutput{:,:};

simTimeIndx = arrayfun(@(x) uniqueIndx(x.simTime), PI.data, 'UniformOutput', false);
[PI.data(1:length(simOutput)).simTimeIndx] = simTimeIndx{:,:};

simOutput = repelem({'nan'}, length(PI.data),1);
for i=1:length(PI.data)
    dataMat = nan(length(par.simTime), size(PI.data(1).dataValue,2));
    destinyRowIndx = ismember(par.simTime,PI.data(i).simTime(PI.data(i).simTimeIndx));
    [originRowIndx, ~] = ismember(PI.data(i).simTime(PI.data(i).simTimeIndx), par.simTime);
    dataMat(destinyRowIndx,:) = PI.data(i).simOutput(originRowIndx,:);
    simOutput{i} = dataMat;
end
% simOutput=arrayfun(@(x)x.simOutput(ismember(x.simTime(x.simTimeIndx),par.simTime),:),PI.data,...
%     'UniformOutput',false);
[PI.data(1:length(simOutput)).('simOutput')]=simOutput{:,:};
simTime = arrayfun(@(x)x.simTime(ismember(x.simTime(x.simTimeIndx),par.simTime),:),PI.data,...
    'UniformOutput',false);
[PI.data(1:length(simOutput)).simTime]=simTime{:,:};

% Input into data array
[PI.data(1:length(dataOutput)).('yOutput')]=dataOutput{:,:};

% Match observed time points
dataOutput=arrayfun(@(x)x.yOutput(ismember(par.simTime,x.dataTime),:),PI.data,...
    'UniformOutput',false);
[PI.data(1:length(dataOutput)).('y_hat')]=dataOutput{:,:};

try
% Get lower and upper boudariestry
sigma = p(setdiff(H.SigmaParams, [H.CellParams(:).OmegaIndex H.IndividualParams(:).OmegaIndex]));
if par.logTransfrom
lb = arrayfun(@(x)quantile(exp(log(x.simOutput)+ randn([size(x.simOutput),par.n_samples]).*sigma),1-par.prob,3),...
    PI.data,'UniformOutput',false);
ub = arrayfun(@(x)quantile(exp(log(x.simOutput)+ randn([size(x.simOutput),par.n_samples]).*sigma),par.prob,3),...
    PI.data,'UniformOutput',false);

else
    lb = arrayfun(@(x)quantile(((x.simOutput)+ randn([size(x.simOutput),par.n_samples]).*sigma),1-par.prob,3),...
    PI.data,'UniformOutput',false);
ub = arrayfun(@(x)quantile(((x.simOutput)+ randn([size(x.simOutput),par.n_samples]).*sigma),par.prob,3),...
    PI.data,'UniformOutput',false);
    
end
[PI.data(1:end).lb]= lb{:,:};
[PI.data(1:end).ub] = ub{:,:};
catch
end
if strcmp(par.output, 'PI')
    output = PI;
elseif strcmp(par.output, 'data')
    output = dataOutput;
else
    output = simOutput;
end
return
