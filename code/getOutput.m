function PI=getOutput(PI,simFun,p,getPhi,normIndx,H,varargin)
par = inputParser;
par.addParameter('prob', 0.95)
par.addParameter('n_samples',1e3)

par.parse(varargin{:})
par=par.Results;
nVar=size(PI.data(1).dataValue,2);

% Generate parameter structure
phi=getPhi(p);

% Simulate model with parameter structure
simdata=simFun(phi);

% Obtain simulation output at pre-designated time points
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

simOutput=arrayfun(@(x)x.simValue./repmat([ones(1,nVar-length(normIndx)) x.simValue(x.simTime==x.dataTime(end),normIndx)],...
        size(x.simValue(:,1),1),1),PI.data,...
        'UniformOutput',false);

% Input into data array
[PI.data(1:length(simOutput)).('simOutput')]=simOutput{:,:};



% Interpolate simulations to match observed time points
simOutput=arrayfun(@(x)x.simOutput(ismember(x.simTime,x.dataTime),:),PI.data,...
    'UniformOutput',false);


% Input into data array
[PI.data(1:length(simOutput)).('y_hat')]=simOutput{:,:};

% Get lower and upper boudaries
sigma = p(setdiff(H.SigmaParams, [H.CellParams(:).OmegaIndex H.IndividualParams(:).OmegaIndex]));


lb = arrayfun(@(x)quantile(exp(log(x.simOutput)+ randn([size(x.simOutput),par.n_samples]).*sigma),1-par.prob,3),...
    PI.data,'UniformOutput',false);
ub = arrayfun(@(x)quantile(exp(log(x.simOutput)+ randn([size(x.simOutput),par.n_samples]).*sigma),par.prob,3),...
    PI.data,'UniformOutput',false);
[PI.data(1:end).lb]= lb{:,:};
[PI.data(1:end).ub] = ub{:,:};

end
