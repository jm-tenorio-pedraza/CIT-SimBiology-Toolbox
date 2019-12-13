function PI = globalSA(sim,PI,observables,varargin)
param = inputParser;
param.addParameter('sigma', 2.3026)
param.addParameter('nsamples', 1e3)
param.addParameter('time', 1:1:PI.tspan(end))

param.parse(varargin{:})

param=param.Results;

parameters = (sim.Parameters.Name);                                    % Extract parameter names from simultation
time = param.time;
if length(param.sigma)>1                                                    % Check if sigma is the default parameter
    sigma = param.sigma;
else
    sigma = repelem(param.sigma, length(parameters),1);
end
inputs = sim.Parameters.Value;
samples = generateLHSSample(parameters, log(inputs), param.nsamples,'Sigma',...
    sigma,'Distribution', 'Normal');

[~,samplerank] = sort(samples);
samplescell{size(samples,1), 1} = [];
samplescell(1:size(samples,1),1) = mat2cell(exp(samples), ones(size(samples,1),1));
PRCC = [];
[PRCC(1:size(samples,1)).samples] = samplescell{:,:};

simulations = arrayfun(@(x)getData((resample(sim(x.samples, time(end),...% array of simulations for each 
    PI.u, time), time))), PRCC, 'UniformOutput', false);

[PRCC(1:size(samples,1)).simulations] = simulations{:,:};
prcc = NaN(length(time)*length(PI.u),length(parameters),length(observables));
for i = 1:length(time)
    firstrow=arrayfun(@(y)arrayfun(@(x)x{1,1}(i,:),y.simulations,'UniformOutput',...
        false),PRCC,'UniformOutput', false)';
firstrow = reshape([firstrow{:}],[],1);
    firstrow = cell2mat(firstrow);
    for j=1:length(PI.u)
        indx =j:length(PI.u):size(firstrow,1);
        [~, staterank] = sort(firstrow(indx,:));
        rho = partialcorri( staterank,samplerank);
        prcc((j-1)*length(time)+i,:,:) = reshape(rho', 1, length(parameters),length(observables));
    end
end
PI.prcc = prcc;
PI.samples = samples;

