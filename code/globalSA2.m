function PI = globalSA2(sim,PI, observables,varargin)
param = inputParser;
param.addParameter('sigma', 2.3026)
param.addParameter('nsamples', 1e3)
param.addParameter('time', 1:1:PI.tspan(end))
param.addParameter('inputs', sim.Parameters.Value)
param.addParameter('variation', 0.9)
param.addParameter('LB', [])
param.addParameter('UB', [])
param.addParameter('samples', [])

param.parse(varargin{:})

param=param.Results;
popParamsIndx = [PI.H.PopulationParams];
parameters = {PI.par(popParamsIndx).name}';                                    % Extract parameter names from simultation
time = param.time;
if isempty(param.samples)
[samples, ~] = generateLHSSample2(PI,param.inputs,  param.nsamples,...
    'Variation', param.variation,'LB', param.LB, 'UB', param.UB);

else
    samples = param.samples;
end
samples_norm = samples - mean(samples);

samplescell{size(samples,1), 1} = [];
samplescell(1:size(samples,1),1) = mat2cell([repelem(samples, length(PI.u),1)],...
    ones(size(samples,1),1)*(length(PI.u)));
PRCC = [];
[PRCC(1:size(samples,1)).samples] = samplescell{:,:};
u = PI.u;
p = gcp('nocreate'); % If no pool, do not create new one.

if isempty(p)
    p = parpool('local')
    pctRunOnAll warning off
else
end
parfor i=1:size(samples,1)
    PRCC(i).simulations = getData(resample(sim(samples(i,:), time(end),...% array of simulations for each 
    u, time), time));
end
% simulations = arrayfun(@(x)getData(resample(sim(x.samples, time(end),...% array of simulations for each 
%     PI.u, time), time)), PRCC, 'UniformOutput', false);

% [PRCC(1:size(samples,1)).simulations] = simulations{:,:};
prcc = NaN(length(time)*length(PI.u),length(parameters),length(observables));
for i = 1:length(time)
    timerow_i=arrayfun(@(y)arrayfun(@(x)x{1,1}(i,:),y.simulations,'UniformOutput',...
        false),PRCC,'UniformOutput', false)';
    timerow_i = reshape([timerow_i{:}],[],1);
    timerow_i = cell2mat(timerow_i);
    for j=1:length(PI.u)
        indx =j:length(PI.u):size(timerow_i,1);
        timerow_j = timerow_i(indx,:);
%         [stateSorted] = sort(firstrow_j,'ascend');
%         staterank = zeros(param.nsamples, length(observables));
%         for k=1:length(observables)
%             [~,staterank(:,k)] = ismember( firstrow_j(:,k),stateSorted(:,k));
%         end
        timerow_j_norm = timerow_j - mean(timerow_j,'omitnan');
        rho = partialcorri(timerow_j_norm,samples_norm,'Rows', 'complete',...
            'Type', 'Spearman');
        prcc((j-1)*length(time)+i,:,:) = reshape(rho', 1, length(parameters),...
            length(observables));
    end
end
PI.prcc = prcc;
PI.samples = samples;
end

