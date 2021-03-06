function [sensmatrix] = getSensitivities(inputs, PI, sim, parameters,...
    observables, time,varargin)
%% inputs is a matrix whose rows are the parameters for each one of the conditions in sim


p = inputParser;
p.addParameter('allVariables', true);
p.addParameter('initialValue', false);
p.addParameter('uniqueGroups', false);
p.addParameter('timedelta', 1);
p.addParameter('parameterdelta', 0.1);

p.parse(varargin{:})
p=p.Results;
p_delta = p.parameterdelta;
PI.sensitivity = struct('simValue', repelem({nan(length(time),...
    length(observables))},size(inputs,1),1));

simdata = sim(inputs);                                                      % Simulate model and resample simulation results
simdata = resample(simdata, time);
[~,y_i,~] = getdata(simdata);
varIndx = ismember(simdata(1).DataNames, observables);                      % Detect desired outputs

try
    [PI.sensitivity(1:end).simValue] = y_i{:,:};                            % for the case when there is more than 1 output
catch
    y_i={y_i};                                                              % for the case when there is only 1 output
    [PI.sensitivity(1:end).simValue] = y_i{:,:};
end  
simValue = arrayfun(@(x) x.simValue(:, varIndx), PI.sensitivity,...         % Input simulated values under no perturbartion
    'UniformOutput', false);
[PI.sensitivity(1:end).simValue] = simValue{:,:};

groups = [PI.data(:).Group];
if ischar((groups))
    groups = {PI.data(:).Group};
end
if p.uniqueGroups
    sensmatrix = NaN(length(observables)*length(time)*length(unique(groups)), length(parameters)); % allocate spcace to the sensitivity matrix 
else
    sensmatrix = NaN(length(observables)*length(time)*size(PI.data,1), length(parameters)); % allocate spcace to the sensitivity matrix 
end
inputs_i = mat2cell(inputs, repelem(1,size(inputs,1)));

[PI.sensitivity(1:end).inputs] = inputs_i{:,:};


for i=1:length(parameters)
    phi_i = inputs;
    phi_i(:,i) = inputs(:,i)*(1+p_delta);
    simdata = sim(phi_i);               % Simulate model and resample simulation results
    simdata = resample(simdata, time);
    [~,y_i,~] = getdata(simdata);
    try
        [PI.sensitivity(1:end).y_i] = y_i{:,:}; % for the case when there is more than 1 output
    catch
        y_i={y_i};                        % for the case when there is only 1 output
         [PI.sensitivity(1:end).y_i] = y_i{:,:};
    end
    y_i = arrayfun(@(x) x.y_i(:, varIndx), PI.sensitivity, 'UniformOutput', false);
    [PI.sensitivity(1:end).y_i] = y_i{:,:};
 
    
    phi_i = mat2cell(phi_i, repelem(1,size(phi_i,1)));
    [PI.sensitivity(1:end).phi] = phi_i{:,:};

    if p.allVariables                       % Check if all output variables are to be considered
    y_i = arrayfun(@(x) (((x.y_i) - (x.simValue)))./...
        (log(x.phi(:,i))-log(x.inputs(:,i))), PI.sensitivity,'UniformOutput',false);
    end
    delta = reshape([y_i{:,:}],[],1);
    sensmatrix(:,i) = real(delta);
end
sensmatrix=sensmatrix*sqrt(p.timedelta/time(end));

if p.initialValue
    sensmatrix = sensmatrix(:,1:end-1);
end

