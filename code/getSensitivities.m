function [sensmatrix] = getSensitivities(inputs, PI,sim,parameters, observables,time,varargin)
p = inputParser;
p.addParameter('allVariables', true);
p.parse(varargin{:})
p=p.Results;

inputs(2,:) =inputs(1,:)*1.1; % Add another row to the initial parameter vector 
phi_i = NaN(1,length(parameters));
sensmatrix = NaN(length(observables)*length(time)*size(PI.data,1), length(parameters)); % allocate spcace to the sensitivity matrix 
for i=1:length(parameters)
    phi_i(1,:) = inputs(1,:);
    phi_i(1,i) = inputs(2,i);
    simdata = sim(phi_i);               % Simulate model and resample simulation results
    simdata = resample(simdata, time);
    [~,y_i,~] = getdata(simdata);
    try
        [PI.data(1:end).y_i] = y_i{:,:}; % for the case when there is more than 1 output
    catch
        y_i={y_i};                        % for the case when there is only 1 output
         [PI.data(1:end).y_i] = y_i{:,:};
    end
    if p.allVariables                       % Check if all output variables are to be considered
    y_i = arrayfun(@(x) ((log(x.y_i) - log(x.simValue)))./...
        (log(phi_i(1,i))-log(inputs(1,i))), PI.data,'UniformOutput',false);
    else
          y_i = arrayfun(@(x) ((log(x.y_i) - log(x.simValue)))./...
        (log(phi_i(1,i))-log(inputs(1,i))), PI.data,'UniformOutput',false);
   
    end
    delta = reshape([y_i{:,:}],[],1);
    sensmatrix(:,i) = delta;
end

