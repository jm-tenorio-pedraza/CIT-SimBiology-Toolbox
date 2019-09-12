function [sensmatrix] = getSensitivities(inputs, PI,sim,parameters, observables,time)

inputs(2,:) =inputs(1,:)*1.1;
phi_i = NaN(1,length(parameters));
sensmatrix = NaN(length(observables)*length(time)*size(PI.data,1), length(parameters));
for i=1:length(parameters)
    phi_i(1,:) = inputs(1,:);
    phi_i(1,i) = inputs(2,i);
    simdata = sim(phi_i);
    simdata = resample(simdata, time);
    [~,y_i,~] = getdata(simdata);
    [PI.data(1:end).y_i] = y_i{:,:};
    y_i = arrayfun(@(x) ((log(x.y_i) - log(x.simValue)))./...
        (log(phi_i(1,i))-log(inputs(1,i))), PI.data,'UniformOutput',false);
    delta = reshape([y_i{:,:}],[],1);
    sensmatrix(:,i) = delta;
end

