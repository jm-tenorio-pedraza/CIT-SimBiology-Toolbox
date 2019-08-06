function [data]=getOutput(data,simFun,p,getPhi)

    % Generate parameter structure
    phi=getPhi(p);
    
    % Simulate model with parameter structure
    simdata=simFun(phi);
    
    % Obtain simulation output at pre-designated time points
    simdata=resample(simdata,0:1:100);
    [T,Y,~]=getdata(simdata);
    
    % Incorporate simulations into data structure array
    [data(1:length(T)).('simTime')]=T{:,:};
    [data(1:length(T)).('simValue')]=Y{:,:};
    
    % Normalizing to final value for DCm, ISC and PDL1
    simOutput=arrayfun(@(x)x.simValue./repmat([ones(1,3) x.simValue(x.simTime==40,4:6)],...
        size(x.simValue,1),1),data,'UniformOutput',false);
    
    % Input into data array
    [data(1:length(simOutput)).('simOutput')]=simOutput{:,:};
    
   
en
rereturn
