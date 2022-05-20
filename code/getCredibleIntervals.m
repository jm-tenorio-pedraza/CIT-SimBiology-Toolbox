function PI = getCredibleIntervals(PI,outputs,params,H,varargin)
inputs = inputParser;
inputs.addParameter('logit_indx', [])
inputs.addParameter('simTime', PI.tspan);
inputs.addParameter('sigma', []);
inputs.addParameter('errorModel', 'additive');
inputs.addParameter('constantVar', 0.1);

inputs.parse(varargin{:})
inputs = inputs.Results;
if isempty(inputs.sigma)
sigma=params(:,setdiff(H.SigmaParams, [H.IndividualParams.OmegaIndex H.CellParams.OmegaIndex H.RespParams.OmegaIndex]));
else
    sigma = inputs.sigma;
end

PI.CI=struct('Name',{PI.data(1:end).Name}','Group',{PI.data(1:end).Group}')';
PI.MSE = nan(length(PI.data), length(outputs));
datasetNames = {PI.output(1:end).Name};
[PI.outputPred(1:length(PI.output)).Name] = datasetNames{:,:};
for i=1:length(outputs)
    output_i=char(outputs(i));
    sigma_i=sigma(:,i);
    
    outputLB=arrayfun(@(x)quantile(x.(output_i),0.025),PI.output,'UniformOutput',false);
    outputUB=arrayfun(@(x)quantile(x.(output_i),0.975),PI.output,'UniformOutput',false);
    logit_transform = @(x) log((x./100)./(1-(x./100)));
    logit_invtransform = @(x) (1./(1+exp(-x)))*100;
    if strcmp(inputs.errorModel,'additive')
    if and(~isempty(inputs.logit_indx), ismember(i, inputs.logit_indx))
        outputPred = arrayfun(@(x)(logit_invtransform(logit_transform(x.(output_i))...
            +randn(size(x.(output_i))).*(sigma_i)))...
            ,PI.output,'UniformOutput',false);
        [PI.outputPred(1:end).(output_i)] = outputPred{:,:};
        
    else
        outputPred = arrayfun(@(x)exp(log(x.(output_i))+...
            randn(size(x.(output_i))).*(sigma_i))...
            ,PI.output,'UniformOutput',false);
        [PI.outputPred(1:end).(output_i)] = outputPred{:,:};
 
    end
    else
        
         outputPred = arrayfun(@(x)((x.(output_i))+...
            randn(size(x.(output_i))).*(inputs.constantVar+sigma_i.*x.(output_i)))...
            ,PI.output,'UniformOutput',false);
        [PI.outputPred(1:end).(output_i)] = outputPred{:,:};
 
    end
    outputMean=arrayfun(@(x) mean(x.(output_i),'omitnan'),PI.outputPred,'UniformOutput',false)';
    outputMedian=arrayfun(@(x) median(x.(output_i),'omitnan'),PI.outputPred,'UniformOutput',false)';
    predLB=arrayfun(@(x)quantile(x.(output_i), 0.025), PI.outputPred,...
        'UniformOutput', false)';
    predUB=arrayfun(@(x)quantile(x.(output_i), 0.975), PI.outputPred,...
        'UniformOutput', false)';
    
    for j=1:length(outputMedian)
        table_ij=table([outputMean{j,:}; outputMedian{j,:}; outputLB{j,:}; outputUB{j,:}; predLB{j,:}; predUB{j,:}],...
            'RowNames', {'Mean' 'Median' 'LB' 'UB' 'Pred_LB' 'Pred_UB'});
        PI.CI(j).(output_i)=table_ij;
    end
    
    %% Add posterior prediction errors
    data = {PI.data(1:end).dataValue};                                      % Add data to outputPred array
    [PI.outputPred(1:end).dataValue] = data{:,:};
    
    dataTimeIndx = arrayfun(@(x) ismember(inputs.simTime, x.dataTime),...   % Match observed time points to simulated time points
        PI.data, 'UniformOutput',false);
    [PI.outputPred(1:end).dataTimeIndx] = dataTimeIndx{:,:};
    
    postPred = arrayfun(@(x) x.(output_i)(:,x.dataTimeIndx),...             % Select the posterior predictions corresponding to the observed time points
        PI.outputPred, 'UniformOutput', false);
    field_i = strjoin({output_i 'predError'}, '_');
    [PI.outputPred(1:end).(field_i)] = postPred{:,:};
    
    predError = arrayfun(@(x) x.dataValue(:,i)' -  x.(field_i),...          % Calculate absolute errors
        PI.outputPred, 'UniformOutput',false);
    [PI.outputPred(1:end).(field_i)] = predError{:,:};
    
    MSE = (arrayfun(@(x) sum((x.dataValue(:,i)' - mean(x.(output_i)(:,x.dataTimeIndx),...
        'omitnan')).^2),...
        PI.outputPred, 'UniformOutput', true));

    PI.MSE(:,i) = MSE';
end

return
