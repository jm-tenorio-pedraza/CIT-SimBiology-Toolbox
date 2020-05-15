function PI = getCredibleIntervals(PI,outputs,params,H,varargin)
inputs = inputParser;
inputs.addParameter('logit_indx', [])

inputs.parse(varargin{:})
inputs = inputs.Results;
sigma=params(:,setdiff(H.SigmaParams, [H.IndividualParams.OmegaIndex H.CellParams.OmegaIndex]));
PI.CI=struct('Name',{PI.data(1:end).Name}','Group',{PI.data(1:end).Group}')';
for i=1:length(outputs)
    output_i=char(outputs(i));
    sigma_i=sigma(:,i);
        outputMean=arrayfun(@(x) mean(x.(output_i),'omitnan'),PI.output,'UniformOutput',false);

    outputMedian=arrayfun(@(x) median(x.(output_i),'omitnan'),PI.output,'UniformOutput',false);
    outputLB=arrayfun(@(x)quantile(x.(output_i),0.025),PI.output,'UniformOutput',false);
    outputUB=arrayfun(@(x)quantile(x.(output_i),0.975),PI.output,'UniformOutput',false);
    logit_transform = @(x) log((x./100)./(1-(x./100)));
    logit_invtransform = @(x) (1./(1+exp(-x)))*100;
    if and(~isempty(inputs.logit_indx), ismember(i, inputs.logit_indx))
        predLB=arrayfun(@(x)quantile(logit_invtransform(logit_transform(x.(output_i))+randn(size(x.(output_i))).*repmat(sigma_i,1,size(x.(output_i),2))),0.025)...
            ,PI.output,'UniformOutput',false);
        predUB=arrayfun(@(x)quantile(logit_invtransform(logit_transform(x.(output_i))+randn(size(x.(output_i))).*repmat(sigma_i,1,size(x.(output_i),2))),0.975)...
            ,PI.output,'UniformOutput',false);
    else
        predLB=arrayfun(@(x)quantile(exp(log(x.(output_i))+randn(size(x.(output_i))).*repmat(sigma_i,1,size(x.(output_i),2))),0.025)...
            ,PI.output,'UniformOutput',false);
        predUB=arrayfun(@(x)quantile(exp(log(x.(output_i))+randn(size(x.(output_i))).*repmat(sigma_i,1,size(x.(output_i),2))),0.975)...
            ,PI.output,'UniformOutput',false);
    end
    
    for j=1:length(outputMedian)
        table_ij=table([outputMean{j,:}; outputMedian{j,:}; outputLB{j,:}; outputUB{j,:}; predLB{j,:}; predUB{j,:}],...
            'RowNames', {'Mean' 'Median' 'LB' 'UB' 'Pred_LB' 'Pred_UB'});
        PI.CI(j).(output_i)=table_ij;
    end
end
