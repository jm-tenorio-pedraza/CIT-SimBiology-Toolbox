function PI = getCredibleIntervals(PI,outputs,params,H)
sigma=params(:,H.SigmaParams);

PI.CI=struct('Name',{PI.data(1:end).Name}','Group',{PI.data(1:end).Group}')';
for i=1:length(outputs)
    output_i=char(outputs(i));
    sigma_i=sigma(:,i);
    outputMedian=arrayfun(@(x) median(x.(output_i),'omitnan'),PI.output,'UniformOutput',false);
    outputLB=arrayfun(@(x)quantile(x.(output_i),0.025),PI.output,'UniformOutput',false);
    outputUB=arrayfun(@(x)quantile(x.(output_i),0.975),PI.output,'UniformOutput',false);
    predLB=arrayfun(@(x)quantile(exp(log(x.(output_i))+randn(size(x.(output_i))).*repmat(sigma_i,1,size(x.(output_i),2))),0.025)...
        ,PI.output,'UniformOutput',false);
    predUB=arrayfun(@(x)quantile(exp(log(x.(output_i))+randn(size(x.(output_i))).*repmat(sigma_i,1,size(x.(output_i),2))),0.975)...
        ,PI.output,'UniformOutput',false);
    for j=1:length(outputMedian)
        table_ij=table([outputMedian{j,:}; outputLB{j,:}; outputUB{j,:}; predLB{j,:}; predUB{j,:}],...
            'RowNames', {'Median' 'LB' 'UB' 'Pred_LB' 'Pred_UB'});
        PI.CI(j).(output_i)=table_ij;
    end
end
