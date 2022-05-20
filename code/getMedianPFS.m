function surv = getMedianPFS(PI,treatments,varargin)
inputs=inputParser;
inputs.addParameter('Method','Brookmeyer')
inputs.addParameter('survivalType','PFS')

inputs.parse(varargin{:})
inputs=inputs.Results;
N = size(PI.output(1).Response,1);
indx = 1:N;
N_p = size(treatments);
if N_p(1)>N_p(2)
    N_p=N_p(1);
else
    N_p = N_p(2);
    treatments=treatments';
end
surv = table(treatments);
for i=1:N_p
    
if strcmp(inputs.Method, 'Brookmeyer')
    if strcmp(inputs.survivalType,'PFS')
        surv{i,2} = log(2)*sum(PI.output(i).T)/sum(~PI.output(i).Censor)/30;

        surv{i,3} = log(2)./exp(log(sum(~PI.output(i).Censor)/sum(PI.output(i).T))+1.96*sum(~PI.output(i).Censor).^(-1/2))/30;
        
        surv{i,4} = log(2)./exp(log(sum(~PI.output(i).Censor)/sum(PI.output(i).T))-1.96*sum(~PI.output(i).Censor).^(-1/2))/30;

    else
        surv{i,2} = log(2)*sum(PI.output(i).OS_T)/sum(~PI.output(i).OS_Censor)/30;

        surv{i,3} = log(2)./exp(log(sum(~PI.output(i).OS_Censor)/sum(PI.output(i).OS_T))+1.96*sum(~PI.output(i).OS_Censor).^(-1/2))/30;
        
        surv{i,4} = log(2)./exp(log(sum(~PI.output(i).OS_Censor)/sum(PI.output(i).OS_T))-1.96*sum(~PI.output(i).OS_Censor).^(-1/2))/30;
    end
end
end
surv.Properties.VariableNames = {'Treatments', 'Median_PFS', 'LB', 'UB'};
return
