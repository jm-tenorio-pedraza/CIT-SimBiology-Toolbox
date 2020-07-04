function surv = getMedianPFS(PI,T, censor, treatments)
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
    
surv{i,2} = log(2)*sum(T(indx))/sum(~censor(indx+N*(i-1)))/30;

surv{:,3} = log(2)./exp(log(1./(surv{:,2}/log(2)))+1.96*sum(~censor(indx+N*(i-1))).^(-1/2));

surv{:,4} = log(2)./exp(log(1./(surv{:,2}/log(2)))-1.96*sum(~censor(indx+N*(i-1))).^(-1/2));
end
surv.Properties.VariableNames = {'Treatments', 'Median_PFS', 'LB', 'UB'};
return
