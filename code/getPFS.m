function [PI, response] = getPFS(PI, groups, cutoff_value)
cutoff_indx = find(PI.tspan==cutoff_value);

response = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1));

for i=1:length(groups)
    pd = any(PI.output(i).SLD(:,1:cutoff_indx)>20,2);
    pr = and((and((PI.output(i).SLD(:,cutoff_indx)<=-30), ~(PI.output(i).SLD(:,cutoff_indx)<-99))),~pd);
    cr = and(((PI.output(i).SLD(:,cutoff_indx)<-99)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
    PI.output(i).PFS = [pd sd pr cr];
    PI.output(i).PFS_St
end
response.Properties.VariableNames = ['Response' groups];
