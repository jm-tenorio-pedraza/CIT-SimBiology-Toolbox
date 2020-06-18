function [PI, response] = getPFS(PI, treatments, varargin)
inputs=inputParser;
inputs.addParameter('TV_0', PI.output(1).TV(:,1))
inputs.addParameter('cutoff_value', PI.tspan(end))

inputs.parse(varargin{:})
inputs=inputs.Results;

cutoff_indx = find(PI.tspan==inputs.cutoff_value);

response = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1));

for i=1:length(treatments)
    pd = and(any(PI.output(i).SLD(:,1:cutoff_indx)>20,2),...
         any(PI.output(i).TV(:,1:cutoff_indx)/(4/3*pi).^(1/3)...
        -inputs.TV_0./(4/3*pi).^(1/3)>0.5,2));
    pr = and((and((PI.output(i).SLD(:,cutoff_indx)<=-30),...
        ~(PI.output(i).SLD(:,cutoff_indx)<-99))),~pd);
    cr = and(((PI.output(i).SLD(:,cutoff_indx)<-99)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
    PI.output(i).PFS = [pd sd pr cr];
%     PI.output(i).PFS_St
end
response.Properties.VariableNames = ['Response' treatments];
