function [PI, response] = getPFS(PI, treatments, varargin)
inputs=inputParser;
inputs.addParameter('TV_0', [])
inputs.addParameter('cutoff_value', PI.tspan(end))
inputs.addParameter('TumorField', 'Tumor');
inputs.parse(varargin{:})
inputs=inputs.Results;

tumorField=inputs.TumorField;
cutoff_indx = find(PI.tspan<=inputs.cutoff_value);

if isempty(inputs.TV_0)
    TV_0=PI.output(1).(tumorField)(1:end,1);
else
    TV_0 = inputs.TV_0;
end

response = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1));

for i=1:length(treatments)
   SLD_rescaled = 2*real((PI.output(i).Tumor(:,cutoff_indx)/(4/3*pi)).^(1/3));
    minTumorSize = 2*real(min((PI.output(i).Tumor(:,cutoff_indx)/(4/3*pi)).^(1/3),[],2));
    SLD_change = real((((SLD_rescaled-minTumorSize)./minTumorSize)*100));
    pd = and(any(SLD_change>20,2),...
        any(SLD_rescaled...
        -minTumorSize>0.5,2));
    pr = and((and(any(PI.output(i).SLD(:,cutoff_indx)<=-30,2), ...
        ~any(PI.output(i).SLD<-99,2))),~pd);
    cr = and((any(PI.output(i).SLD(:,cutoff_indx)<-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    response{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
    PI.output(i).PFS = [pd sd pr cr];
%     PI.output(i).PFS_St
end
response.Properties.VariableNames = ['Response' treatments];
