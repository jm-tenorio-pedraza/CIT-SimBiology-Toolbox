function [PI, ORR] = getORR(PI, treatments,varargin)
inputs=inputParser;
inputs.addParameter('TV_0', PI.output(1).TV(:,1))
inputs.parse(varargin{:})
inputs=inputs.Results;

ORR = table({'PD'; 'SD'; 'PR'; 'CR'},zeros(4,1));
colors = linspecer(4);
colors_i = zeros(size(PI.output(1).SLD,1),3);

for i=2:length(treatments)
    ORR{:,end+1} = zeros(4,1);
end
ORR.Properties.VariableNames = ['Response' treatments];
ORR = ORR;
for i=1:length(treatments)
    pd = and(any(PI.output(i).SLD>20,2),...
        any(PI.output(i).TV/(4/3*pi).^(1/3)...
        -inputs.TV_0./(4/3*pi).^(1/3)>0.5,2));
    pr = and((and(any(PI.output(i).SLD<=-30,2), ...
        ~any(PI.output(i).SLD<-99,2))),~pd);
    cr = and((any(PI.output(i).SLD<-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    ORR{1:4,i+1} = [mean(pd); mean(sd); mean(pr); mean(cr)];
   
    pd_sigma =and(any(PI.output(i).SLD_sigma>20,2),...
        any(PI.output(i).TV_sigma/(4/3*pi).^(1/3)...
        -inputs.TV_0./(4/3*pi).^(1/3)>0.5,2));
    pr_sigma = and((and(any(PI.output(i).SLD_sigma<=-30,2), ...
        ~any(PI.output(i).SLD_sigma<-99,2))),~pd);
    cr_sigma = and((any(PI.output(i).SLD_sigma<-99,2)),~pd);
    sd_sigma = and(and(~pd_sigma,~pr_sigma), ~cr_sigma);
%     ORR{1:4,i+1} = [mean(pd_sigma); mean(sd_sigma); mean(pr_sigma); mean(cr_sigma)];
    
    colors_i(pd,:) = repmat(colors(1,:),sum(pd),1);
    colors_i(sd,:) = repmat(colors(2,:),sum(sd),1);
    colors_i(pr,:) = repmat(colors(3,:),sum(pr),1);
    colors_i(cr,:) = repmat(colors(4,:),sum(cr),1);
    PI.output(i).Response = [pd sd pr cr];
    PI.output(i).Response_sigma = [pd_sigma sd_sigma pr_sigma cr_sigma];
    PI.output(i).colors = colors_i;

end
return
