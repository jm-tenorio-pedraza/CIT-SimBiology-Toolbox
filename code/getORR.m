function [PI, ORR] = getORR(PI, treatments,varargin)
inputs=inputParser;
inputs.addParameter('TV_0', []);
inputs.addParameter('TumorField', 'Tumor');
inputs.addParameter('cutoff_value', PI.tspan(end))
inputs.addParameter('N', 1)
inputs.parse(varargin{:})
inputs=inputs.Results;
tumorField=inputs.TumorField;
cutoff_indx = find(PI.tspan<=inputs.cutoff_value);

if isempty(inputs.TV_0)
    TV_0=PI.output(1).(tumorField)(1:end,1);
else
    TV_0 = inputs.TV_0;
end
colors = linspecer(4);
colors_i = zeros(size(PI.output(1).SLD,1),3);
for i=1:length(treatments)
    SLD = 2*real((PI.output(i).Tumor(:,cutoff_indx)/(4/3*pi)).^(1/3));
% %     minTumorSize = 2*real(min((PI.output(i).Tumor(:,cutoff_indx)/(4/3*pi)).^(1/3),[],2));
    minTumorSize = 2*real(((PI.output(i).Tumor(:,1)/(4/3*pi)).^(1/3)));

    SLD_change = real(((SLD-minTumorSize)./minTumorSize)*100);
    pd = and(any(SLD_change>20,2),...
        any(SLD-minTumorSize>0.5,2));
    pr = and((and(any(PI.output(i).SLD(:,cutoff_indx)<=-30,2), ...
        ~any(PI.output(i).SLD(:,cutoff_indx)<-99,2))),~pd);
    cr = and((any(PI.output(i).SLD(:,cutoff_indx)<=-99,2)),~pd);
    sd = and(and(~pd,~pr), ~cr);
    
    SLD_sigma_rescaled = 2*real((PI.output(i).Tumor_sigma(:,cutoff_indx)/(4/3*pi)).^(1/3));
    SLD_sigma_change = real((SLD_sigma_rescaled)./minTumorSize)*100;
    pd_sigma =and(any(SLD_sigma_change>20,2),...
        any(SLD_sigma_change>0.5,2));
    pr_sigma = and((and(any(PI.output(i).SLD_sigma(:,cutoff_indx)<=-30,2), ...
        ~any(PI.output(i).SLD_sigma(:,cutoff_indx)<-99,2))),~pd);
    cr_sigma = and((any(PI.output(i).SLD_sigma(:,cutoff_indx)<-99,2)),~pd);
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

ORR2 = cell2mat(arrayfun(@(x)mean(x.Response), PI.output,'UniformOutput', false));
ORR2_sigma =  cell2mat(arrayfun(@(x)mean(x.Response_sigma), PI.output,'UniformOutput', false));
ORR2_SE = sqrt(ORR2_sigma.*(1-ORR2_sigma)/inputs.N);
ORR2_CI_UB = round(ORR2+1.96*ORR2_SE,3);
ORR2_CI_LB = round(ORR2-1.96*ORR2_SE,3);

ORR = table({'PD'; 'SD'; 'PR'; 'CR'}, 'VariableNames', {'Response'});

for i=1:length(treatments)
    indx = (i-1)*2+2;
    ORR{:,indx} = ORR2(i,:)';
    ORR{:,indx+1} = [cellfun(@num2str, num2cell(ORR2_CI_LB(i,:)), 'UniformOutput', false)',...
        cellfun(@num2str, num2cell(ORR2_CI_UB(i,:)), 'UniformOutput', false)'];
    ORR.Properties.VariableNames(indx) = treatments(i);
    ORR.Properties.VariableNames(indx+1) = {strjoin({treatments{i} 'CI'}, '_')};
end



return
