function PI = getSurvivalTime(PI,groups,varargin)
inputs=inputParser;
inputs.addParameter('TV_0', []);
inputs.addParameter('TumorField', 'Tumor');
inputs.addParameter('cutoff_value', PI.tspan(end))
inputs.addParameter('survivalType', 'PFS')
inputs.addParameter('CC', [])

inputs.addParameter('N',1)
inputs.parse(varargin{:})
inputs=inputs.Results;
tumorField=inputs.TumorField;

if isempty(inputs.TV_0)
    TV_0=PI.output(1).(tumorField)(1:end,1);
else
    TV_0 = inputs.TV_0;
end
N=inputs.N;
cutoff_indx = find(PI.tspan<=inputs.cutoff_value);
for i=1:length(groups)
    if strcmp(inputs.survivalType,'PFS')
        SLD_rescaled = 2*real((PI.output(i).Tumor(:,cutoff_indx)/(4/3*pi)).^(1/3));
%         minTumorSize = 2*real(min((PI.output(i).Tumor(:,cutoff_indx)/(4/3*pi)).^(1/3),[],2));
        minTumorSize = 2*real(((PI.output(i).Tumor(:,1)/(4/3*pi)).^(1/3)));

        SLD_change = real(((SLD_rescaled-minTumorSize)./minTumorSize)*100);
%        pd = and(any(SLD_change>20,2),...
%              any(SLD_rescaled...
%             -minTumorSize>0.5,2));
        SLD = [];
        SLD_change = mat2cell(SLD_change,ones(size(SLD_change,1),1));

        minTV = num2cell(minTumorSize);
        [SLD(1:size(SLD_change,1)).SLD_change] = SLD_change{:,:};
                [SLD(1:size(SLD_change,1)).minTV] = minTV{:,:};

        pd = arrayfun(@(x)detectPD(x.SLD_change,PI.tspan(cutoff_indx),'minTV',x.minTV),...
            SLD,'UniformOutput',true)';
        T=pd;
        survival = PI.tspan(cutoff_indx)<=pd;
        atRisk = sum(survival);
        deaths = [0 diff(sum(survival))*(-1)];
        atRisk = atRisk+deaths;
        S_t = cumprod(1-deaths./atRisk);
        V_t = S_t.^2.*cumsum(deaths./(atRisk.*(atRisk-deaths)));
        censor = T==inputs.cutoff_value;
        PI.output(i).T = T;
        PI.output(i).S_t = S_t;
        PI.output(i).V_t = V_t;
        PI.output(i).Censor = censor;
    else
        pd = rowfun(@(x)detectPD(x,PI.tspan(cutoff_indx),'type','OS','CC',inputs.CC),...
            table(PI.output(i).(tumorField)(:,cutoff_indx)));
        T=pd{:,1};
        survival = PI.tspan(cutoff_indx)<=pd{:,1};
        atRisk = sum(survival);
        deaths = [0 diff(sum(survival))*(-1)];
        atRisk = atRisk+deaths;
        S_t = cumprod(1-deaths./atRisk);
        V_t = S_t.^2.*cumsum(deaths./(atRisk.*(atRisk-deaths)));
        censor = T==inputs.cutoff_value;
        PI.output(i).OS_T = T;
        PI.output(i).OS_S_t = S_t;
        PI.output(i).OS_V_t = V_t;
        PI.output(i).OS_Censor = censor;
    end
end
