function PI = getSLD(PI,sigma, varargin)
inputs = inputParser;
inputs.addParameter('TV_0', []);
inputs.addParameter('TumorField', 'Tumor');
inputs.parse(varargin{:});
inputs = inputs.Results;
tumorField=inputs.TumorField;
if isempty(inputs.TV_0)
    TV_0=PI.output(1).(tumorField)(1:end,1);
else
    TV_0 = inputs.TV_0;
end
N = size(inputs.TV_0,1);

SLD = arrayfun(@(x)real(((((x.(tumorField))./TV_0)).^(1/3))*100-100), PI.output,...
    'UniformOutput', false);
[PI.output(1:end).SLD] = SLD{:,:};

Tumor_sigma = arrayfun(@(x)real(exp(log(x.(tumorField))+randn(N,...
    size(x.(tumorField),2)).*sigma(:,1))), PI.output, 'UniformOutput', false);
[PI.output(1:end).Tumor_sigma] = Tumor_sigma{:,:};
try
    CD8_sigma = arrayfun(@(x)exp(log(x.CD8)+randn(N,...
        size(x.CD8,2)).*sigma(:,2)), PI.output, 'UniformOutput', false);
    [PI.output(1:end).CD8_sigma] = CD8_sigma{:,:};
catch
end
SLD_sigma = arrayfun(@(x)real((((((x.Tumor_sigma)./TV_0)).^(1/3))*100-100)),...
    PI.output, 'UniformOutput', false);
[PI.output(1:end).SLD_sigma] = SLD_sigma{:,:};
return

