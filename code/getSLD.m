function PI = getSLD(PI,sigma, varargin)
inputs = inputParser;
inputs.addParameter('TV_0', PI.output(1).TV(:,1));
inputs.parse(varargin{:});
inputs = inputs.Results;

N = size(inputs.TV_0,1);

SLD = arrayfun(@(x)(x.TV./inputs.TV_0).^(1/3)*100-100, PI.output,...
    'UniformOutput', false);
[PI.output(1:end).SLD] = SLD{:,:};

TV_sigma = arrayfun(@(x)exp(log(x.TV)+randn(N,...
        size(x.TV,2)).*sigma(:,1)), PI.output, 'UniformOutput', false);
    [PI.output(1:end).TV_sigma] = TV_sigma{:,:};

CD8_sigma = arrayfun(@(x)exp(log(x.CD8)+randn(N,...
        size(x.CD8,2)).*sigma(:,2)), PI.output, 'UniformOutput', false);
    [PI.output(1:end).CD8_sigma] = CD8_sigma{:,:};

SLD_sigma = arrayfun(@(x)((x.TV_sigma./inputs.TV_0).^(1/3))*100-100,...
    PI.output, 'UniformOutput', false);
[PI.output(1:end).SLD_sigma] = SLD_sigma{:,:};
return

