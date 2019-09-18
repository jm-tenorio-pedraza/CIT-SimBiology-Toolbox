function var = createVariant(PI,H,variantname,varargin)

p=inputParser;
p.parse(varargin{:});
p=p.Results;

var = sbiovariant(variantname);

parameters = {PI.par(H.PopulationParams).name};

parameters([H.IndividualParams(:).EtaIndex]) = {H.IndividualParams(:).name};

values = num2cell([PI.par(H.PopulationParams).finalValue]);
content = [];
[content(1:length(parameters)).params] = parameters{:,:};
[content(1:end).values] = values{:,:};

content = arrayfun(@(x) {'parameter', x.params, 'value', x.values},...
    content, 'UniformOutput',false);

var.Content = content;
return


