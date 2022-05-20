function [params, lhssample] = generateLHSSample2(PI,parameters,nsamples, varargin)
p=inputParser;
p.addParameter('Distribution','Uniform');
p.addParameter('Variation',0.9);
p.addParameter('LB',[]);
p.addParameter('UB', []);
p.parse(varargin{:});
p=p.Results;

popParamsIndx = [PI.H.PopulationParams];
n_p = length(popParamsIndx);
if size(parameters,1)>size(parameters,2)
    parameters = parameters';
end
lhssample = lhsdesign(nsamples, n_p);
if isempty(p.LB)
    variation = repelem(p.Variation, 1, length(parameters));

    lb = parameters.*variation;
    ub = parameters./variation;
else
    lb = par.LB;
    ub = par.UB;
end
params = lb+lhssample.*(ub-lb);

return




