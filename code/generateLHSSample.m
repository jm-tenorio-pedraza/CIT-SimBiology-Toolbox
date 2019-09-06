function params = generateLHSSample(parameters, meanValues,nsamples, varargin)

p=inputParser;
p.addParameter('Distribution','Normal');
p.addParameter('Sigma',repelem(1, 1,length(meanValues)));

p.parse(varargin{:});
p=p.Results;

n_p = length(parameters);

lhssample = lhsdesign(nsamples, n_p);
if strcmp(p.Distribution, 'Normal')
    params = norminv(lhssample, repmat(meanValues, nsamples,1),repmat(p.Sigma,nsamples,1));
elseif strcmp(p.Distribution, 'Uniform')
    params = unifinv(lhssample, repmat(meanValues*0.9,nsamples, 1), repmat(meanValues*1.1, nsamples,1));
end

return




