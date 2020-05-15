function errors = getErrors(PI, sigma,varargin)
p=inputParser;
p.addParameter('log',true)
p.parse(varargin{:})
p=p.Results;
if p.log
errors=arrayfun(@(x)reshape((log(x.y_hat)-log(x.dataValue)).^2./...% squared residuals
    (2*(sigma.^2))+...% normalized by their variance
    (log(sigma*sqrt(pi*2))),1,[]),...% adding their estimated variance
    PI.data,'UniformOutput',false);
else
    errors=arrayfun(@(x)reshape(((x.y_hat)-(x.dataValue)).^2./...% squared residuals
    (2*(sigma.^2))+...% normalized by their variance
    (log(sigma*sqrt(pi*2))),1,[]),...% adding their estimated variance
    PI.data,'UniformOutput',false);

end
errors=cellfun(@(x)x(~isnan(x)),errors,'UniformOutput',false);
errors=[errors{:,:}];
return