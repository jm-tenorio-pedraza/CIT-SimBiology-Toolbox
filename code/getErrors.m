function errors = getErrors(PI, sigma,varargin)
p=inputParser;
p.addParameter('log',true)
p.parse(varargin{:})
p=p.Results;
if p.log
error=arrayfun(@(x)reshape((log(x.y_hat)-log(x.dataValue)).^2./...% squared residuals
    (2*(sigma.^2))+...% normalized by their variance
    (log(sigma*sqrt(pi*2))),1,[]),...% adding their estimated variance
    PI.data,'UniformOutput',false);
else
    error=arrayfun(@(x)reshape(((x.y_hat)-(x.dataValue)).^2./...% squared residuals
    (2*(sigma.^2))+...% normalized by their variance
    (log(sigma*sqrt(pi*2))),1,[]),...% adding their estimated variance
    PI.data,'UniformOutput',false);

end
error=cellfun(@(x)x(~isnan(x)),error,'UniformOutput',false);
errors=[error{:,:}];
return