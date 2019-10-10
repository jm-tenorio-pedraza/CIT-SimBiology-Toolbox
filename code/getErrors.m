function errors = getErrors(PI, sigma)

error=arrayfun(@(x)reshape((log(x.y_hat)-log(x.dataValue)).^2./...% squared residuals
    (2*(sigma.^2))+...% normalized by their variance
    (log(sigma*sqrt(pi*2))),1,[]),...% adding their estimated variance
    PI.data,'UniformOutput',false);

error=cellfun(@(x)x(~isnan(x)),error,'UniformOutput',false);
errors=[error{:,:}];
return