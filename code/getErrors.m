function errors = getErrors(PI,varargin)
p=inputParser;
p.addParameter('log',true)
p.addParameter('indivData',false)
p.parse(varargin{:})
p=p.Results;
normIndx=ismember(PI.varDist,'Normal');
betaIndx=ismember(PI.varDist,'Beta');
if p.indivData
    if p.log
        errors=arrayfun(@(x)reshape((log(x.y_hat(:,normIndx))-log(x.dataValue(:,normIndx))).^2./...% squared residuals
            (2*(x.sigma(:,normIndx).^2))+...% normalized by their variance
            (log(x.sigma(:,normIndx)*sqrt(pi*2))),1,[]),...% adding their estimated variance
            PI.IndivData,'UniformOutput',false);
       
    else
        errors=arrayfun(@(x)reshape(((x.y_hat(:,normIndx))-(x.dataValue(:,normIndx))).^2./...% squared residuals
            (2*(x.sigma(:,normIndx).^2))+...% normalized by their variance
            (log(x.sigma(:,normIndx)*sqrt(pi*2))),1,[]),...% adding their estimated variance
            PI.IndivData,'UniformOutput',false);
    end
    if sum(betaIndx)>0
        alpha= arrayfun(@(x)((1-x.y_hat(:,betaIndx)/100)./(x.sigma(:,betaIndx).^2)-1./x.y_hat(:,betaIndx)/100).*(x.y_hat(:,betaIndx).^2),PI.IndivData,'UniformOutput', false);
        [PI.IndivData(1:end).('alpha')] = alpha{:,:};
        beta=arrayfun(@(x)x.alpha.*(1./(x.y_hat(:,betaIndx)/100)-1), PI.IndivData,'UniformOutput',false);
        [PI.IndivData(1:end).('beta')] = beta{:,:};

        betaloglik=arrayfun(@(x)reshape(log(betapdf(x.dataValue(:,betaIndx),...
            x.alpha, x.beta)),1,[]),...% adding their estimated variance
            PI.IndivData,'UniformOutput',false);
        errors=[errors; betaloglik];
    end
else
    if p.log
        errors=arrayfun(@(x)reshape((log(x.y_hat(:,normIndx))-log(x.dataValue(:,normIndx))).^2./...% squared residuals
            (2*(x.sigma(:,normIndx).^2))+...% normalized by their variance
            (log(x.sigma(:,normIndx)*sqrt(pi*2))),1,[]),...% adding their estimated variance
            PI.data,'UniformOutput',false);
    else
        errors=arrayfun(@(x)reshape(((x.y_hat(:,normIndx))-(x.dataValue(:,normIndx))).^2./...% squared residuals
            (2*(x.sigma(:,normIndx).^2))+...% normalized by their variance
            (log(x.sigma(:,normIndx)*sqrt(pi*2))),1,[]),...% adding their estimated variance
            PI.data,'UniformOutput',false);
    end
    if sum(betaIndx)>0
        alpha= arrayfun(@(x)(((1-x.y_hat(:,betaIndx)/100)./x.sigma(:,betaIndx).^2)-1./x.y_hat(:,betaIndx)).*x.y_hat(:,betaIndx).^2,PI.data,'UniformOutput', false);
        [PI.data(1:length(simOutput)).('alpha')] = alpha{:,:};
        beta=arrayfun(@(x)x.alpha.*(1./(x.y_hat(:,betaIndx)/100)-1), PI.data,'UniformOutput',false);
        [PI.data(1:length(simOutput)).('beta')] = beta{:,:};

        betaloglik=arrayfun(@(x)reshape(log(betapdf(x.dataValue(:,betaIndx),...
            x.alpha, x.beta)),1,[]),...% adding their estimated variance
            PI.data,'UniformOutput',false);
        errors=[errors; betaloglik];
    end
end
errors=cellfun(@(x)x(~isnan(x)),errors,'UniformOutput',false);
errors=[errors{:,:}];
return