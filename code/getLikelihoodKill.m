function likelihood= getLikelihoodKill(p,PI,model,varargin)
par=inputParser;
par.addParameter('free','E')
par.parse(varargin{:})
par=par.Results;
if strcmp(par.free,'E')
    E = num2cell(p(3:end-1));
    kill = [repelem({p(1)},12,1);repelem({p(2)},7,1)];
    [PI.data(1:end).E] = E{:,:};
else
    kill = num2cell(p(1:end-1));
end
[PI.data(1:end).kill] =kill{:,:};

sim = arrayfun(@(x)model(x.kill,x.kpro,x.E,x.T_0,x.dataTime),PI.data,'UniformOutput',false);
[PI.data(1:end).simValues] = sim{:,:};

likelihood = arrayfun(@(x)reshape(-(log(x.dataValue)-log(x.simValues)).^2/(2*p(end)^2)...
    -1/2*log(2*pi*p(end)^2),1,[]),PI.data, 'UniformOutput',false);

likelihood = sum([likelihood{:,:}]);

likelihood_ces = arrayfun(@(x)reshape(likelihood_censored(log(x.simValues(end)),...
    log(x.dataValue(end)),p(end),x.censored),1,[]),PI.data,'UniformOutput',false);
likelihood_ces=[likelihood_ces{:}];
likelihood = likelihood+sum(likelihood_ces);
return
