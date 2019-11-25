function residuals= getResidualsKill(p,PI,model,varargin)
par=inputParser;
par.addParameter('free','E')
par.parse(varargin{:})
par=par.Results;
if strcmp(par.free,'E')
    E = num2cell(p(3:end));
    kill = [repelem({p(1)},12,1);repelem({p(2)},7,1)];
    [PI.data(1:end).E] = E{:,:};
else
    kill = num2cell(p(1:end));
end
[PI.data(1:end).kill] =kill{:,:};

sim = arrayfun(@(x)model(x.kill,x.kpro,x.E,x.T_0,x.dataTime),PI.data,'UniformOutput',false);
[PI.data(1:end).simValues] = sim{:,:};

residuals = arrayfun(@(x)reshape(((x.dataValue)-(x.simValues)),1,[]),PI.data, 'UniformOutput',false);

residuals = [residuals{:,:}];
return
