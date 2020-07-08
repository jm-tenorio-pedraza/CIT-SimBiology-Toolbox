function [PI,T, censor] = getSurvivalTime(PI,groups,varargin)
inputs=inputParser;
inputs.addParameter('TV_0', PI.output(1).TV(:,1))
inputs.addParameter('cutoff_value', PI.tspan(end))
inputs.addParameter('N', size(PI.output(1).TV,1))

inputs.parse(varargin{:})
inputs=inputs.Results;

N=inputs.N;
cutoff_indx = find(PI.tspan==inputs.cutoff_value);
for i=1:length(groups)
    rowindx = (i-1)*N;
    pd = rowfun(@(x)detectPD(x,PI.tspan(1:cutoff_indx)),table(PI.output(i).SLD(:,1:cutoff_indx)));
    T(rowindx+1:rowindx+N,1)=pd{:,1};
    survival = PI.tspan(1:cutoff_indx)<=pd{:,1};
    atRisk = sum(survival);
    deaths = [0 diff(sum(survival))*(-1)];
    atRisk = atRisk+deaths;
    S_t = cumprod(1-deaths./atRisk);
    V_t = S_t.^2.*cumsum(deaths./(atRisk.*(atRisk-deaths)));
    
    PI.output(i).T = T;
    PI.output(i).S_t = S_t;
    PI.output(i).V_t = V_t;
end
censor = T==inputs.cutoff_value;