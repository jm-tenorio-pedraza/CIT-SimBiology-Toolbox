function plotPRCC(PI,parameters,observables,varargin)

par=inputParser;
par.addParameter('n',1)
par.addParameter('time',1:1:PI.tspan(end))
par.parse(varargin{:})
par=par.Results;
time = par.time;

%% Plot results
ncol = ceil(sqrt(length(PI.u)));
nrow = ceil(length(PI.u)/ncol);
colors = linspecer(length(parameters));
styles = repmat({'-'; '--'}, ceil(length(parameters)/2),1);
for i = 1:length(observables)
figure
    for k=1:length(PI.u)
        index = (k-1)*length(time)+1:(k-1)*length(time)+length(time);
        subplot(nrow, ncol, k)
        h= plot(time, PI.prcc(index,:,i));
        for j = 1:length(parameters)
            
            line_j =  styles{j};
            set(h(j), 'Color', colors(j,:),'LineWidth', 2,'LineStyle', line_j)
        end
        
        title(PI.data(k).Name)
        xlabel('Time [days]')
        ylabel('Partial correlation coefficient')
        
        ylim([-1 1])
        
    end
    legend(parameters)
end

legend(parameters)

return
