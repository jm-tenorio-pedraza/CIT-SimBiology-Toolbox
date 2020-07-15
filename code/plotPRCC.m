function ranking = plotPRCC(PI,parameters,observables,varargin)

par=inputParser;
par.addParameter('n',1)
par.addParameter('time',1:1:PI.tspan(end))
par.addParameter('timeUnit', 'days')
par.addParameter('output', 'individual')

par.parse(varargin{:})
par=par.Results;
time = par.time;

%% Plot results
colors = linspecer(length(parameters));
styles = repmat({'-'; '--'}, ceil(length(parameters)/2),1);
ranking = table(repelem({'nan'},length(parameters),1));
prcc_mean = nan(length(time), length(parameters));
for i = 1:length(observables)

if strcmp(par.output, 'individual')
    figure
    ncol = ceil(sqrt(length(PI.u)));
    nrow = ceil(length(PI.u)/ncol);

    for k=1:length(PI.u)
        index = (k-1)*length(time)+1:(k-1)*length(time)+length(time);
        subplot(nrow, ncol, k)
        hold on
        h= plot(time, PI.prcc(index,:,i));
        for j = 1:length(parameters)
            line_j =  styles{j};
            set(h(j), 'Color', colors(j,:),'LineWidth', 2,'LineStyle', line_j)
        end
%         plot(time,repelem(0.2, size(time,1), size(time,2)), '--r', 'LineWidth', 2)
%         plot(time,repelem(-0.2, size(time,1), size(time,2)), '--r', 'LineWidth', 2)

        title(PI.data(k).Name, 'interpreter', 'none')
        xlabel(strjoin({'Time [' par.timeUnit ']'}, ''))
        ylabel('Partial correlation coefficient')
        
        ylim([-1 1])
        
    end
    legend(parameters)
else
    if i==1
        figure
    end
    ncol = ceil(sqrt(length(observables)));
    nrow = ceil(length(observables)/ncol);
    subplot(nrow,ncol,i)
    hold on
        for j = 1:length(parameters)
            prcc_j= reshape(PI.prcc(:,j,i), [], length(PI.u));
            prcc_mean(:,j) = mean(prcc_j,2);
            h= plot(time, prcc_mean(:,j));

            line_j =  styles{j};
            set(h, 'Color', colors(j,:),'LineWidth', 2,'LineStyle', line_j)
        end
        [~, prcc_mean_indx] = sort(abs(mean(prcc_mean)),'descend');
        ranking(:,i) = parameters(prcc_mean_indx);
        title(PI.observablesPlot(i), 'interpreter', 'tex')
        xlabel(strjoin({'Time [' par.timeUnit ']'}, ''))
        ylabel('Partial correlation coefficient')
        
        ylim([-1 1])
        
end
    
end
ranking.Properties.VariableNames = observables;

legend(parameters)

return
