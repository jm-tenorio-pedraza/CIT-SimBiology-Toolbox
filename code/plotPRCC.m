function ranking = plotPRCC(PI,parameters,observables,varargin)

par=inputParser;
par.addParameter('n',1)
par.addParameter('time',1:1:PI.tspan(end))
par.addParameter('timeUnit', 'days')
par.addParameter('output', 'individual')
par.addParameter('kpi', 'mean')
par.addParameter('newFig', true)
par.addParameter('plotIndx', 1:length(observables))
par.addParameter('ncol', ceil(sqrt(length(observables))))
par.addParameter('nrow', ceil(length(observables)/ceil(sqrt(length(observables)))))

par.parse(varargin{:})
par=par.Results;
time = par.time;

%% Plot results
colors = linspecer(length(parameters));
styles = repmat({'-'; '--'; '-.'; ':'}, ceil(length(parameters)/4),1);
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
    if i==1 && par.newFig
        figure
    end
    ncol = par.ncol;
    nrow = par.nrow;
    plotIndx = par.plotIndx;
    subplot(nrow,ncol,plotIndx(i))
    hold on
        for j = 1:length(parameters)
            prcc_j= reshape(PI.prcc(:,j,i), [], length(PI.u));
            prcc_mean(:,j) = mean(prcc_j,2);
            h= plot(time, prcc_mean(:,j));
            
            line_j =  styles{j};
            set(h, 'Color', colors(j,:),'LineWidth', 2,'LineStyle', line_j)
        end
        plot(time, ones(1,length(time))*0.5, 'Linewidth', 2, 'Color', 'black');
        plot(time, ones(1,length(time))*(-0.5), 'Linewidth', 2, 'Color', 'black');
        if strcmp(par.kpi, 'mean')
        [~, prcc_mean_indx] = sort(mean(abs(prcc_mean)),'descend');
        elseif strcmp(par.kpi,'max')
        [~, prcc_mean_indx] = sort(max(abs(prcc_mean)),'descend');

        end
        ranking(:,i) = parameters(prcc_mean_indx);
        title(PI.observablesPlot(i), 'interpreter', 'tex')
        xlabel(strjoin({'Time [' par.timeUnit ']'}, ''))
        ylabel('Partial correlation coefficient')
        ylim([-1 1])
end
    ax = gca;
    ax.FontSize = 14;
end
ranking.Properties.VariableNames = observables;

legend(parameters)

return
