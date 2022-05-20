function ranking = plotPRCC2(PI,parameters,observables,varargin)

par=inputParser;
par.addParameter('n',1)
par.addParameter('time',1:1:PI.tspan(end))
par.addParameter('timeUnit', 'days')
par.addParameter('output', 'individual')
par.addParameter('kpi', 'mean')
par.addParameter('newFig', true)
par.addParameter('plotIndx', 1:length(observables))
par.addParameter('nfig',1)
par.addParameter('nrow', 2)
par.addParameter('ncol', 3)
par.addParameter('nParams', 5)
par.parse(varargin{:})
par=par.Results;
time = par.time;

%% Plot results
colors = linspecer(par.nParams);
styles = repmat({'-'; '--'; '-.'; ':'}, ceil(length(parameters)/4),1);
ranking = table(repelem({'nan'},length(parameters),1));
prcc_mean = nan(length(time), length(parameters));
plotIndx = 1;

for i = 1:length(observables)
    if i==1 && par.newFig
        figure
    end
    if plotIndx==7
        plotIndx=1;
        figure
    end
    ncol = par.ncol;
    nrow = par.nrow;
    subplot(nrow,ncol,plotIndx)
    hold on
    for j = 1:length(parameters)
        prcc_j= reshape(PI.prcc(:,j,i), [], length(PI.u));
        prcc_mean(:,j) = mean(prcc_j,2);
    end
    if strcmp(par.kpi, 'mean')
        [~, prcc_mean_indx] = sort(mean(abs(prcc_mean)),'descend');
    elseif strcmp(par.kpi,'max')
        [~, prcc_mean_indx] = sort(max(abs(prcc_mean)),'descend');
    else
        [~, prcc_mean_indx] = sort(sum(abs(prcc_mean)),'descend');
        
    end
    ranking(:,i) = parameters(prcc_mean_indx);
    for j = 1:par.nParams
        paramIndx = prcc_mean_indx(j);
        
        h= plot(time, prcc_mean(:,paramIndx));
        line_j =  styles{j};
        set(h, 'Color', colors(j,:),'LineWidth', 2,'LineStyle', line_j)
    end
    plot(time, ones(1,length(time))*0.5, 'Linewidth', 2, 'Color', 'black');
    plot(time, ones(1,length(time))*(-0.5), 'Linewidth', 2, 'Color', 'black');
    
    h=gca;
    h= h.Children(end:-1:1);
    legend(h((1:par.nParams)),parameters(prcc_mean_indx(1:par.nParams)))
    
    title(PI.observablesPlot(i), 'interpreter', 'tex')
    xlabel(strjoin({'Time [' par.timeUnit ']'}, ''))
    ylabel('Partial correlation coefficient')
    ylim([-1 1])
    ax = gca;
    ax.FontSize = 14;
    plotIndx=plotIndx+1;
end
ranking.Properties.VariableNames = observables;


return
