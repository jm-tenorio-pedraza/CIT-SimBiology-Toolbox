function ranking = plotPRCC3(PI,parameters,observables,varargin)

par=inputParser;
par.addParameter('n',1)
par.addParameter('time',1:1:PI.tspan(end))
par.addParameter('timeUnit', 'days')
par.addParameter('output', 'individual')
par.addParameter('kpi', 'mean')
par.addParameter('newFig', true)
par.addParameter('plotIndx', 1:length(observables))
par.addParameter('plotTreatments', true)
par.addParameter('nfig',1)
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
nsim=length(PI.data);
nrow=nsim;
ncol = par.ncol;
group = [PI.data(1:end).Group];
group=cellfun(@(x)regexprep(x,'MOC1_',''),group,'UniformOutput',false);
for i = 1:length(observables)
    for k=1:nsim
        simIndx = (k-1)*length(time)+1;
        if plotIndx>ncol 
            plotIndx=1;
            figure('Position',[10 10 1500 700])
            tiledlayout(nrow,ncol, 'TileSpacing','compact','Padding','compact')
        end
        subplotIndx = plotIndx+(ncol)*(k-1);
        if  par.newFig && k*i==1
            figure('Position',[10 10 1500 700])
            tiledlayout(nrow,ncol, 'TileSpacing','compact','Padding','compact')

        end
%         subplot(nrow,ncol,subplotIndx)
        nexttile(subplotIndx)
        hold on
        for j = 1:length(parameters)
            prcc_j= (PI.prcc(simIndx:(simIndx+length(time)-1),j,i));
            prcc_mean(:,j) = (prcc_j);
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
        hChildren= h.Children(end:-1:1);
        legend(hChildren((1:par.nParams)),parameters(prcc_mean_indx(1:par.nParams)),'Location','bestoutside')
        if subplotIndx<=ncol
        title(PI.observablesPlot(i), 'interpreter', 'tex')
        end
        if k==nsim
        xlabel(strjoin({'Time [' par.timeUnit ']'}, ''))
        else
                    h.XTickLabel ={};

        end
        
        if ismember(subplotIndx, 1:(ncol):(nrow*ncol))
            ylabel(group(k),'interpreter','none')
        else
            h.YTickLabel ={};
        end
       
        ylim([-1 1])
        h.FontSize = 12;
        h.LineWidth=2;
    end
    plotIndx=plotIndx+1;

end
ranking.Properties.VariableNames = observables;


return
