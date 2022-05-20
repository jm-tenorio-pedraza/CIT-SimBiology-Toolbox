function plotPosteriorPredictions2(PI,outputs,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end
p=inputParser;
p.addParameter('simTime',PI.tspan);
p.addParameter('plot', 'data');
p.addParameter('indx', 1:length(PI.output));
p.addParameter('central', 'median')
p.addParameter('newFig', true)
p.addParameter('TimeUnit', 'days')
p.addParameter('interpreter', 'tex');
p.addParameter('color', 'group');
p.addParameter('YScale', 'log');
p.addParameter('indivData', false) 
p.addParameter('treatments', []) 

p.parse(varargin{:});
p=p.Results;
if isempty(p.treatments)
    treatments = unique({PI.output(1:end).Treatment},'stable');
else
    treatments=p.treatments;
end
cells = unique({PI.output(1:end).Cell},'stable');
response = unique({PI.output(1:end).Response});
nrow=length(treatments);
ncol = length(cells);
nresp = length(response);
treatments_labels = cellfun(@(x)strrep(x,'_',' + '),treatments,'UniformOutput',false);
for i=1:length(outputs)
    figure('Position',[10 10 1500 700])
    tiledlayout(nrow,ncol, 'TileSpacing','compact','Padding','compact')
    output_i = outputs{i};
    tileIndx = 1;
    % Determine consistent data ranges for all simulations
    observerIndx = arrayfun(@(x)ismember(x.observables,output_i),PI.IndivData,'UniformOutput',false);
    [PI.IndivData(1:end).ObserverIndx] = observerIndx{:,:};
    try
        minX = min(arrayfun(@(x) min(x.(output_i){'Pred_LB',:}), PI.CI, 'UniformOutput', true));
%         minX = min(arrayfun(@(x) min(x.dataValue(:,x.ObserverIndx)), PI.IndivData, 'UniformOutput', true));
%         maxX = max(arrayfun(@(x) max(x.dataValue(:,x.ObserverIndx)), PI.IndivData, 'UniformOutput', true));
        maxX = max(arrayfun(@(x) max(x.(output_i){'Pred_UB',:}), PI.CI, 'UniformOutput', true));
    catch
    end
    for j=1:nrow
        
        treat_j = treatments(j);
        for k=1:ncol
            
            cell_k = cells(k);
            for l=1:nresp
                resp_l = response(l);
                indx = and(and(ismember({PI.output(1:end).Treatment},treat_j),...
                    ismember({PI.output(1:end).Cell},cell_k)),ismember({PI.output(1:end).Response},resp_l));
                CI = {PI.CI(indx).(output_i)};
                CI_names = {PI.CI(indx).Name};
                colors = linspecer(nresp);
                if ~isempty(CI)
                    nexttile(tileIndx)
                    ax=gca;
                    hold on
                    grid on
                    for m = 1:length(CI)
                        ci = CI{m};
                        ci_data = [ci{'UB',:},ci{'LB',:}(end:-1:1)];
                        m_data = ci{'Median',:};
                        pi_data= [ci{'Pred_UB',:},ci{'Pred_LB',:}(end:-1:1)];
                        dataIndx = ismember({PI.IndivData(1:end).Name},CI_names{m});
                        % Determine whether nan's are observed to eliminate them
                        simTime = [1:1:length(m_data) length(m_data):-1:1];
                        ci_nan=or(isnan(ci_data),ci_data==0);
                        m_nan=or(isnan(m_data), m_data==0);
                        pi_nan=or(isnan(pi_data), pi_data==0);
                        
                        pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));       % Plot PIs
                        ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));           % Credible interval
                        sim=plot(simTime(~m_nan),m_data(~m_nan));
                        dat=arrayfun(@(x)plot(x.dataTime,x.dataValue(:,x.ObserverIndx)),...
                            PI.IndivData(dataIndx),'UniformOutput',false);
                        % Line specs
                        col_i=colors(l,:);
                        pi.FaceColor=col_i;
                        pi.EdgeColor=col_i;
                        pi.FaceAlpha=0.2;
                        pi.LineStyle='none';
                        
                        ci.FaceColor=col_i;
                        ci.EdgeColor=col_i;
                        ci.FaceAlpha=0.3;
                        ci.LineStyle='-';
                        ci.LineWidth=1;
                        
                        sim.Color=col_i;
                        
                        for w=1:length(dat)
                            dat{w}.Color =col_i;
                            dat{w}.LineStyle='none';
                            dat{w}.Marker='d';
                            dat{w}.MarkerFaceColor=col_i;
                            dat{w}.MarkerEdgeColor=col_i;
                        end
                    end
                else
                    
                end
            end
            % Axes specs
            if k==1
                ylabel(({treatments_labels{j}; PI.variableUnits{i}}),...
                    'interpreter', 'none','fontweight','bold','fontsize',12)
            else
                ax.YTickLabels = {};
            end
            if j == nrow
                xlabel(strjoin({'Time [' p.TimeUnit ']'},''),'fontweight','bold','fontsize',12)
                lagTime = ceil(length(m_data)/10);
                simTimeLabel = cellfun(@(x)num2str(x),num2cell(0:lagTime:length(m_data)),'UniformOutput',false);
                ax.XTick=0:lagTime:length(m_data);
                ax.XTickLabels = simTimeLabel;
            else
                ax.XTickLabels = {};
            end
            if j==1
                title(cells(k),'interpreter',p.interpreter,'fontweight','bold','fontsize',16)
            end
            if maxX<1
                ylim([(minX), (maxX)])
            else
                ylim([floor(minX), ceil(maxX)])
            end
            set(ax, 'YScale', p.YScale)
                        tileIndx = tileIndx +1;

        end
    end
end



