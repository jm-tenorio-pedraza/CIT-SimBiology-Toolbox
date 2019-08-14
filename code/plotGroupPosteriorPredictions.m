function plotGroupPosteriorPredictions(PI,outputs,params,H,varargin)
if nargin<2
    error('GWMCMC:toofewinputs','AMCMC requires atleast 2 inputs.')
end

p=inputParser;
p.addParameter('simTime',0:1:100);
p.parse(varargin{:});
p=p.Results;

treatments=unique({PI.CI.Group});
treatment_colors=linspecer(length(treatments));

ncol=length(treatments);
nrow=length(outputs);

fig_indx=1;
figure


for i=1:length(outputs)
    var_i=char(outputs(i));

    for j=1:length(treatments)
        subplot(nrow,ncol,fig_indx)
        hold on
        % Obtain treatment indexes
        treatment_bool=ismember({PI.CI.Group}, treatments(j));
        
        % Extract output from corresponding treatment
        output_j = arrayfun(@(x) x.(var_i),PI.output(treatment_bool),'UniformOutput',false);
        output_j =cell2mat(output_j);
        
        % Calculate credible interval
        ci_ub=quantile(output_j,0.975);
        ci_lb=quantile(output_j,0.025);
        
        % Calculate prediction interval
        sigma=repmat(params(:,H.SigmaParams(i)),sum(treatment_bool),size(output_j,2));

        pi_ub=quantile(exp(log(output_j)+randn(size(output_j)).*sigma),0.975);
        pi_lb=quantile(exp(log(output_j)+randn(size(output_j)).*sigma),0.025);
        
        % Create datasets for plotting
        ci_data=[ci_ub, ci_lb(end:-1:1)];
        pi_data=[pi_ub, pi_lb(:,end:-1:1)];
        
        m_data=median(output_j);
        simTime=[p.simTime p.simTime(end:-1:1)];
        
        % Removing nan and zeros
        pi_nan=or(isnan(pi_data), pi_data==0);
        ci_nan=or(isnan(ci_data),ci_data==0);
        m_nan=or(isnan(m_data), m_data==0);
        
        % Plotting prediction interval
        pi=patch('XData', simTime(~pi_nan),'YData',pi_data(~pi_nan));

        % Plotting credible interval
        ci=patch('XData', simTime(~ci_nan),'YData',ci_data(~ci_nan));

        % Plotting simulation median
        sim=plot(simTime(~m_nan),m_data(~m_nan));
         
        
        col_i=treatment_colors(j,:);
        
        % Plotting data
        treatment_indx=find(ismember({PI.data.Group}, treatments(j)));
        for k=treatment_indx
            dat=plot(PI.data(k).dataTime, PI.data(k).dataValue(:,i));
            dat.LineStyle='none';
            dat.MarkerFaceColor=col_i;
            dat.MarkerEdgeColor=col_i;
            dat.Marker='d';
        end

        % Specifications of Prediction intervals
        pi.FaceColor=col_i;
        pi.EdgeColor=col_i;
        pi.FaceAlpha=0.2;
        pi.LineStyle='none';
        
        % Specifications of credible intervals
        ci.FaceColor=col_i;
        ci.EdgeColor=col_i;
        ci.FaceAlpha=0.4;
        ci.LineStyle='--';
        ci.LineWidth=1;
        
        % Specifications of median response
        sim.Color=col_i;
        sim.LineWidth= 2;
        
        % Axes labels
        title(treatments(j),'interpreter','none')
        ylabel(PI.variableUnits(i))
        xlabel('Time [days]')
        grid on
        hold on
        if max(m_data)/min(m_data)>100
            set(gca,'YScale','log')
        else
        end
        fig_indx=fig_indx+1;
        
        % Legends
        % Fraction of observations outside the credible regions
        n_obs=sum(arrayfun(@(x)length(x.dataValue(~isnan(x.dataValue(:,i)),i)),PI.data(treatment_indx)));
        err_c=sum(arrayfun(@(x) sum(or(x.dataValue(:,i)'<ci_lb(ismember(p.simTime,x.dataTime)),...
            x.dataValue(:,i)'>ci_ub(ismember(p.simTime,x.dataTime)))), PI.data(treatment_indx)));
        err_p=sum(arrayfun(@(x) sum(or(x.dataValue(:,i)'<pi_lb(ismember(p.simTime,x.dataTime)),...
            x.dataValue(:,i)'>pi_ub(ismember(p.simTime,x.dataTime)))), PI.data(treatment_indx)));
        text(100,1,strjoin({'% OOCI:' num2str(err_c/n_obs*100)},' '))
        text(100,2,strjoin({'% OOPI:' num2str(err_p/n_obs*100)},' '))
                

    end
end
