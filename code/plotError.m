function plotError(sigma,PI,colIndx)

n_sim=size(PI.data,1);
n_col=ceil(sqrt(n_sim));
n_row=ceil(n_sim/n_col);
treatments=([PI.data(:).Group]);
if length(treatments)~=n_sim
    treatments = unique({PI.data(:).Group});
else
    treatments = unique(treatments);
end
treatment_colors=linspecer(length(treatments));
figure
for i=1:n_sim
    subplot(n_row,n_col,i)
    hold on
    if ~all(isnan(PI.data(i).dataValue(:,colIndx)))
    dat=plot(PI.data(i).dataTime, PI.data(i).dataValue(:,colIndx) - ...
        PI.data(i).y_hat(:,colIndx));
    else
        continue
    end
    zeroPlot = plot(PI.data(i).dataTime, zeros(size(PI.data(i).dataTime)),'-k');
    ub = plot(PI.data(i).dataTime, ones(size(PI.data(i).dataTime))*0.2,'--k');
    lb = plot(PI.data(i).dataTime, -ones(size(PI.data(i).dataTime))*0.2,'--k');

    title(PI.data(i).Name)
    
    col_i=treatment_colors(ismember(treatments,PI.data(i).Group),:);
    sim.Color=col_i;
    dat.LineStyle='none';
    dat.MarkerFaceColor=col_i;
    dat.MarkerEdgeColor=col_i;

    dat.Marker='d';
     zeroPlot.LineStyle = 'none';
    
    
    ylabel(PI.variableUnits(colIndx))

    
    % Change axis if there is large variability in output
%     if std(log10(PI.data(i).simOutput(:,colIndx)))>1
%         set(gca,'YScale','log')
%     end
       %set(gca,'YScale','log')
% ylim([-1 1])
end
return