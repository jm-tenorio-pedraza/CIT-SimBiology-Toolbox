function plotSimOutput(PI,colIndx)

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
    sim=plot(PI.data(i).simTime, PI.data(i).simOutput(:,colIndx));
    dat=plot(PI.data(i).dataTime, PI.data(i).dataValue(:,colIndx));
try
    X = [PI.data(i).simTime; PI.data(i).simTime(end:-1:1)];
    Y = [PI.data(i).ub(:,colIndx); PI.data(i).lb(end:-1:1,colIndx)];
    nanindx = isnan(Y);
    error = patch('XData', X(~nanindx),...
        'YData',Y(~nanindx));
   
catch
    error = [];
end
    title(PI.data(i).Name)
    
    col_i=treatment_colors(ismember(treatments,PI.data(i).Group),:);
    sim.Color=col_i;
    dat.LineStyle='none';
    dat.MarkerFaceColor=col_i;
    dat.MarkerEdgeColor=col_i;

    dat.Marker='d';
     error.LineStyle = 'none';
    error.FaceColor = col_i;
    error.FaceAlpha = 0.2;
    
    ylabel(PI.variableUnits(colIndx))

    
    % Change axis if there is large variability in output
    if std(log10(PI.data(i).simOutput(:,colIndx)))>1
        set(gca,'YScale','log')
    end
        set(gca,'YScale','log')

end
