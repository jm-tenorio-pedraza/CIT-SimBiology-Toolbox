function plotSimOutput(data,colIndx)

n_sim=size(data,1);
n_col=ceil(sqrt(n_sim));
n_row=ceil(n_sim/n_col);
treatments=unique({data.Group});
treatment_colors=linspecer(length(treatments));
figure
for i=1:n_sim
    subplot(n_row,n_col,i)
    hold on
    sim=plot(data(i).simTime, data(i).simOutput(:,colIndx));
    dat=plot(data(i).dataTime, data(i).dataValue(:,colIndx));
    title(data(i).Name)
    
    col_i=treatment_colors(ismember(treatments,data(i).Group),:);
    sim.Color=col_i;
    dat.LineStyle='none';
    dat.MarkerFaceColor=col_i;
    dat.Marker='d';
    
    % Change axis if there is large variability in output
    if std(log(data(i).simOutput(:,colIndx)))>2
        set(gca,'YScale','log')
    end
end
