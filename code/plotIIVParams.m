function plotIIVParams(postSamples, PI,varargin)
inputs = inputParser;
inputs.addParameter('name', {PI.par(:).name})
inputs.addParameter('cellnames', {'MOC1' 'MOC2'});
inputs.addParameter('newFig', true);
inputs.addParameter('dim', false);
inputs.addParameter('n_row', 1);
inputs.addParameter('n_col', 1);
inputs.addParameter('figIndx', 1:length(PI.data));
inputs.addParameter('panel', true);
inputs.parse(varargin{:})
inputs=inputs.Results;

name_eta = inputs.name([PI.H.IndividualParams(:).EtaIndex]);
name_beta = inputs.name([PI.H.CellParams(:).EtaIndex]);
name_beta = inputs.name([PI.H.RespParams(:).EtaIndex]);

indivColors = linspecer(length(PI.data));
if ~isempty(PI.H.IndividualParams(1).EtaIndex)
try
    if inputs.newFig
        figure
    end
    if inputs.dim
        n_col = inputs.n_col;
        n_row = inputs.n_row;
    else
        n_col = ceil(sqrt(length(PI.H.IndividualParams)));
        n_row = ceil(length(PI.H.IndividualParams)/n_col);
    end

    for j=1:length(PI.H.IndividualParams)
        subplot(n_row, n_col,inputs.figIndx(j))
        hold on
        for i=1:length(PI.data)
            histogram((postSamples(:,PI.H.IndividualParams(j).Index(i))),...
                'FaceColor',indivColors(i,:),'FaceAlpha', 0.7, 'Normalization', 'probability')
        end
        legend({PI.data(:).Name},'interpreter', 'none')
        ylabel('prob')
        xlabel('Deviations wrt mean parameter [log scale]')
        title(strjoin({'Inter-individual variation in', name_eta{j}},' '))
    end
catch
    close gcf
end
end
if ~isempty(PI.H.CellParams(1).EtaIndex)
try
    cellColors = linspecer(length(PI.H.CellTypes));
    if inputs.newFig
        figure
    end
    if inputs.dim
        n_col = inputs.n_col;
        n_row = inputs.n_row;
    else
        n_col = ceil(sqrt(length(PI.H.CellParams)));
        n_row = ceil(length(PI.H.CellParams)/n_col);
    end
    if inputs.panel
        figIndx = inputs.figIndx;
    else
        figIndx = 1:length(PI.H.CellParams);
    end
     for j=1:length(PI.H.CellParams)
        subplot(n_row, n_col,figIndx(j))
        hold on
        for i=1:length(PI.H.CellParams(j).Index)
            histogram((postSamples(:,PI.H.CellParams(j).Index(i))),...
                'FaceColor',cellColors(i,:),'FaceAlpha', 0.5, 'Normalization',...
                'probability')    
        end
        cellIndx = 1:3;
        cellIndx = PI.H.CellIndx*cellIndx';
        cellIndx = unique(cellIndx, 'stable');
        legend(PI.H.CellTypes(cellIndx),'interpreter', 'none')
        ylabel('prob')
        xlabel('Deviations wrt mean parameter [log scale]')
        title(strjoin({'Inter-model variation in', name_beta{j}},' '))
       % set(gca, 'XScale', 'log')
     end
catch
        close gcf

end
end
if ~isempty(PI.H.RespParams(1).EtaIndex)
try
    cellColors = linspecer(length(PI.H.RespTypes));
    if inputs.newFig
        figure
    end
    if inputs.dim
        n_col = inputs.n_col;
        n_row = inputs.n_row;
    else
        n_col = ceil(sqrt(length(PI.H.RespParams)));
        n_row = ceil(length(PI.H.RespParams)/n_col);
    end
    if inputs.panel
        figIndx = inputs.figIndx;
    else
        figIndx = 1:length(PI.H.RespParams);
    end
     for j=1:length(PI.H.RespParams)
        subplot(n_row, n_col,figIndx(j))
        hold on
        for i=1:length(PI.H.RespParams(j).Index)
            histogram((postSamples(:,PI.H.RespParams(j).Index(i))),...
                'FaceColor',cellColors(i,:),'FaceAlpha', 0.5, 'Normalization',...
                'probability')    
        end
        cellIndx = 1:size(PI.H.RespIndx,2);
        cellIndx = PI.H.RespIndx*cellIndx';
        cellIndx = unique(cellIndx, 'stable');
        legend(PI.H.RespTypes(cellIndx),'interpreter', 'none')
        ylabel('prob')
        xlabel('Deviations wrt mean parameter [log scale]')
        title(strjoin({'Inter-model variation in', name_beta{j}},' '))
       % set(gca, 'XScale', 'log')
     end
catch
        close gcf

end
end
return
% histogram(exp(postSamples(:,PI.H.IndividualParams(1).EtaIndex)),...
%     'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')

