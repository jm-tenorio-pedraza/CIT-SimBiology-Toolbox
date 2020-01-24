function plotIIVParams(postSamples, PI,varargin)
inputs = inputParser;
inputs.addParameter('name', {PI.par(:).name})
inputs.addParameter('cellnames', {'MOC1' 'MOC2'});
inputs.parse(varargin{:})
inputs=inputs.Results;

name_eta = inputs.name([PI.H.IndividualParams(:).EtaIndex]);
name_beta = inputs.name([PI.H.CellParams(:).EtaIndex]);

try
    for j=1:length(PI.H.IndividualParams)
        figure
        hold on
        
        for i=1:length(PI.data)
            histogram((postSamples(:,PI.H.IndividualParams(j).Index(i))),...
                'FaceColor',PI.data(i).colors,'FaceAlpha', 0.7, 'Normalization', 'probability')
            
        end
        legend({PI.data(:).Name},'interpreter', 'none')
ylabel('prob')
xlabel('Deviations wrt mean parameter [log scale]')
title(strjoin({'Inter-individual variation in', name_eta{j}},' '))
    end
catch
    close gcf
    

    
end
try
     for j=1:length(PI.H.CellParams)
        figure
        hold on
        
        for i=1:length(PI.H.CellParams(j).Index)
            histogram((postSamples(:,PI.H.CellParams(j).Index(i))),...
                'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization',...
                'probability')
            
        end
        legend(inputs.cellnames,'interpreter', 'none')
        ylabel('prob')
        xlabel('Deviations wrt mean parameter [log scale]')
        title(strjoin({'Inter-model variation in', name_beta{j}},' '))
       % set(gca, 'XScale', 'log')
     end
catch
        close gcf

end
return
% histogram(exp(postSamples(:,PI.H.IndividualParams(1).EtaIndex)),...
%     'FaceColor',PI.data(i).colors,'FaceAlpha', 0.5, 'Normalization', 'probability')

