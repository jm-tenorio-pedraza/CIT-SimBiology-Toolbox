function [corrInOut] = plotInputToOutput(inputs,outputs,PI,groups,names,varargin)
par = inputParser;
par.addParameter('plotOutput', false)
par.parse(varargin{:})
par=par.Results;
ncol = ceil(sqrt(size(inputs,2)));
nrow = ceil(size(inputs,2)/ncol);
colors = linspecer(length(groups));
corrInOut = zeros(length(groups),size(inputs,2));
if par.plotOutput
    figure('Name', 'Input-output associations', 'Renderer', 'painters', 'Position', [10 10 1200 800])

for k=1:size(outputs)
    output_k = char(outputs{k});
    for i=1:size(inputs,2)
        subplot(nrow,ncol,i)
        hold on
        for j=1:length(groups)
            h=plot(inputs(:,i), PI.output(j).(output_k)(:,end),'o');
            h(1:end).MarkerFaceColor = repmat(colors(j,:),length(h),1);
          
        end
        set(gca,'XScale', 'log')
        ylim([-100 100])
        title(names{i},'interpreter','none')
    end
end
legend(groups)
end
for k=1:size(outputs)
    output_k = char(outputs{k});
        for j=1:length(groups)
            nanIndx = isnan(PI.output(j).(output_k)(:,end));
            corrInOut(j,:) = partialcorri((inputs(~nanIndx,:)), PI.output(j).(output_k)(~nanIndx,end),...
                'Type','Spearman');
        end
        corrInOut = array2table(corrInOut);
        corrInOut.Properties.VariableNames = names;
        corrInOut.Properties.RowNames = groups;
end
