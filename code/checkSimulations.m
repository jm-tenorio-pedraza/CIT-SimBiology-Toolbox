function PI = checkSimulations(PI,Theta,varargin)
inputs=inputParser;
inputs.addParameter('TV_0', [])
inputs.addParameter('cutoff_value', PI.tspan(end))
inputs.addParameter('TumorField', 'Tumor');
inputs.parse(varargin{:})
inputs=inputs.Results;
indxMat = nan(size(PI.output(1).Tumor,1),length(PI.observablesPlot));
indx = repelem([false], size(indxMat,1),1);
for k=1:length(PI.output)
for i=1:length(PI.observablesPlot)
    output_i = PI.observablesPlot{i};
    nanIndx = all(isnan(PI.output(k).(output_i)),2);
    imIndx = all(~isreal(PI.output(k).(output_i)),2);
    infIndx=all(isinf(PI.output(k).(output_i)),2);
    indx = or(or(nanIndx,or(imIndx,infIndx)),indx);
%     indx=or(sum(or(or(isinf(PI.output(k).(output_i)), ~isreal(PI.output(k).(output_i))),isnan(PI.output(k).(output_i))),2)>0,indx);
end
end
for k=1:length(PI.output)
for i=1:length(PI.observablesPlot)
    output_i = PI.observablesPlot{i};
    PI.output(k).(output_i)=PI.output(k).(output_i)(~indx,:);
end
end
PI.Theta=Theta(~indx,:);
