function PI = initialValues(PI,TV_0,Theta,varargin)
inputs=inputParser;
inputs.addParameter('NR','remove')
inputs.parse(varargin{:})
inputs=inputs.Results;
%% Identify point when TV exceeds TV_0
% Convert matrix to cell
TV=mat2cell(PI.output.Tumor,ones(size(PI.output.Tumor,1),1));
TV_0=num2cell(TV_0);
% Create structure to compare initial TV to the time course of TV for each
% sample
TV_struct=[];
[TV_struct(1:size(PI.output.Tumor,1)).TV]=TV{:,:};
[TV_struct(1:size(PI.output.Tumor,1)).TV_0]=TV_0{:,:};
TV_indx=arrayfun(@(x)find(x.TV>=x.TV_0,1),TV_struct,'UniformOutput',false)';
% What to do with samples where TV_0 is not observed
% Samples where TV_0 is not observed use the last time point
if strcmp(inputs.NR,'remove')
    TV_indx_empty=cellfun(@(x) isempty(x),TV_indx, 'UniformOutput',true);
    for i=1:length(PI.observablesPlot)
        output_i=PI.observablesPlot{i};
        mat_i = PI.output(1).(output_i)(~TV_indx_empty,:);
        PI.output(1).(output_i)=mat_i;
    end
    PI.theta=Theta(~TV_indx_empty,:);
    TV_indx=TV_indx(~TV_indx_empty);
else
    TV_indx_empty=cellfun(@(x) isempty(x),TV_indx);
    TV_indx(TV_indx_empty) = {length(PI.tspan)};
    % Create dummy matrix
end
    TV_indx=cellfun(@(x)1:length(PI.tspan)==x,TV_indx, 'UniformOutput',false);

    TV_indx=cell2mat(TV_indx);

%% Select initial values
initialValues = nan(size(PI.theta,1),length(PI.observablesPlot));
for i=1:length(PI.observablesPlot)
        output_i=PI.observablesPlot{i};
        initial_i = sum(PI.output(1).(output_i).*TV_indx,2);
        initialValues(:,i) = initial_i;
end

%% Identify samples with NaNs in initial values
naIndx=or(or(or(isnan(initialValues), ~isreal(initialValues)),...
    isinf(initialValues)),initialValues<0);
nanIndx=sum(naIndx,2)>0;

% Remove samples with NaN
theta_red = PI.theta(~nanIndx,:);
for i=1:length(PI.observablesPlot)
        output_i=PI.observablesPlot{i};
        PI.output.(output_i)=PI.output.(output_i)(~nanIndx,:);
end
PI.initialValues = initialValues(~nanIndx,:);
PI.theta=theta_red;
end

