%% Subset HumanPI
N = 133;
HumanPI_Red = HumanPI;
for i=2:length(HumanPI)
    sampleIndx = randsample(size(HumanPI(i).Theta,1),N);
    HumanPI_Red(i).Theta = HumanPI(i).Theta(sampleIndx,:);
%     name = {HumanPI(i).PI.output(1:end).Name};
%     [HumanPI_Red(i).PI.output(1:end).Name] =  name{:,:};
    for j=1:length(PI.observablesPlot)
        output_i = PI.observablesPlot{j};
        simoutput = arrayfun(@(x)x.(output_i)(sampleIndx,:),HumanPI(i).PI.output,'UniformOutput',false);
        [HumanPI_Red(i).PI.output(1:end).(output_i)]= simoutput{:,:};
    end
end
%% Sigma
sigma_indx= ismember({Meta(1).Struct.PI.par(:).name},'sigma_Tumor');
sigma = arrayfun(@(x)exp(randsample(Meta(1).Struct.PI.postSamples(:,sigma_indx), ...
    size(x.Theta,1),true,ones(size(Meta(1).Struct.PI.postSamples,1),1))),HumanPI_Red,'UniformOutput',false);
[HumanPI_Red(1:end).Sigma]=sigma{:,:};

%%
for i=1:length(HumanPI_Red)
   plotSimulations(HumanPI_Red(i).PI,'YScale','log')
end
%% Obtain SLD change and PFS
treatments = {'Control' 'antiPDL1' 'antiCTLA4' 'antiPDL1_antiCTLA4'};
cutoff=24*30;
% Calculate SLD and response %
for i =1:length(HumanPI_Red)
    HumanPI_Red(i).PI=getSLD(HumanPI_Red(i).PI, HumanPI_Red(i).Sigma, 'TV_0',...
        HumanPI_Red(i).PI.output(1).Tumor(1:end,1));
    [HumanPI_Red(i).PI, HumanPI_Red(i).ORR] = getORR(HumanPI_Red(i).PI, treatments,'TV_0',...
        HumanPI_Red(i).PI.output(1).Tumor(1:end,1),'N',size( HumanPI_Red(i).Theta,1),'cutoff_value',24*30);
    [HumanPI_Red(i).PI, HumanPI_Red(i).Response] = getPFS(HumanPI_Red(i).PI, treatments,'cutoff_value',cutoff);
    HumanPI_Red(i).PI = getSurvivalTime(HumanPI_Red(i).PI,...
        treatments,'N',size(HumanPI_Red(i).Theta,1),'cutoff_value',cutoff);
    HumanPI_Red(i).PI = getSurvivalTime(HumanPI_Red(i).PI,...
        treatments,'N',size(HumanPI_Red(i).Theta,1),'tumorField','Tumor','survivalType','OS','CC',...
        75);
end
%% Plot SLD and PFS
for i = 1:length(HumanPI_Red)
    plotORR(HumanPI_Red(i).PI, treatments,'output', {'SLD'})
end

for i = 1:length(HumanPI_Red)
%     figure
    plotSurvivalFunction(HumanPI_Red(i).PI,30*24,treatments)
end
%% Comparing ORR
ORR = {HumanPI_Red(1:end).Response};
parametrizations=num2cell(1:length(HumanPI_Red));
parametrizationNames =   cellfun(@(x)strjoin({'Par' num2str(x)},''),parametrizations,'UniformOutput',false);
{'Par1' 'Par2' 'Par3' 'Par4' 'Par5' 'Par6' 'Par7'...
    'Par8' 'Par9' 'Par10' 'Par11' 'Par12' 'Par13'};
controlORR = array2table(cell2mat(cellfun(@(x) x{:,2}, ORR,...
    'UniformOutput', false)), 'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames);
antiPDL1ORR = array2table(cell2mat(cellfun(@(x) x{:,3}, ORR, ...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames );
antiCTLA4ORR = array2table(cell2mat(cellfun(@(x) x{:,4}, ORR,...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames) ;
antiPDL1_antiCTLA4ORR = array2table(cell2mat(cellfun(@(x) x{:,5}, ORR, ...
    'UniformOutput', false)),'RowNames', {'PD' 'SD' 'PR' 'CR'}, 'VariableNames',...
  parametrizationNames);
antiPDL1ORR{:,:}-controlORR{:,:}
antiCTLA4ORR{:,:}-controlORR{:,:}
antiPDL1_antiCTLA4ORR{:,:}-controlORR{:,:}
%% 
pfs = arrayfun(@(x)getMedianPFS(x.PI,treatments),HumanPI,'UniformOutput',false);
[HumanPI(1:end).PFS]=pfs{:,:};