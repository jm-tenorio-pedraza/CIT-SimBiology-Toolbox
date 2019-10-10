function PK = postProcessPKTable(PK)
Time = round(PK{1:size(PK,1)/2,1});
Mean = PK{1:size(PK,1)/2,2:2:end};
SD = PK{size(PK,1)/2+1:end,2:2:end} - Mean;

Names = PK.Properties.VariableNames(2:2:end);
M = nan(size(Time,1),2*size(Mean,2));
M(:,1:2:end) = Mean;
M(:,2:2:end) = SD;

N=repelem({'a'},1,size(M,2));
N(1:2:length(Names)*2) = Names;
N(2:2:end) = cellfun(@(x)strjoin({'SD', num2str(x)},''),...
    num2cell(1:size(SD,2)),'UniformOutput',false);

M = [Time, M];
PK = array2table([M]);

PK.Properties.VariableNames = ['Time' N];
return

