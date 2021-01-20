function PK = postProcessPKTable(PK, varargin)
par = inputParser;
par.addParameter('roundDigits', 0)
par.addParameter('VarNames', {})
par.parse(varargin{:})
par=par.Results;
Time = round(PK{1:size(PK,1)/2,1}, par.roundDigits);
Mean = PK{1:size(PK,1)/2,2:2:end};
SD = PK{size(PK,1)/2+1:end,2:2:end} - Mean;
if isempty(par.VarNames)
    Names = PK.Properties.VariableNames(2:2:end);
else 
    Names = par.VarNames;
end
M = nan(size(Time,1),2*size(Mean,2));
M(:,1:2:end) = Mean;
M(:,2:2:end) = SD;

N=repelem({'a'},1,size(M,2));
N(1:2:length(Names)*2) = Names;
N(2:2:end) = cellfun(@(x)strjoin({'SD', num2str(x)},''),...
    num2cell(1:size(SD,2)),'UniformOutput',false);

M = [Time, M];
PK = array2table([M]);
if isempty(par.VarNames)
    PK.Properties.VariableNames = ['Time' N];
else
    PK.Properties.VariableNames = Names;
end
return

