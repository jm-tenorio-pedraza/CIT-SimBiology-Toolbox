function varNames = getVarNames(PI, stateVar)
if ~isempty(PI.H.IndividualParams(1).Index)
        indivSigmaNames=arrayfun(@(x)strjoin({'omega', x.name}, '_'),PI.H.IndividualParams,'UniformOutput',false)';
else
    indivSigmaNames = [];
end
if ~isempty(PI.H.CellParams(1).Index)
    cellSigmaNames=arrayfun(@(x)strjoin({'psi', x.name}, '_'),PI.H.CellParams,'UniformOutput',false)';
else
    cellSigmaNames = [];
end
try
varNames = [cellSigmaNames; indivSigmaNames];
varNames(end+1:end+length(stateVar)) =  cellfun(@(x) strjoin({'sigma', x}, '_'),...
    stateVar','UniformOutput', false);
catch
    varNames=cellfun(@(x) strjoin({'sigma', x}, '_'),...
    stateVar','UniformOutput', false);
end

if size(varNames,2)>size(varNames,1)
    varNames = varNames';
end
return
