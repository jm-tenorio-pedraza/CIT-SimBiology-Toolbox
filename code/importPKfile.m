function PK = importPKfile(filename,VariableNames,VariableTypes,dataLines)


%% Input handling

% If dataLines is not specified, define defaults
if nargin < 4
    dataLines = [1, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", length(VariableNames));

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableTypes = VariableTypes;
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
PK = readtable(filename, opts);
PK.Properties.VariableNames = VariableNames;
end