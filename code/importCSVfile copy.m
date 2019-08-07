function dataTable = importCSVfile(filename, VariableNames,VariableTypes, dataLines)
%IMPORTFILE1 Import data from a .csv file

%% Input handling
if nargin<3
    error('importCSVfile requires a VariableNames and a VariableTypes cell object')
end
% If dataLines is not specified, define defaults
if nargin < 4
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", length(VariableNames));

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = VariableNames;
opts.VariableTypes = VariableTypes;
% opts = setvaropts(opts, 4, "TrimNonNumeric", true);
% opts = setvaropts(opts, 4, "ThousandsSeparator", ",");
opts = setvaropts(opts, 4, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
dataTable = readtable(filename, opts);

end