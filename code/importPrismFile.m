function dataMatrix = importPrismFile(filename, dataLines, varNames,startVar)
%% Input handling

%% Setup the Import Options

% Specify column names and types
if startVar==0
    n_vars=length(varNames);
    opts = delimitedTextImportOptions("NumVariables", n_vars);

    opts.VariableNames = varNames;
    opts.VariableTypes = repelem({'double'},n_vars);
else
    n_vars=length(varNames)+1;
    opts = delimitedTextImportOptions("NumVariables", n_vars);

    opts.VariableNames = ["Baseline" varNames];
    opts = setvaropts(opts, startVar, "WhitespaceRule", "preserve");
    opts.VariableTypes = ['string' repelem({'double'},n_vars-1)];
end
% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

opts.SelectedVariableNames = varNames;

opts = setvaropts(opts, 2:n_vars, "TrimNonNumeric", true);
opts = setvaropts(opts, 2:n_vars, "DecimalSeparator", ",");
opts = setvaropts(opts, 2:n_vars, "ThousandsSeparator", ".");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
dataMatrix = readtable(filename, opts);

end