function dataTable = importExcelFile(workbookFile, sheetName, dataLines)
%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, 12];
end

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":L" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["Time", "CD8", "CD107a", "Treg", "DC", "MDSC", "PDL1_Tumor_Rel", "PDL1_Immune_Rel", "DC_Rel", "MDSC_Rel", "Dose", "Group"];
opts.SelectedVariableNames = ["Time", "CD8", "CD107a", "Treg", "DC", "MDSC", "PDL1_Tumor_Rel", "PDL1_Immune_Rel", "DC_Rel", "MDSC_Rel", "Dose", "Group"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical"];
opts = setvaropts(opts, 12, "EmptyFieldRule", "auto");

% Import the data
dataTable = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":L" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    dataTable = [dataTable; tb]; %#ok<AGROW>
end

end