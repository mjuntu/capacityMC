function sekvenssikestotKOOSTETODELLISETAJATS4 = importNeurofile(workbookFile, sheetName, dataLines)
%IMPORTFILE1 Import data from a spreadsheet
%  SEKVENSSIKESTOTKOOSTETODELLISETAJATS4 = IMPORTFILE1(FILE) reads data
%  from the first worksheet in the Microsoft Excel spreadsheet file
%  named FILE.  Returns the data as a table.
%
%  SEKVENSSIKESTOTKOOSTETODELLISETAJATS4 = IMPORTFILE1(FILE, SHEET)
%  reads from the specified worksheet.
%
%  SEKVENSSIKESTOTKOOSTETODELLISETAJATS4 = IMPORTFILE1(FILE, SHEET,
%  DATALINES) reads from the specified worksheet for the specified row
%  interval(s). Specify DATALINES as a positive scalar integer or a
%  N-by-2 array of positive scalar integers for dis-contiguous row
%  intervals.
%
%  Example:
%  sekvenssikestotKOOSTETODELLISETAJATS4 = importfile1("/Users/mikaeljuntunen/Documents_Mikael/Studies/KTM/Masters thesis/DATA/sekvenssikestot_KOOSTE_TODELLISET_AJAT.xlsx", "Neuro todelliset ajat", [1, 381]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 06-Dec-2023 18:07:04

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [1, 381];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 10);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":J" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["Tunniste", "Lyhenne", "Tutkimus", "Sekvenssi", "Aika", "DRB", "DRG", "CS", "WAVE", "NopeinKiihdytys"];
opts.VariableTypes = ["categorical", "categorical", "categorical", "categorical", "double", "double", "string", "string", "string", "double"];

% Specify variable properties
opts = setvaropts(opts, ["DRG", "CS", "WAVE"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Tunniste", "Lyhenne", "Tutkimus", "Sekvenssi", "DRG", "CS", "WAVE"], "EmptyFieldRule", "auto");

% Import the data
sekvenssikestotKOOSTETODELLISETAJATS4 = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":J" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    sekvenssikestotKOOSTETODELLISETAJATS4 = [sekvenssikestotKOOSTETODELLISETAJATS4; tb]; %#ok<AGROW>
end

end