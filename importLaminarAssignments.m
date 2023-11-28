function officialLaminarAssignmentbmcBRFS = importLaminarAssignments(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  OFFICIALLAMINARASSIGNMENTBMCBRFS = IMPORTFILE(FILE) reads data from
%  the first worksheet in the Microsoft Excel spreadsheet file named
%  FILE.  Returns the data as a table.
%
%  OFFICIALLAMINARASSIGNMENTBMCBRFS = IMPORTFILE(FILE, SHEET) reads from
%  the specified worksheet.
%
%  OFFICIALLAMINARASSIGNMENTBMCBRFS = IMPORTFILE(FILE, SHEET, DATALINES)
%  reads from the specified worksheet for the specified row interval(s).
%  Specify DATALINES as a positive scalar integer or a N-by-2 array of
%  positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  officialLaminarAssignmentbmcBRFS = importfile("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "Sheet1", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 27-Nov-2023 15:31:31

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["SessionProbe", "Probe11stFold4c", "Probe12ndFold", "Probe21stFold", "Probe22ndFold"];
opts.VariableTypes = ["string", "double", "categorical", "categorical", "categorical"];

% Specify variable properties
opts = setvaropts(opts, "SessionProbe", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["SessionProbe", "Probe12ndFold", "Probe21stFold", "Probe22ndFold"], "EmptyFieldRule", "auto");

% Import the data
officialLaminarAssignmentbmcBRFS = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    officialLaminarAssignmentbmcBRFS = [officialLaminarAssignmentbmcBRFS; tb]; %#ok<AGROW>
end

end