%% Laminar - stimulus type comparison


clear

%% Setup
sessionLabel = '211008_B_bmcBRFS001';

% Directories
codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\Brock Carlson\Documents\MATLAB\FormattedDataOutputs';

cd(outputDir)
load(strcat('sortedData_',sessionLabel,'.mat'))

% Variables
cond = 1; % Using most excitatory stimulus, condition 1, 'Simult. Dioptic. PO'
probeLength = size(IDX(cond).LFP_bb{1,1},2);
timeLength = length(STIM(1).sdftm);
trlLength = size(IDX(cond).LFP_bb,1);
xAxisTime = STIM(1).sdftm;
idxPS = 12;
idxNS = 11;

%% Create timetable
% Organize matrix, (15ch,timeseries) for each stimulus
PS = nan(15,1001);
NS = nan(15,1001);

% Create matrix of all trials (ch x time x trial)
psTrlLength = size(IDX(12).correctTrialIndex,1);
for psTrl = 1:psTrlLength
    PS_CHxMSxTRL(:,:,psTrl) = IDX(idxPS).MUAe{psTrl,2}(:,:);
end
%trl avg
PS_avg = mean(PS_CHxMSxTRL,3);

%Now convert PS_avg (a double array) into a table
% The goal is to have 32 variables, each as a column, representing a
% different depth. 
youAreHere

Time = [milliseconds(-200):milliseconds(1):milliseconds(800)];
PS_TT = table2timetable(PS_table,'RowTimes',Time);
NS_TT = table2timetable(NS_table,'RowTimes',Time);

%% stackedplot()
stackedplot(PS_TT,NS_TT)