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

%% Create timetable
% Organize matrix, (15ch,timeseries) for each stimulus
PS = nan(15,1001);
NS = nan(15,1001);

% Create matrix of all trials (ch x time x trial)
psTrlLength = size(IDX(12).correctTrialIndex,1);
for ch = 1:15
    for psTrl = 1:psTrlLength
        PS_CHxMSxTRL(ch,:,psTrl) = IDX(12).MUAe{psTrl,2}(:,ch);
    end
end
%trl avg

Depth1_NS = [98;97.5;97.9;98.1;97.9];
Depth2_NS = [120;111;119;117;116];
NS_table = table(Depth1,Depth2);

Time = [milliseconds(-200):milliseconds(1):milliseconds(800)];
PS_TT = table2timetable(PS_table,'RowTimes',Time);
NS_TT = table2timetable(NS_table,'RowTimes',Time);

%% stackedplot()
stackedplot(PS_TT,NS_TT)