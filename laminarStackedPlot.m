%% Laminar - stimulus type comparison


clear

%% Setup
sessionLabel = '221202_J_bmcBRFS001';

% Directories
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(outputDir)

%% load data
load(strcat('sortedData_',sessionLabel,'.mat'))

%% Variables
cond = 1; % Using most excitatory stimulus, condition 1, 'Simult. Dioptic. PO'
probeLength = size(IDX(cond).LFP_bb{1,1},2);
timeLength = length(STIM(1).sdftm);
trlLength = size(IDX(cond).LFP_bb,1);
xAxisTime = STIM(1).sdftm;
% idxps = 9;
% idxns = 10;
idxps = 20;
idxns = 19;
granBtm = 27; % channel corresponding to the bottom of layer 4c

% granBtm = 10; % channel corresponding to the bottom of layer 4c

%% Calculate V1 ch boundaries
v1Top = granBtm - 9;
v1Btm = granBtm + 5;
v1Ch = v1Top:v1Btm;
v1Ch = v1Btm:-1:v1Top;


%% Create timetable
% Create matrix of all trials (time x ch x trial)
% IDX(CONDITION).MUAe{TRIAL,TRIAL PERIOD}(TIME,CHANNEL)
psTrlLength = size(IDX(idxps).correctTrialIndex,1);
for psTrl = 1:psTrlLength
    ps_msXchXtrl(:,:,psTrl) = IDX(idxps).MUAe{psTrl,2}(:,v1Ch); % MUA output is time x ch
end
nsTrlLength = size(IDX(idxns).correctTrialIndex,1);
for nsTrl = 1:nsTrlLength
    ns_msXchXtrl(:,:,nsTrl) = IDX(idxns).MUAe{nsTrl,2}(:,v1Ch); % MUA output is time x ch
end
%trl avg
ps_avg = mean(ps_msXchXtrl,3);
ps_table = array2table(ps_avg);
ns_avg = mean(ns_msXchXtrl,3);
ns_table = array2table(ns_avg);

%Now convert ps_avg (a double array) into a table (input to timetable must
%be a table)
% The goal is to have 32 variables, each as a column, representing a
% different depth. 
columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
Time = milliseconds(-200):milliseconds(1):milliseconds(800);
ps_TT = table2timetable(ps_table,'RowTimes',Time);
ps_TT2 = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
ns_TT = table2timetable(ns_table,'RowTimes',Time);
ns_TT2 = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);


%% stackedplot()
stackedplot(ps_TT2,ns_TT2)