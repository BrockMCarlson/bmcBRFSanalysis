%% laminarStackedPlot_eachColumn.m
% Brock Carlson
% October 31st, 2023
% 
% The goal of this script is to create a .png for every cortical column
% that experienced bmcBRFS. 

% The 2 x 4 subplot (tiles, i = 1:8) is to show the first 800 ms and then
% second 800 of bmcBRFS conditions. Rows 1 and 2 are identical except for
% their physical stimulus reversals. Either row 1 or row 2 will show
% comparisons of the maximally different stimulus preferences.
% Rows are constructured (1x2)x(1x2) = (BRFS) vs (Physical alternation)
% 2 different colors lines are presented for each laminar electrode. Each
% color represents a different adaptation paradigm. Howver, is should be
% noted that in the case of the BRFS test flash there is no physcial
% difference in the retinal stimulation. The only difference between the
% two clors at that point is their history, and, accordingly, the
% perceptual state of the subject. 


%% Setup
sessionLabel = '211008_B_bmcBRFS001';

% Directories
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(outputDir)

% load data
load(strcat('sortedData_',sessionLabel,'.mat'))

%% Variables
clearvars -except IDX STIM sessionLabel
close all

v1Ch = 1:25;
% v1Ch = 7:31;

% figure
figure
sgtitle(sessionLabel,'Interpreter','none')
%% Calculate V1 ch boundaries
% % % granBtm = 23; % channel corresponding to the bottom of layer 4c
% % % v1Top = granBtm - 9;
% % % v1Btm = granBtm + 5;
v1Ch = 1:25;
% v1Ch = 7:31;
for i = 1:25 % only 25 rows allowed for stackedplot()
    columnNames{i} = strcat('ch_',num2str(v1Ch(i)));
end

% Select conditions
idxps = 9;
idxns = 10;

% First 800 ms
subplot(1,4,1)
epoch = 1;

    % Create timetable
    % Create matrix of all trials (time x ch x trial)
    % IDX(CONDITION).LFP_bb{TRIAL,TRIAL PERIOD}(TIME,CHANNEL)
    psTrlLength = size(IDX(idxps).correctTrialIndex,1);
    for psTrl = 1:psTrlLength
        ps_msXchXtrl(:,:,psTrl) = IDX(idxps).LFP_bb{psTrl,epoch}(:,v1Ch); % MUA output is time x ch
    end
    nsTrlLength = size(IDX(idxns).correctTrialIndex,1);
    for nsTrl = 1:nsTrlLength
        ns_msXchXtrl(:,:,nsTrl) = IDX(idxns).LFP_bb{nsTrl,epoch}(:,v1Ch); % MUA output is time x ch
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

    Time = milliseconds(-200):milliseconds(1):milliseconds(800);
    ps_TT = table2timetable(ps_table,'RowTimes',Time);
    ps_TT2 = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
    ns_TT = table2timetable(ns_table,'RowTimes',Time);
    ns_TT2 = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);
    
    % stackedplot()
    stackedplot(ps_TT2,ns_TT2)
    title('BRFS - monocular adaptation')



% Second 800 ms
subplot(1,4,2)
epoch = 2;


    % Create timetable
    % Create matrix of all trials (time x ch x trial)
    % IDX(CONDITION).LFP_bb{TRIAL,TRIAL PERIOD}(TIME,CHANNEL)
    psTrlLength = size(IDX(idxps).correctTrialIndex,1);
    for psTrl = 1:psTrlLength
        ps_msXchXtrl(:,:,psTrl) = IDX(idxps).LFP_bb{psTrl,epoch}(:,v1Ch); % MUA output is time x ch
    end
    nsTrlLength = size(IDX(idxns).correctTrialIndex,1);
    for nsTrl = 1:nsTrlLength
        ns_msXchXtrl(:,:,nsTrl) = IDX(idxns).LFP_bb{nsTrl,epoch}(:,v1Ch); % MUA output is time x ch
    endns_msXchXtrl(:,:,nsTrl) = IDX(idxns).LFP_bb{nsTrl,2}(:,v1Ch); % MUA output is time x ch
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
    Time = milliseconds(600):milliseconds(1):milliseconds(1600);
    ps_TT = table2timetable(ps_table,'RowTimes',Time);
    ps_TT2 = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
    ns_TT = table2timetable(ns_table,'RowTimes',Time);
    ns_TT2 = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);
    
    
    % stackedplot()
    stackedplot(ps_TT2,ns_TT2)
    title('BRFS - test flash')


% Select conditions
idxps = 20;
idxns = 19;

% First 800 ms
subplot(1,4,3)
epoch = 1;

    % Create timetable
    % Create matrix of all trials (time x ch x trial)
    % IDX(CONDITION).LFP_bb{TRIAL,TRIAL PERIOD}(TIME,CHANNEL)
    psTrlLength = size(IDX(idxps).correctTrialIndex,1);
    for psTrl = 1:psTrlLength
        ps_msXchXtrl(:,:,psTrl) = IDX(idxps).LFP_bb{psTrl,epoch}(:,v1Ch); % MUA output is time x ch
    end
    nsTrlLength = size(IDX(idxns).correctTrialIndex,1);
    for nsTrl = 1:nsTrlLength
        ns_msXchXtrl(:,:,nsTrl) = IDX(idxns).LFP_bb{nsTrl,epoch}(:,v1Ch); % MUA output is time x ch
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
    Time = milliseconds(-200):milliseconds(1):milliseconds(800);
    ps_TT = table2timetable(ps_table,'RowTimes',Time);
    ps_TT2 = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
    ns_TT = table2timetable(ns_table,'RowTimes',Time);
    ns_TT2 = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);
    
    % stackedplot()
    stackedplot(ps_TT2,ns_TT2)
    title('Physical alternation - monocular adaptation')



% Second 800 ms
subplot(1,4,4)
epoch = 2;


    % Create timetable
    % Create matrix of all trials (time x ch x trial)
    % IDX(CONDITION).LFP_bb{TRIAL,TRIAL PERIOD}(TIME,CHANNEL)
    psTrlLength = size(IDX(idxps).correctTrialIndex,1);
    for psTrl = 1:psTrlLength
        ps_msXchXtrl(:,:,psTrl) = IDX(idxps).LFP_bb{psTrl,epoch}(:,v1Ch); % MUA output is time x ch
    end
    nsTrlLength = size(IDX(idxns).correctTrialIndex,1);
    for nsTrl = 1:nsTrlLength
        ns_msXchXtrl(:,:,nsTrl) = IDX(idxns).LFP_bb{nsTrl,epoch}(:,v1Ch); % MUA output is time x ch
    endns_msXchXtrl(:,:,nsTrl) = IDX(idxns).LFP_bb{nsTrl,2}(:,v1Ch); % MUA output is time x ch
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
    Time = milliseconds(600):milliseconds(1):milliseconds(1600);
    ps_TT = table2timetable(ps_table,'RowTimes',Time);
    ps_TT2 = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
    ns_TT = table2timetable(ns_table,'RowTimes',Time);
    ns_TT2 = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);
    
    
    % stackedplot()
    stackedplot(ps_TT2,ns_TT2)
    title('Physical Alternation - test flash')




