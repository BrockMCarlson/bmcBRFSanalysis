% notes to move

%% Notes
%
%
%

% What is currently working: 



%% Step 1 = load IDX
% Step 2 = chose NS and PS for BRFS
% Step 3, trial average
ps = 12; % 9;
ns = 11 ; % 10;
j = ps;
for trlNum = 1:length(IDX(j).correctTrialIndex)
    MUAe_ps(1:32,1:1001,trlNum) =  IDX(j).MUAe{trlNum,2}(1:1001,1:32)';
end
ps_avg = median(MUAe_ps,3);

k = ns;
for trlNum = 1:length(IDX(k).correctTrialIndex)
    MUAe_ns(1:32,1:1001,trlNum) =  IDX(k).MUAe{trlNum,2}(1:1001,1:32)';
end
ns_avg = median(MUAe_ns,3);



%% Step 4, Show variance with SEM
% Average data
ps_avg = mean(DATAOUT_ps,3,"omitmissing"); % Avrage acros penetrations (other steps were done with median
ns_avg = mean(DATAOUT_ns,3,"omitmissing");

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
% % ps_SEM = std(DATAOUT_ps,0,3,"omitmissing")/sqrt(size(DATAOUT_ps,3)); 
% % ns_SEM = std(DATAOUT_ns,0,3,"omitmissing")/sqrt(size(DATAOUT_ns,3)); 
ps_SEM = std(DATAOUT_ps,0,3,"omitmissing")/sqrt(140); 
ns_SEM = std(DATAOUT_ns,0,3,"omitmissing")/sqrt(140); 


% smooth data
ps_smooth_mean = smoothdata(ps_avg,1,"gaussian",20);
ns_smooth_mean = smoothdata(ns_avg,1,"gaussian",20);
ps_smooth_sem = smoothdata(ps_SEM,1,"gaussian",20);
ns_smooth_sem = smoothdata(ns_SEM,1,"gaussian",20);

% convert to table
ps_table_1_mean = array2table(ps_smooth_mean(:,1:16));
ns_table_1_mean = array2table(ns_smooth_mean(:,1:16));
ps_table_2_mean = array2table(ps_smooth_mean(:,17:32));
ns_table_2_mean = array2table(ns_smooth_mean(:,17:32));

ps_table_1_sem = array2table(ps_smooth_sem(:,1:16));
ns_table_1_sem = array2table(ns_smooth_sem(:,1:16));
ps_table_2_sem = array2table(ps_smooth_sem(:,17:32));
ns_table_2_sem = array2table(ns_smooth_sem(:,17:32));

%Now convert ps_avg (a double array) into a table (input to timetable must
%be a table)
% The goal is to have 32 variables, each as a column, representing a
% different depth. 
Time = milliseconds(-200):milliseconds(1):milliseconds(1600);
% electrodes 1-16
elName_1 = {'el1','el2','el3','el4','el5','el6','el7','el8','el9',...
    'el10','el11','el12','el13','el14','el15','el16'};
ps_TT_1_mean = table2timetable(ps_table_1_mean,'RowTimes',Time);
preferredStimFlash_1_mean = renamevars(ps_TT_1_mean,ps_TT_1_mean.Properties.VariableNames,elName_1);
ns_TT_1_mean = table2timetable(ns_table_1_mean,'RowTimes',Time);
nullStimFlash_1_mean = renamevars(ns_TT_1_mean,ns_TT_1_mean.Properties.VariableNames,elName_1);
% electrodes 17-32
elName_2 = {'el17','el18','el19','el20','el21','el22','el23','el24',...
    'el25','el26','el27','el28','el29','el30','el31','el32'};
ps_TT_2_mean = table2timetable(ps_table_2_mean,'RowTimes',Time);
preferredStimFlash_2_mean = renamevars(ps_TT_2_mean,ps_TT_2_mean.Properties.VariableNames,elName_2);
ns_TT_2_mean = table2timetable(ns_table_2_mean,'RowTimes',Time);
nullStimFlash_2_mean = renamevars(ns_TT_2_mean,ns_TT_2_mean.Properties.VariableNames,elName_2);    

% SEM
ps_TT_1_sem = table2timetable(ps_table_1_sem,'RowTimes',Time);
preferredStimFlash_1_sem = renamevars(ps_TT_1_sem,ps_TT_1_sem.Properties.VariableNames,elName_1);
ns_TT_1_sem = table2timetable(ns_table_1_sem,'RowTimes',Time);
nullStimFlash_1_sem = renamevars(ns_TT_1_sem,ns_TT_1_sem.Properties.VariableNames,elName_1);
ps_TT_2_sem = table2timetable(ps_table_2_sem,'RowTimes',Time);
preferredStimFlash_2_sem = renamevars(ps_TT_2_sem,ps_TT_2_sem.Properties.VariableNames,elName_2);
ns_TT_2_sem = table2timetable(ns_table_2_sem,'RowTimes',Time);
nullStimFlash_2_sem = renamevars(ns_TT_2_sem,ns_TT_2_sem.Properties.VariableNames,elName_2);    

% Mean +/- SEM
% preferredStimFlash_1_sem
% nullStimFlash_1_sem
ps_mps_1 = preferredStimFlash_1_mean + preferredStimFlash_1_sem; %pref stim -- mean plus sem 1 
ps_mms_1 = preferredStimFlash_1_mean - preferredStimFlash_1_sem; %pref stim -- mean minus sem 1 
ns_mps_1 = nullStimFlash_1_mean + nullStimFlash_1_sem; %pref stim -- mean plus sem 1 
ns_mms_1 = nullStimFlash_1_mean - nullStimFlash_1_sem; %pref stim -- mean minus sem 1 

ps_mps_2 = preferredStimFlash_2_mean + preferredStimFlash_2_sem; %pref stim -- mean plus sem 1 
ps_mms_2 = preferredStimFlash_2_mean - preferredStimFlash_2_sem; %pref stim -- mean minus sem 1 
ns_mps_2 = nullStimFlash_2_mean + nullStimFlash_2_sem; %pref stim -- mean plus sem 1 
ns_mms_2 = nullStimFlash_2_mean - nullStimFlash_2_sem; %pref stim -- mean minus sem 1 


% Step 5, stackedplot()
close all
f = figure;
tl = tiledlayout(3,2);

nexttile(tl,[3 1])
s_1 = stackedplot(preferredStimFlash_1_mean,nullStimFlash_1_mean,...
    ps_mps_1,ps_mms_1,ns_mps_1,ns_mms_1);
s_1.LineWidth = 1;

nexttile(tl,[3 1])
s_2 = stackedplot(preferredStimFlash_2_mean,nullStimFlash_2_mean,...
    ps_mps_2,ps_mms_2,ns_mps_2,ns_mms_2);
s_2.LineWidth = 1;
