%% fig4
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur


datetime

%% Setup
clear
% Directories

% codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures');
cd(codeDir)
% outDir = 'S:\formattedDataOutputs';
outDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\figures_240402';
% dataDir = 'S:\bmcBRFS_sortedData_Nov23';
dataDir = 'D:\sortedData_240229';
cd(dataDir)
            size(array_ofMonoc4,3)];
% officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);


%% load DATAOUT
cd(outDir)
load("DATAOUT.mat")

%% Align to sink bottom
This is where you need to work

%% plot grand averages of whole probe


    % convert 0 to NaN
    ps_NaNmatrix = DATAOUT.MUAe.ps_PercentC;
    ps_NaNmatrix(ps_NaNmatrix==0) = NaN;
    ns_NaNmatrix = DATAOUT.MUAe.ns_PercentC;
    ns_NaNmatrix(ns_NaNmatrix==0) = NaN;

    ps_grandAvg = median(ps_NaNmatrix,3,"omitmissing");
    ns_grandAvg = median(ns_NaNmatrix,3,"omitmissing");

    % Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
    ps_SEM = std(ps_NaNmatrix,0,3)/sqrt(size(ps_NaNmatrix,3)); 
    ns_SEM = std(ns_NaNmatrix,0,3)/sqrt(size(ns_NaNmatrix,3)); 

    % smooth data
    ps_smooth_grandAvg = smoothdata(ps_grandAvg,1,"gaussian",20);
    ns_smooth_grandAvg = smoothdata(ns_grandAvg,1,"gaussian",20);

    % convert to table
    ps_table_grandAvg = array2table(ps_smooth_grandAvg);
    ns_table_grandAvg = array2table(ns_smooth_grandAvg);

    % columnNames
    columnNames = {'1','2','3','4','5','6','7','8','9','10','11','12',...
        '13','14','15','16','17','18','19','20','21','22','23','24',...
        '25','26','27','28','29','30','31','32'};

    %Now convert ps_avg (a double array) into a table (input to timetable must
    %be a table)
    % The goal is to have 32 variables, each as a column, representing a
    % different depth. 
    Time = milliseconds(-200):milliseconds(1):milliseconds(1600);
    ps_TT_grandAvg = table2timetable(ps_table_grandAvg,'RowTimes',Time);
    preferredStimFlash_grandAvg = renamevars(ps_TT_grandAvg,ps_TT_grandAvg.Properties.VariableNames,columnNames);
    ns_TT_grandAvg = table2timetable(ns_table_grandAvg,'RowTimes',Time);
    nullStimFlash_grandAvg = renamevars(ns_TT_grandAvg,ns_TT_grandAvg.Properties.VariableNames,columnNames);


    % stackedplot()
    close all
    stk_grandAvg = figure;
    set(stk_grandAvg,"Position",[1000 60.3333 560 1.2933e+03])
    s_grandAvg = stackedplot(preferredStimFlash_grandAvg,nullStimFlash_grandAvg);
    titleText_grandAvg = {'_grandAvg'};
    s_grandAvg.Title = titleText_grandAvg;
    s_grandAvg.LineWidth = 1;
    % s_grandAvg.DisplayLabels = ["% Change"];

    %save fig
    cd(outDir)
    figName_grandAvg = strcat('stackedPlot_','_grandAvg_','.png');
    saveas(stk_grandAvg,figName_grandAvg)


%% plot grand averages of laminar compartments



    % convert 0 to NaN
    ps_NaNmatrix = DATAOUT.MUAe.ps_PercentC;
    ps_NaNmatrix(ps_NaNmatrix==0) = NaN;
    ns_NaNmatrix = DATAOUT.MUAe.ns_PercentC;
    ns_NaNmatrix(ns_NaNmatrix==0) = NaN;

    ps_grandAvg = median(ps_NaNmatrix,3,"omitmissing");
    ns_grandAvg = median(ns_NaNmatrix,3,"omitmissing");

    % Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
    ps_SEM = std(ps_NaNmatrix,0,3)/sqrt(size(ps_NaNmatrix,3)); 
    ns_SEM = std(ns_NaNmatrix,0,3)/sqrt(size(ns_NaNmatrix,3)); 

    %Now convert grand average array into laminar compartments
    ps_S = median(ps_grandAvg(:,1:5),2); 
    ps_G = median(ps_grandAvg(:,6:10),2); 
    ps_I = median(ps_grandAvg(:,11:15),2); 
    ns_S = median(ns_grandAvg(:,1:5),2); 
    ns_G = median(ns_grandAvg(:,6:10),2); 
    ns_I = median(ns_grandAvg(:,11:15),2); 


    % smooth data
    ps_S_smooth = smoothdata(ps_S,1,"gaussian",20);
    ps_G_smooth = smoothdata(ps_G,1,"gaussian",20);
    ps_I_smooth = smoothdata(ps_I,1,"gaussian",20);
    ns_S_smooth = smoothdata(ns_S,1,"gaussian",20);
    ns_G_smooth = smoothdata(ns_G,1,"gaussian",20);
    ns_I_smooth = smoothdata(ns_I,1,"gaussian",20);

    % Ok, the data is together for plotting, now lets run statistics on
    % each laminar compartment to see if perceptual modulation occurs. The
    % goal here is to run a t-test to see if the average response between
    % 1200 and 1600ms significantly differs between ps and ns
    useIdx = squeeze(~isnan(ps_NaNmatrix(1,1,:))); 
    tInput_ps_S = reshape(squeeze(median(ps_NaNmatrix(1400:1801,1:5,useIdx),1)),[],1);
    tInput_ps_G = reshape(squeeze(median(ps_NaNmatrix(1400:1801,6:10,useIdx),1)),[],1);
    tInput_ps_I = reshape(squeeze(median(ps_NaNmatrix(1400:1801,11:15,useIdx),1)),[],1);
    tInput_ns_S = reshape(squeeze(median(ns_NaNmatrix(1400:1801,1:5,useIdx),1)),[],1);
    tInput_ns_G = reshape(squeeze(median(ns_NaNmatrix(1400:1801,6:10,useIdx),1)),[],1);
    tInput_ns_I = reshape(squeeze(median(ns_NaNmatrix(1400:1801,11:15,useIdx),1)),[],1);

    [h_S,p_S] = ttest2(tInput_ps_S,tInput_ns_S);
    [h_G,p_G] = ttest2(tInput_ps_G,tInput_ns_G);
    [h_I,p_I] = ttest2(tInput_ps_I,tInput_ns_I);

    % tiledLayout plot
    close all
    tm_full = -200:1600; % 1801 total timepoints
    lamCom = figure;
    set(gcf,"Position",[1000 503 560 734.6667])
    t = tiledlayout(3,1);
    nexttile
        plot(tm_full,ps_S_smooth,tm_full,ns_S_smooth,'LineWidth',1.5)
        vline(0)
        vline(800)
        xregion(1200,1600)
        ylabel({'Supragranular','% Change'})
        title({strcat('Significant perceptual modulation?_',string(h_S)),strcat('pVal =',string(p_S))},'Interpreter','none')
    nexttile
        plot(tm_full,ps_G_smooth,tm_full,ns_G_smooth,'LineWidth',1.5)
        vline(0)
        vline(800)
        xregion(1200,1600)
        ylabel({'Granular','% Change'})
        title({strcat('Significant perceptual modulation?_',string(h_G)),strcat('pVal =',string(p_G))},'Interpreter','none')
    nexttile
        plot(tm_full,ps_I_smooth,tm_full,ns_I_smooth,'LineWidth',1.5)
        vline(0)
        vline(800)
        xregion(1200,1600)
        ylabel({'Supragranular','% Change'})
        title({strcat('Significant perceptual modulation?_',string(h_I)),strcat('pVal =',string(p_I))},'Interpreter','none')

    titleText = {'Grand average of 135 multi-units per laminar compartment'};
    title(t,titleText,'Interpreter','none')

    %save fig
    cd(outDir)
    figName_lamCom = strcat('laminarCompartment_','_grandAvg_','.png');
    saveas(lamCom,figName_lamCom)



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

%% Laminar align data
cd('C:\Users\Brock Carlson\Box\Manuscripts\Maier')
load('DATAOUT.mat')
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
aligned_ps = nan(1801,64,size(officLamAssign,1));
aligned_ns = nan(1801,64,size(officLamAssign,1));
for i = 1:size(officLamAssign,1)
    % Calculate V1 ch boundaries
    granBtm = officLamAssign.Probe11stFold4c(i); % channel corresponding to the bottom of layer 4c
    v1Top_old = granBtm-9;
    v1Btm_old = granBtm+5;
    v1Ch_old = v1Top_old:v1Btm_old;
    
    numPenetrations = 30;
    alignDiff = 32-granBtm;
    v1Top_new = v1Top_old+alignDiff;
    v1Btm_new = v1Btm_old+alignDiff;
    v1Ch_new = v1Top_new:v1Btm_new;
    
    aligned_ps(:,v1Ch_new,i) = DATAOUT_ps(:,v1Ch_old,i);
    aligned_ns(:,v1Ch_new,i) = DATAOUT_ns(:,v1Ch_old,i);
end



%% Step 4, Show variance with SEM
% Average data
ps_avg = mean(DATAOUT_ps,3,"omitmissing"); % Avrage acros penetrations (other steps were done with median
ns_avg = mean(DATAOUT_ns,3,"omitmissing");

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
% % ps_SEM = std(DATAOUT_ps,0,3,"omitmissing")/sqrt(size(DATAOUT_ps,3)); 
% % ns_SEM = std(DATAOUT_ns,0,3,"omitmissing")/sqrt(size(DATAOUT_ns,3)); 
ps_SEM = std(DATAOUT_ps,0,3,"omitmissing")/sqrt(135); 
ns_SEM = std(DATAOUT_ns,0,3,"omitmissing")/sqrt(135); 


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

