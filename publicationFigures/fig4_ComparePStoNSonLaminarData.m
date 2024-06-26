%% fig4
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur


datetime

%% Setup
clear
% Directories

% % % % codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
% % % codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures');
% % % cd(codeDir)
% % % % outDir = 'S:\formattedDataOutputs';
% % % outDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\figures_240404';
% % % % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
% % % dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
% % % cd(dataDir)
% % % % officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
% % % officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);


% Directories

codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
cd(codeDir)
outDir = 'S:\formattedDataOutputs';
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

%% For loop
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)



%% load DATAOUT
cd('C:\Users\Brock Carlson\Box\Manuscripts\Maier')
load("DATAOUT.mat")

%% Laminar align data
aligned_100_ps = nan(1801,100,size(officLamAssign,1));
aligned_100_ns = nan(1801,100,size(officLamAssign,1));
for i = 1:size(officLamAssign,1)

    % Calculate V1 ch boundaries
    granBtm = officLamAssign.Probe11stFold4c(i); % channel corresponding to the bottom of layer 4c

    v1Top_new = 50-granBtm+1;
    v1Btm_new = 50+(32-granBtm);
    v1Ch_new = v1Top_new:v1Btm_new;
    
    aligned_100_ps(:,v1Ch_new,i) = DATAOUT_ps(:,1:32,i);
    aligned_100_ns(:,v1Ch_new,i) = DATAOUT_ns(:,1:32,i);
end
% Now we cut down to just 15 channels
aligned_ps = aligned_100_ps(:,41:55,:);
aligned_ns = aligned_100_ns(:,41:55,:);



%% grand averages oand SEM
ps_grandAvg = mean(aligned_ps,3,"omitmissing");
ns_grandAvg = mean(aligned_ns,3,"omitmissing");

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
penetrationNumber = sum(~isnan(aligned_ps(1,1,:)));
ps_SEM = std(aligned_ps,0,3,"omitmissing")/sqrt(penetrationNumber); 
ns_SEM = std(aligned_ns,0,3,"omitmissing")/sqrt(penetrationNumber); 

% Smooth Data
ps_smooth_mean = smoothdata(ps_grandAvg,1,"gaussian",20);
ns_smooth_mean = smoothdata(ns_grandAvg,1,"gaussian",20);
ps_smooth_sem = smoothdata(ps_SEM,1,"gaussian",20);
ns_smooth_sem = smoothdata(ns_SEM,1,"gaussian",20);

% Mean +/- SEM
ps_mps = ps_smooth_mean + ps_smooth_sem; %pref stim -- mean plus sem 1 
ps_mms = ps_smooth_mean - ps_smooth_sem; %pref stim -- mean minus sem 1 
ns_mps = ns_smooth_mean + ns_smooth_sem; %pref stim -- mean plus sem 1 
ns_mms = ns_smooth_mean - ns_smooth_sem; %pref stim -- mean minus sem 1 

%% Now convert grand average array into laminar compartments
ps_S_mean = mean(ps_smooth_mean(:,1:5),2); 
ps_G_mean = mean(ps_smooth_mean(:,6:10),2); 
ps_I_mean = mean(ps_smooth_mean(:,11:15),2); 
ns_S_mean = mean(ns_smooth_mean(:,1:5),2); 
ns_G_mean = mean(ns_smooth_mean(:,6:10),2); 
ns_I_mean = mean(ns_smooth_mean(:,11:15),2); 

% Mean plus sem
ps_S_mps = mean(ps_mps(:,1:5),2); 
ps_G_mps = mean(ps_mps(:,6:10),2); 
ps_I_mps = mean(ps_mps(:,11:15),2); 
ns_S_mps = mean(ns_mps(:,1:5),2); 
ns_G_mps = mean(ns_mps(:,6:10),2); 
ns_I_mps = mean(ns_mps(:,11:15),2); 

% Mean minus sem
ps_S_mms = mean(ps_mms(:,1:5),2); 
ps_G_mms = mean(ps_mms(:,6:10),2); 
ps_I_mms = mean(ps_mms(:,11:15),2); 
ns_S_mms = mean(ns_mms(:,1:5),2); 
ns_G_mms = mean(ns_mms(:,6:10),2); 
ns_I_mms = mean(ns_mms(:,11:15),2); 

%% statistics
% Ok, the data is together for plotting, now lets run statistics on
% each laminar compartment to see if perceptual modulation occurs. The
% goal here is to run a t-test to see if the average response between
% 1200 and 1600ms significantly differs between ps and ns
useIdx = squeeze(~isnan(aligned_ps(1,1,:))); 
tInput_ps_S_trans = reshape(squeeze(mean(aligned_ps(1050:1200,1:5,useIdx),1)),[],1);
tInput_ps_G_trans = reshape(squeeze(mean(aligned_ps(1050:1200,6:10,useIdx),1)),[],1);
tInput_ps_I_trans = reshape(squeeze(mean(aligned_ps(1050:1200,11:15,useIdx),1)),[],1);
tInput_ns_S_trans = reshape(squeeze(mean(aligned_ns(1050:1200,1:5,useIdx),1)),[],1);
tInput_ns_G_trans = reshape(squeeze(mean(aligned_ns(1050:1200,6:10,useIdx),1)),[],1);
tInput_ns_I_trans = reshape(squeeze(mean(aligned_ns(1050:1200,11:15,useIdx),1)),[],1);
tInput_ps_S_susta = reshape(squeeze(mean(aligned_ps(1400:1801,1:5,useIdx),1)),[],1);
tInput_ps_G_susta = reshape(squeeze(mean(aligned_ps(1400:1801,6:10,useIdx),1)),[],1);
tInput_ps_I_susta = reshape(squeeze(mean(aligned_ps(1400:1801,11:15,useIdx),1)),[],1);
tInput_ns_S_susta = reshape(squeeze(mean(aligned_ns(1400:1801,1:5,useIdx),1)),[],1);
tInput_ns_G_susta = reshape(squeeze(mean(aligned_ns(1400:1801,6:10,useIdx),1)),[],1);
tInput_ns_I_susta = reshape(squeeze(mean(aligned_ns(1400:1801,11:15,useIdx),1)),[],1);

[h_S_trans,p_S_trans] = ttest2(tInput_ps_S_trans,tInput_ns_S_trans);
[h_G_trans,p_G_trans] = ttest2(tInput_ps_G_trans,tInput_ns_G_trans);
[h_I_trans,p_I_trans] = ttest2(tInput_ps_I_trans,tInput_ns_I_trans);
[h_S_susta,p_S_susta] = ttest2(tInput_ps_S_susta,tInput_ns_S_susta);
[h_G_susta,p_G_susta] = ttest2(tInput_ps_G_susta,tInput_ns_G_susta);
[h_I_susta,p_I_susta] = ttest2(tInput_ps_I_susta,tInput_ns_I_susta);


%% Figure generation! 
% tiledLayout plot
close all
tm_full = -200:1600; % 1801 total timepoints
lamCom = figure;
set(gcf,"Position",[1000 123.6667 757.6667 1.1140e+03])
t = tiledlayout(3,1);
nexttile
    plot(tm_full,ps_S_mean,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_S_mps,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_S_mms,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_S_mean,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_S_mps,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_S_mms,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Supragranular','% Change'})
    title({strcat('Significant transient modulation?_',string(h_S_trans),'_pVal =',string(p_S_trans)),...
       strcat('Significant sustained modulation?_',string(h_S_susta),'_pVal =',string(p_S_susta))},'Interpreter','none')
nexttile
    plot(tm_full,ps_G_mean,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_G_mps,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_G_mms,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_G_mean,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_G_mps,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_G_mms,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Granular','% Change'})
    title({strcat('Significant transient modulation?_',string(h_G_trans),'_pVal =',string(p_G_trans)),...
       strcat('Significant sustained modulation?_',string(h_G_susta),'_pVal =',string(p_G_susta))},'Interpreter','none')
nexttile
    plot(tm_full,ps_I_mean,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_I_mps,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_I_mms,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_I_mean,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_I_mps,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_I_mms,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':');
    ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Infragranular','% Change'})
    title({strcat('Significant transient modulation?_',string(h_I_trans),'_pVal =',string(p_I_trans)),...
       strcat('Significant sustained modulation?_',string(h_I_susta),'_pVal =',string(p_I_susta))},'Interpreter','none')

titleText = {'Grand average of 140 multi-units (28 penetrations) per laminar compartment'};
title(t,titleText,'Interpreter','none')

%save fig
cd(outDir)
figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.png');
saveas(lamCom,figName_lamCom)
figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.svg');
saveas(lamCom,figName_lamCom)




