%% fig4
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur


datetime

%% Setup
clear
% Directories

codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
% % % codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures');
% % % cd(codeDir)
outDir = 'S:\formattedDataOutputs';
% % % outDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\figures_240404';
% % % % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
% % % % dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
% % % cd(dataDir)
% % % % officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
% % % officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);


% Directories
%Old run - on previous laminar alignment
% % codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
% % cd(codeDir)
% % outDir = 'S:\formattedDataOutputs';
% % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
% % cd(dataDir)
% % officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

%New run - new laminar alignment
% (here is hoping that the results did not change)
% codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures');
cd(codeDir)
% outDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\fig4_MUA';
% dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\sortedData_240229';
% cd(dataDir)
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

%% For loop
% % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
% % dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\sortedData_240229';
% % 
% % cd(dataDir)



%% load DATAOUT
cd('C:\Users\Brock Carlson\Box\Manuscripts\Maier')
load("DATAOUT.mat")

%% Laminar align data
aligned_100_ps = nan(1801,100,size(officLamAssign,1));
aligned_100_ns = nan(1801,100,size(officLamAssign,1));
for i = 1:size(officLamAssign,1)
    if i==5 || i==7
        continue
    end
    

    % Calculate V1 ch boundaries
    granBtm = officLamAssign.stFold4c(i); % channel corresponding to the bottom of layer 4c

    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
    if strcmp(string(officLamAssign.ChToUse(i)),"1:32")
        probeTop = 1;
        probeBtm = 32;
    elseif strcmp(string(officLamAssign.ChToUse(i)),"33:64")
        probeTop = 33;
        probeBtm = 64;
    end
    sinkToTopOfProbe = granBtm-probeTop+1;
    sinkToBtmOfProbe = probeBtm-granBtm;
    v1Top_new = 50-sinkToTopOfProbe+1;
    v1Btm_new = 50+sinkToBtmOfProbe;
    v1Ch_new = v1Top_new:v1Btm_new;
    
    aligned_100_ps(:,v1Ch_new,i) = DATAOUT_ps(:,:,i);
    aligned_100_ns(:,v1Ch_new,i) = DATAOUT_ns(:,:,i);
end
% Now we cut down to just 15 channels
aligned_ps = aligned_100_ps(:,41:55,:);
aligned_ns = aligned_100_ns(:,41:55,:);



%% grand averages and SEM on laminar compartments
% Expecting output of 1801x125 (25 penetrations * 5 ch per compartment)
useIdx = squeeze(~isnan(aligned_ps(1,1,:))); 
reshape_ps_S = reshape(aligned_ps(:,1:5,useIdx),[1801,125]);
reshape_ps_G = reshape(aligned_ps(:,6:10,useIdx),[1801,125]);
reshape_ps_I = reshape(aligned_ps(:,11:15,useIdx),[1801,125]);
reshape_ns_S = reshape(aligned_ns(:,1:5,useIdx),[1801,125]);
reshape_ns_G = reshape(aligned_ns(:,6:10,useIdx),[1801,125]);
reshape_ns_I = reshape(aligned_ns(:,11:15,useIdx),[1801,125]);

grandAvg_ps_S = smoothdata(mean(reshape_ps_S,2,"omitmissing"),1,"gaussian",20);
grandAvg_ps_G = smoothdata(mean(reshape_ps_G,2,"omitmissing"),1,"gaussian",20);
grandAvg_ps_I = smoothdata(mean(reshape_ps_I,2,"omitmissing"),1,"gaussian",20);
grandAvg_ns_S = smoothdata(mean(reshape_ns_S,2,"omitmissing"),1,"gaussian",20);
grandAvg_ns_G = smoothdata(mean(reshape_ns_G,2,"omitmissing"),1,"gaussian",20);
grandAvg_ns_I = smoothdata(mean(reshape_ns_I,2,"omitmissing"),1,"gaussian",20);


% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
SEM_ps_S = smoothdata(std(reshape_ps_S,0,2,"omitmissing")/sqrt(125),1,"gaussian",20);
SEM_ps_G = smoothdata(std(reshape_ps_G,0,2,"omitmissing")/sqrt(125),1,"gaussian",20);
SEM_ps_I = smoothdata(std(reshape_ps_I,0,2,"omitmissing")/sqrt(125),1,"gaussian",20);
SEM_ns_S = smoothdata(std(reshape_ns_S,0,2,"omitmissing")/sqrt(125),1,"gaussian",20);
SEM_ns_G = smoothdata(std(reshape_ns_G,0,2,"omitmissing")/sqrt(125),1,"gaussian",20);
SEM_ns_I = smoothdata(std(reshape_ns_I,0,2,"omitmissing")/sqrt(125),1,"gaussian",20);


% Mean +/- SEM
avgPlusSEM_ps_S = grandAvg_ps_S + SEM_ps_S; %pref stim -- mean plus sem 1 
avgMinusSEM_ps_S = grandAvg_ps_S - SEM_ps_S; %pref stim -- mean minus sem 1 
avgPlusSEM_ps_G = grandAvg_ps_G + SEM_ps_G; %pref stim -- mean plus sem 1 
avgMinusSEM_ps_G = grandAvg_ps_G - SEM_ps_G; %pref stim -- mean minus sem 1 
avgPlusSEM_ps_I = grandAvg_ps_I + SEM_ps_I; %pref stim -- mean plus sem 1 
avgMinusSEM_ps_I = grandAvg_ps_I - SEM_ps_I; %pref stim -- mean minus sem 1 
avgPlusSEM_ns_S = grandAvg_ns_S + SEM_ns_S; %pref stim -- mean plus sem 1 
avgMinusSEM_ns_S = grandAvg_ns_S - SEM_ns_S; %pref stim -- mean minus sem 1 
avgPlusSEM_ns_G = grandAvg_ns_G + SEM_ns_G; %pref stim -- mean plus sem 1 
avgMinusSEM_ns_G = grandAvg_ns_G - SEM_ns_G; %pref stim -- mean minus sem 1 
avgPlusSEM_ns_I = grandAvg_ns_I + SEM_ns_I; %pref stim -- mean plus sem 1 
avgMinusSEM_ns_I = grandAvg_ns_I - SEM_ns_I; %pref stim -- mean minus sem 1 



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
    plot(tm_full,grandAvg_ps_S,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,avgPlusSEM_ps_S,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,avgMinusSEM_ps_S,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,grandAvg_ns_S,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,avgPlusSEM_ns_S,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,avgMinusSEM_ns_S,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Supragranular','% Change'})
    title({strcat('Significant transient modulation?_',string(h_S_trans),'_pVal =',string(p_S_trans)),...
       strcat('Significant sustained modulation?_',string(h_S_susta),'_pVal =',string(p_S_susta))},'Interpreter','none')
nexttile
    plot(tm_full,grandAvg_ps_G,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,avgPlusSEM_ps_G,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,avgMinusSEM_ps_G,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,grandAvg_ns_G,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,avgPlusSEM_ns_G,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,avgMinusSEM_ns_G,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Granular','% Change'})
    title({strcat('Significant transient modulation?_',string(h_G_trans),'_pVal =',string(p_G_trans)),...
       strcat('Significant sustained modulation?_',string(h_G_susta),'_pVal =',string(p_G_susta))},'Interpreter','none')
nexttile
    plot(tm_full,grandAvg_ps_I,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,avgPlusSEM_ps_I,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,avgMinusSEM_ps_I,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,grandAvg_ns_I,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,avgPlusSEM_ns_I,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,avgMinusSEM_ns_I,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Infragranular','% Change'})
    title({strcat('Significant transient modulation?_',string(h_I_trans),'_pVal =',string(p_I_trans)),...
       strcat('Significant sustained modulation?_',string(h_I_susta),'_pVal =',string(p_I_susta))},'Interpreter','none')

titleText = {'Grand average of 150 multi-units (29 penetrations) per laminar compartment'};
title(t,titleText,'Interpreter','none')

%save fig
cd(outDir)
figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.png');
saveas(lamCom,figName_lamCom)
figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.svg');
saveas(lamCom,figName_lamCom)




