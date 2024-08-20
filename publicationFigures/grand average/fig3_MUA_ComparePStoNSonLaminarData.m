%% fig3
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur

%% Setup
disp('start time')
datetime
close
clearvars -except MUA_trials
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures\grand average';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir\fig3_MUA';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures\grand average';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\fig3_MUA';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)
cd(dataDir)

if ~exist('MUA_trials','var')
    tic
    load('MUA_trials.mat') % format is MUA_trials{penetration,1}{cond,1}{trial,flash}
    toc
end



%%
averageMUAMatrix_BRFSps = nan(2001,15,31);
averageLFPMatrix_BRFSns = nan(2001,15,31);
for penetration = 1:size(MUA_trials,1)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');

    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',fileToLoad))
        continue
    end
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    if strcmp(string(officLamAssign.ChToUse(penetration)),"1:32")
        if any(v1Ch > 32) || any(v1Ch < 1)
            warning('skipping session without full column for now')
            disp(fileToLoad)
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('skipping session without full column for now')
            disp(fileToLoad)
            continue
        end
    end
    
    %% Obtain Monocular preference
    % bl Sub at average level (better for plotting)
    clear array_ofMonoc1 array_ofMonoc2 array_ofMonoc3 array_ofMonoc4
    
    % Monocular
    monoc_1 = [5, 11, 13, 19]; % PO RightEye
    monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    
    % convert from cell to double and combine monocular conditions
    count = 0;
    for cond = monoc_1
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    % Test for monocular preference with anova
    % The anova function ignores NaN values, <undefined> values, empty 
    % characters, and empty strings in y. If factors or tbl contains NaN or 
    % <undefined> values, or empty characters or strings, the function ignores 
    % the corresponding observations in y. The ANOVA is balanced if each factor 
    % value has the same number of observations after the function disregards 
    % empty or NaN values. Otherwise, the function performs an unbalanced 
    % ANOVA. (We will be doing an unbalanced ANOVA).
    
    % preallocate array with nan - running an unbalanced ANOVA)
    A = [size(array_ofMonoc1,3),...
        size(array_ofMonoc2,3),...
        size(array_ofMonoc3,3),...
        size(array_ofMonoc4,3)];
    maxTrls = max(A);
    minTrls = min(A);
    tuned = nan(length(v1Ch),1);
    prefMonoc = nan(length(v1Ch),1);
    
    % FOR each individual unit 
    for i = 1:length(v1Ch)
        % preallocate y based on trial count
        y = nan(maxTrls,4);
        % create an array of each trial's median response for each monoc 
        % condition (trl x 4 monoc)
        y(1:size(array_ofMonoc1,3),1) = median(squeeze(array_ofMonoc1(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc2,3),2) = median(squeeze(array_ofMonoc2(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc3,3),3) = median(squeeze(array_ofMonoc3(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc4,3),4) = median(squeeze(array_ofMonoc4(200:450,i,:)),1)'; % median of each trial after stim onset
        % Now we perform baseline subtraction
        % First we get the blAvg for this contact across all trials
        y_bl(1,1) = median(median(squeeze(array_ofMonoc1(100:200,i,:)),1)');
        y_bl(1,2) = median(median(squeeze(array_ofMonoc2(100:200,i,:)),1)');
        y_bl(1,3) = median(median(squeeze(array_ofMonoc3(100:200,i,:)),1)');
        y_bl(1,4) = median(median(squeeze(array_ofMonoc4(100:200,i,:)),1)');
        y_blSub = y - median(y_bl);
        p = anova1(y_blSub,[],'off');
        if p < .05
            tuned(i,1) = true;
        else
            tuned(i,1) = false;
        end
        % Now we find the maximum response
        [M,maxRespIdx] = max(median(y_blSub,1,"omitmissing"));
        prefMonoc(i,1) = maxRespIdx;
        % And assign the null condition
        if maxRespIdx == 1
            nullRespIdx = 4;
        elseif maxRespIdx == 2
            nullRespIdx = 3;
        elseif maxRespIdx == 3
            nullRespIdx = 2;
        elseif maxRespIdx == 4
            nullRespIdx = 1;
        end
        nullMonoc(i,1) = nullRespIdx;
    end
   

    for i = 1:length(v1Ch)
        % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
        % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
        % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
        % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
        if prefMonoc(i,1) == 1
            prefCondOnFlash(i,1) = 12;
            nullCondOnFlash(i,1) = 11;
        elseif prefMonoc(i,1) == 2
            prefCondOnFlash(i,1) = 9;
            nullCondOnFlash(i,1) = 10;
        elseif prefMonoc(i,1) == 3
            prefCondOnFlash(i,1) = 10;
            nullCondOnFlash(i,1) = 9;
        elseif prefMonoc(i,1) == 4
            prefCondOnFlash(i,1) = 11;
            nullCondOnFlash(i,1) = 12;
        end

    end
   

    %% MUA in percent change
    % BRFS pref vs null
    % Get the number of trials for the chosen condition
    for i = 1:length(v1Ch)
        numTrials_ps = size(MUA_trials{penetration,1}{prefCondOnFlash(i,1),1},1);
        MUAflashOut_ps = nan(2001,numTrials_ps);
        for trl = 1:numTrials_ps
            MUAflashOut_ps(:,trl) = MUA_trials{penetration,1}{prefCondOnFlash(i,1),1}{trl,1}(:,v1Ch(i));
        end
        
        % median across trials
        ps_avg = median(MUAflashOut_ps,2,"omitmissing"); % input is (tm,trl)    
        % Calculate as Percent Change
        %              X(t) - avgBl
        % %Ch = 100 * -------------
        %                 avgBl
        psBl = median(ps_avg(1:200,:));
        if psBl == 0
            psBl = .1;
        end
        ps_PercentC = 100*((ps_avg-psBl)./psBl);
        averageMUAMatrix_BRFSps(:,i,penetration) = ps_PercentC;

    
        % Get the number of trials for the chosen condition
        numTrials_ns = size(MUA_trials{penetration,1}{nullCondOnFlash(i,1) ,1},1);
        MUAflashOut_ns = nan(2001,numTrials_ns);
        for trl = 1:numTrials_ns
            MUAflashOut_ns(:,trl) = MUA_trials{penetration,1}{nullCondOnFlash(i,1),1}{trl,1}(:,v1Ch(i));
        end
        % Average across trials and save output
        ns_avg = median(MUAflashOut_ns,2,"omitmissing"); % input is (tm,trl)    
        nsBl = median(ns_avg(1:200,:));
        if nsBl == 0
            nsBl = .1;
        end
        ns_PercentC = 100*((ns_avg-nsBl)./nsBl);  
        averageMUAMatrix_BRFSns(:,i,penetration) = ns_PercentC; % Average across trl. averageMUAMatrix is (tm x ch x x penetration)
    end
    disp(strcat('Done with file number: ',string(penetration)))
end


%% Organize into compartments, median and std across contacts+penetration
% reshape
useIdx = squeeze(~isnan(averageMUAMatrix_BRFSps(1,1,:))); 
ps_S = reshape(averageMUAMatrix_BRFSps(:,1:5,useIdx),[2001,135]);
ps_G = reshape(averageMUAMatrix_BRFSps(:,6:10,useIdx),[2001,135]);
ps_I = reshape(averageMUAMatrix_BRFSps(:,11:15,useIdx),[2001,135]);
ns_S = reshape(averageMUAMatrix_BRFSns(:,1:5,useIdx),[2001,135]);
ns_G = reshape(averageMUAMatrix_BRFSns(:,6:10,useIdx),[2001,135]);
ns_I = reshape(averageMUAMatrix_BRFSns(:,11:15,useIdx),[2001,135]);

% Average across penetrations
ps_S_avg = smoothdata(mean(ps_S,2,"omitmissing"),"gaussian",50);
ps_G_avg = smoothdata(mean(ps_G,2,"omitmissing"),"gaussian",50);
ps_I_avg = smoothdata(mean(ps_I,2,"omitmissing"),"gaussian",50);
ns_S_avg = smoothdata(mean(ns_S,2,"omitmissing"),"gaussian",50);
ns_G_avg = smoothdata(mean(ns_G,2,"omitmissing"),"gaussian",50);
ns_I_avg = smoothdata(mean(ns_I,2,"omitmissing"),"gaussian",50);
% % ps_S_avg = median(ps_S,2,"omitmissing");
% % ps_G_avg = median(ps_G,2,"omitmissing");
% % ps_I_avg = median(ps_I,2,"omitmissing");
% % ns_S_avg = median(ns_S,2,"omitmissing");
% % ns_G_avg = median(ns_G,2,"omitmissing");
% % ns_I_avg = median(ns_I,2,"omitmissing");

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
contactNum = size(ps_S,2);
ps_S_sem = std(ps_S,0,2,"omitmissing")./sqrt(contactNum); 
ps_G_sem = std(ps_G,0,2,"omitmissing")./sqrt(contactNum); 
ps_I_sem = std(ps_I,0,2,"omitmissing")./sqrt(contactNum); 
ns_S_sem = std(ns_S,0,2,"omitmissing")./sqrt(contactNum); 
ns_G_sem = std(ns_G,0,2,"omitmissing")./sqrt(contactNum); 
ns_I_sem = std(ns_I,0,2,"omitmissing")./sqrt(contactNum); 

% Mean +/- SEM
ps_S_avgPlusSEM = ps_S_avg + ps_S_sem; %pref stim -- median plus sem 1 
ps_S_avgMinusSEM = ps_S_avg - ps_S_sem; %pref stim -- median minus sem 1 
ps_G_avgPlusSEM = ps_G_avg + ps_G_sem; %pref stim -- median plus sem 1 
ps_G_avgMinusSEM = ps_G_avg - ps_G_sem; %pref stim -- median minus sem 1 
ps_I_avgPlusSEM = ps_I_avg + ps_I_sem; %pref stim -- median plus sem 1 
ps_I_avgMinusSEM = ps_I_avg - ps_I_sem; %pref stim -- median minus sem 1 

ns_S_avgPlusSEM = ns_S_avg + ns_S_sem; 
ns_S_avgMinusSEM = ns_S_avg - ns_S_sem;  
ns_G_avgPlusSEM = ns_G_avg + ns_G_sem; 
ns_G_avgMinusSEM = ns_G_avg - ns_G_sem; 
ns_I_avgPlusSEM = ns_I_avg + ns_I_sem; 
ns_I_avgMinusSEM = ns_I_avg - ns_I_sem; 

%% statistics
% Ok, the data is together for plotting, now lets run statistics on
% each laminar compartment to see if perceptual modulation occurs. The
% goal here is to run a t-test to see if the average response between
% 1200 and 1600ms significantly differs between ps and ns
% % % % useIdx = squeeze(~isnan(averageMUAMatrix_BRFSps(1,1,:))); 
% % % % tInput_ps_S_trans = reshape(squeeze(mean(averageMUAMatrix_BRFSps(1050:1200,1:5,useIdx),1)),[],1);
% % % % tInput_ps_G_trans = reshape(squeeze(mean(averageMUAMatrix_BRFSps(1050:1200,6:10,useIdx),1)),[],1);
% % % % tInput_ps_I_trans = reshape(squeeze(mean(averageMUAMatrix_BRFSps(1050:1200,11:15,useIdx),1)),[],1);
% % % % tInput_ns_S_trans = reshape(squeeze(mean(averageMUAMatrix_BRFSns(1050:1200,1:5,useIdx),1)),[],1);
% % % % tInput_ns_G_trans = reshape(squeeze(mean(averageMUAMatrix_BRFSns(1050:1200,6:10,useIdx),1)),[],1);
% % % % tInput_ns_I_trans = reshape(squeeze(mean(averageMUAMatrix_BRFSns(1050:1200,11:15,useIdx),1)),[],1);
% % % % tInput_ps_S_susta = reshape(squeeze(mean(averageMUAMatrix_BRFSps(1400:1801,1:5,useIdx),1)),[],1);
% % % % tInput_ps_G_susta = reshape(squeeze(mean(averageMUAMatrix_BRFSps(1400:1801,6:10,useIdx),1)),[],1);
% % % % tInput_ps_I_susta = reshape(squeeze(mean(averageMUAMatrix_BRFSps(1400:1801,11:15,useIdx),1)),[],1);
% % % % tInput_ns_S_susta = reshape(squeeze(mean(averageMUAMatrix_BRFSns(1400:1801,1:5,useIdx),1)),[],1);
% % % % tInput_ns_G_susta = reshape(squeeze(mean(averageMUAMatrix_BRFSns(1400:1801,6:10,useIdx),1)),[],1);
% % % % tInput_ns_I_susta = reshape(squeeze(mean(averageMUAMatrix_BRFSns(1400:1801,11:15,useIdx),1)),[],1);
% % % % 
% % % % [h_S_trans,p_S_trans] = ttest2(tInput_ps_S_trans,tInput_ns_S_trans);
% % % % [h_G_trans,p_G_trans] = ttest2(tInput_ps_G_trans,tInput_ns_G_trans);
% % % % [h_I_trans,p_I_trans] = ttest2(tInput_ps_I_trans,tInput_ns_I_trans);
% % % % [h_S_susta,p_S_susta] = ttest2(tInput_ps_S_susta,tInput_ns_S_susta);
% % % % [h_G_susta,p_G_susta] = ttest2(tInput_ps_G_susta,tInput_ns_G_susta);
% % % % [h_I_susta,p_I_susta] = ttest2(tInput_ps_I_susta,tInput_ns_I_susta);


%% Figure generation! 
% tiledLayout plot
% close all
tm_full = -200:1800; % 1801 total timepoints
lamCom = figure;
set(gcf,"Position",[1000 123.6667 757.6667 1.1140e+03])
t = tiledlayout(3,1);
nexttile
    plot(tm_full,ps_S_avg,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_S_avgPlusSEM,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_S_avgMinusSEM,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_S_avg,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_S_avgPlusSEM,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_S_avgMinusSEM,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(833)
    vline(1633)
    % xregion(850,1000)
    % xregion(1200,1600)
    title('Supragranular')
    set(gca,'xtick',[])
    box("off")
    legend('Preferred stimulus','Null stimulus')
nexttile
    plot(tm_full,ps_G_avg,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_G_avgPlusSEM,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_G_avgMinusSEM,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_G_avg,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_G_avgPlusSEM,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_G_avgMinusSEM,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    ylim([0 40])
    vline(0)
    vline(833)
    vline(1633)
    % xregion(850,1000)
    % xregion(1200,1600)
    title('Granular')
    set(gca,'xtick',[])
    box("off")
nexttile
    plot(tm_full,ps_I_avg,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_I_avgPlusSEM,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_I_avgMinusSEM,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_I_avg,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_I_avgPlusSEM,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_I_avgMinusSEM,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':');
    ylim([0 40])
    vline(0)
    vline(833)
    vline(1633)
    % xregion(850,1000)
    % xregion(1200,1600)
    ylabel({'% change from baseline'})
    xlabel('Time (ms)')
    title('Infragranular')
    box("off")

titleText = {'MEAN MUA of 135 electrodes (27 penetrations) per laminar compartment'};
title(t,titleText,'Interpreter','none')


%save fig
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
       disp('alright, saving figure to plotdir')
        cd(plotDir)
        figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.png');
        saveas(lamCom,figName_lamCom)
        figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.svg');
        saveas(lamCom,figName_lamCom)
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end


