%% bmcBRFS_CSD_grandAverage_loadMAT
% Are there observable differences between trial-types with LFP CSD?
% initialize variables
clearvars -except LFP_trials
close all


%% For loop
plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir\CSDFigs';
dataDir = 'S:\TrialTriggeredLFPandMUA';
% % dataDir = 'D:\sortedData_240229';
% % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
% % tic %takes just over 2 min
if ~exist('LFP_trials','var')
    load('LFP_trials.mat') % format is LFP_trials{1,penetration}{cond,1}{trial,flash}
end
% % toc
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

averageCSDMatrix_diopDichop = nan(1201,15,2,size(officLamAssign,1));
averageCSDMatrix_BRFS = nan(1201,15,2,size(officLamAssign,1));


for penetration = 1:size(LFP_trials,2)
    
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
        for trl = 1:size(LFP_trials{1,penetration}{cond,1},1)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(LFP_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:size(LFP_trials{1,penetration}{cond,1},1)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(LFP_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:size(LFP_trials{1,penetration}{cond,1},1)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(LFP_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:size(LFP_trials{1,penetration}{cond,1},1)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(LFP_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); 
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
        
    overallPref = mode(prefCondOnFlash);
    if overallPref == 9
        overallNull = 10;
    elseif overallPref == 10
        overallNull = 9;
    elseif overallPref == 11
        overallNull = 12;
    elseif overallPref == 12
        overallNull = 11;
    end

    %% CSD matrix across channels
    chForCSDcalc = v1Ch(1)-1:v1Ch(end)+1;
    % Dioptic vs dichoptic
    count = 0;
    for conditionNumber = [1 3]        
        % Get the number of trials for the chosen condition
        numTrials = size(LFP_trials{1,penetration}{conditionNumber,1},1);
        CSDbinocOut = nan(1201,length(v1Ch),numTrials);
        for trl = 1:numTrials
            LFPbinocTrl = LFP_trials{1,penetration}{conditionNumber,1}{trl,1}(:,chForCSDcalc); % adding ch above and below for CSD calc (to use as vaknin pad)
            blLFPbinoc = mean(LFPbinocTrl(100:200,:),1);
            LFPbinocBlSub = LFPbinocTrl-blLFPbinoc;
            CSDbinocTrl = calcCSD_classic(LFPbinocBlSub);
            CSDbinocOut(:,:,trl) = CSDbinocTrl(:,2:16); % limit to origional V1 ch lim
        end
        % Average across trials and save output
        count = count + 1; % for pref vs null
        averageCSDMatrix_diopDichop(:,:,count,penetration) = mean(CSDbinocOut,3); % Average across trl. averageCSDMatrix is (ch1 x ch2 x cond x penetration)
    end

    % BRFS pref vs null
    count = 0;
    for conditionNumber = [overallPref overallNull]        
        % Get the number of trials for the chosen condition
        numTrials = size(LFP_trials{1,penetration}{conditionNumber,1},1);
        CSDflashOut = nan(1201,length(v1Ch),numTrials);
        for trl = 1:numTrials
            LFPflashTrl = LFP_trials{1,penetration}{conditionNumber,1}{trl,2}(:,chForCSDcalc);
            blLFPflash = mean(LFPflashTrl(100:200,:),1);
            LFPflashBlSub = LFPflashTrl-blLFPflash;
            CSDflashTrl = calcCSD_classic(LFPflashBlSub);
            CSDflashOut(:,:,trl) = CSDflashTrl(:,2:16);
        end
        % Average across trials and save output
        count = count + 1; % for pref vs null
        averageCSDMatrix_BRFS(:,:,count,penetration) = mean(CSDflashOut,3); % Average across trl. averageCSDMatrix is (ch1 x ch2 x cond x penetration)
    end


    disp(strcat('Done with file number: ',string(penetration)))
end


%% plot dioptic vs dichoptic
grandAverageCSD_diopDichop = mean(averageCSDMatrix_diopDichop,4,"omitmissing"); % average across penetration


f = figure;
set(f,"Position",[345.6667 256.3333 1698 981.3333])
% Visualize the CSD matrix
% dioptic
ax(1) = subplot(2,5,1);
imagesc(-200:1000,1:15,grandAverageCSD_diopDichop(:,:,1)');
oldcmap = colormap(ax(1),'jet');
newcmap = colormap(flipud(oldcmap));
colormap(ax(1),newcmap)
clim([-3000 3000]);
xlabel('Time (ms)');
ylabel('Channel');
cb = colorbar(); 
ylabel(cb,'(nA/mm)^3','FontSize',12)
xl = xline(0,'--','Stimulus onset','LineWidth',3);
title('Dioptic');

%dichoptic
ax(2) = subplot(2,5,2);
imagesc(-200:1000,1:15,grandAverageCSD_diopDichop(:,:,2)');
colormap(ax(2),newcmap)
clim([-3000 3000]);
xlabel('Time (ms)');
cb = colorbar(); 
xl = xline(0,'--','Stimulus onset','LineWidth',3);
title('Dichoptic');

% Visualize the raw diff
differenceMatrix_diopDichoip = ...
    abs(grandAverageCSD_diopDichop(:,:,1)')-abs(grandAverageCSD_diopDichop(:,:,2)');
ax(3) = subplot(2,5,3);
imagesc(-200:1000,1:15,differenceMatrix_diopDichoip);
colormap(ax(3),'bone');
clim([-1500 1500]);
xlabel('Time (ms)');
cb = colorbar(); 
xl = xline(0,'--','Stimulus onset','LineWidth',3);
title('|Difference|: |dioptic-dichoptic|');



% Difference plot - with tStat
usePenetration = [1:4 6 8:25 27:29 31];
h = nan(15,24);
p = nan(15,24);
ci = nan(15,24,2);
for ch = 1:15
    CSDMatrix1 = squeeze(averageCSDMatrix_diopDichop(:,ch,1,usePenetration)); 
    CSDMatrix2 = squeeze(averageCSDMatrix_diopDichop(:,ch,2,usePenetration)); 
    for bin = 1:24
        tm = [1:50]+(50*(bin-1));
        CSD_binned1(bin,:) = abs(mean(CSDMatrix1(tm,:),1));   % output should be 24xlength(usePenetration)
        CSD_binned2(bin,:) = abs(mean(CSDMatrix2(tm,:),1));
        [h(ch,bin),p(ch,bin),ci(ch,bin,:),stats(ch,bin)] =...
            ttest2(CSD_binned1(bin,:),CSD_binned2(bin,:),'Alpha',0.05,'Tail','left');
        tStat_1(ch,bin) = stats(ch,bin).tstat;

    end
end
ax(4) = subplot(2,5,4);
imagesc(tStat_1);
colormap(ax(4),'bone');
cb = colorbar(); 
xl = xline(3.5,'--','Stimulus onset','LineWidth',3);
xlabel('Index of 50ms bin');
title('left-tail tScoreMap of |difference|');

% hypothesis test logical
ax(5) = subplot(2,5,5);
imagesc(h);
colormap(ax(5),'bone');
e = colorbar(); 
e.Label.String = "h = 1 or 0";
e.Label.Rotation = 270;
xl = xline(3.5,'--','Stimulus onset','LineWidth',3);
xlabel('Index of 50ms bin');
title('|Dichoptic > Dioptic| Hypothesis test');




%% plot BRFS
grandAverageCSD_BRFS = median(averageCSDMatrix_BRFS,4,"omitmissing"); % average across penetration

% Visualize the CSD matrix
ax(6) = subplot(2,5,6);
imagesc(-200:1000,1:15,grandAverageCSD_BRFS(:,:,1)');
colormap(ax(6),newcmap)
clim([-1500 1500]);
xlabel('Time (ms)');
ylabel('channel');
cb = colorbar(); 
ylabel(cb,'(nA/mm)^3','FontSize',12)
xl = xline(0,'--','Stimulus onset','LineWidth',3);
xl = xline(800,'--','Stimulus offset','LineWidth',3);
title('Preferred stimulus BRFS flash');


ax(7) = subplot(2,5,7);
imagesc(-200:1000,1:15,grandAverageCSD_BRFS(:,:,2)');
colormap(ax(7),newcmap)
clim([-1500 1500]);
xlabel('Time (ms)');
ylabel('channel');
cb = colorbar(); 
ylabel(cb,'(nA/mm)^3','FontSize',12)
xl = xline(0,'--','Stimulus onset','LineWidth',3);
xl = xline(800,'--','Stimulus offset','LineWidth',3);
title('Non-preferred stimulus BRFS flash');

% Visualize the raw diff
differenceMatrix_diopDichoip = ...
    abs(grandAverageCSD_BRFS(:,:,1)')-abs(grandAverageCSD_BRFS(:,:,2)');
ax(8) = subplot(2,5,8);
imagesc(-200:1000,1:15,differenceMatrix_diopDichoip);
colormap(ax(8),'bone');
clim([-1000 1000]);
xlabel('Time (ms)');
cb = colorbar(); 
xl = xline(0,'--','Stimulus onset','LineWidth',3);
title('|Difference|: |PS-NS|');



% Difference plot - with tStat
usePenetration = [1:4 6 8:25 27:29 31];
h = nan(15,24);
p = nan(15,24);
ci = nan(15,24,2);
for ch = 1:15
    CSDMatrix1 = squeeze(averageCSDMatrix_BRFS(:,ch,1,usePenetration)); 
    CSDMatrix2 = squeeze(averageCSDMatrix_BRFS(:,ch,2,usePenetration)); 
    for bin = 1:24
        tm = [1:50]+(50*(bin-1));
        CSD_binned1(bin,:) = abs(mean(CSDMatrix1(tm,:),1));   % output should be 24xlength(usePenetration)
        CSD_binned2(bin,:) = abs(mean(CSDMatrix2(tm,:),1));
        [h(ch,bin),p(ch,bin),ci(ch,bin,:),stats(ch,bin)] =...
            ttest2(CSD_binned1(bin,:),CSD_binned2(bin,:),'Alpha',0.05,'Tail','right');
        tStat_1(ch,bin) = stats(ch,bin).tstat;

    end
end
ax(9) = subplot(2,5,9);
imagesc(tStat_1);
colormap(ax(9),'bone');
cb = colorbar(); 
xl = xline(3.5,'--','Stimulus onset','LineWidth',3);
xlabel('Index of 50ms bin');
title('right-tail tScoreMap of |difference|');

% hypothesis test logical
ax(10) = subplot(2,5,10);
imagesc(h);
colormap(ax(10),'bone');
e = colorbar(); 
e.Label.String = "h = 1 or 0";
e.Label.Rotation = 270;
xl = xline(3.5,'--','Stimulus onset','LineWidth',3);
xlabel('Index of 50ms bin');
title('|PS > NS| Hypothesis test');



sgtitle('CSD Penetration Average')


%% Save output
cd(plotDir)
saveName = strcat('CSDPenetrationAvg_prefFromLFP.png');
saveas(f,saveName) 
saveName = strcat('CSDPenetrationAvg_prefFromLFP.png');
saveas(f,saveName) 
