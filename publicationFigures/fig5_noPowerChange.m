%% bmcBRFS_coherenceMatrix
% Are there observable differences between trial-types with LFP coherence?
% initialize variables
clearvars -except LFP_trials
tic
workingPC = 'home'; % options: 'home', 'office'

%% Setup
disp('start time')
datetime
if strcmp(workingPC,'home')
    codeDir     = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)


cd(dataDir)
if ~exist('LFP_trials','var')
    load('LFP_trials.mat') % format is LFP_trials{penetration,1}{cond,1}{trial,1}(2001tm x 65Ch)
end

averageCoherenceMatrix_diopDichop = nan(15,15,2,size(officLamAssign,1));
averageCoherenceMatrix_BRFS = nan(15,15,2,size(officLamAssign,1));

%% for loop
for penetration = 1:size(LFP_trials,1)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    penetrationFileName = strcat('trialTriggeredData_',probeName(1:19),'.mat');
    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',penetrationFileName))
        continue
    end
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
    if strcmp(string(officLamAssign.ChToUse(penetration)),"1:32")
        if any(v1Ch > 32) || any(v1Ch < 1)
            warning('skipping session without full column for now')
            disp(penetrationFileName)
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('skipping session without full column for now')
            disp(penetrationFileName)
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
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
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

    %% Coherence matrix across channels
    tm_coher = 1289:1800; % Time window of data. Last 512ms of trial. 


    %% PSD
    % Compute average PSD for each conditions 
    count = 0;
    for cond = [1 3 overallPref overallNull]
        count = count + 1; % for condition index
        % format is LFP_trials{penetration,1}{cond,1}{trial}(tm,ch)
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1);
            DATA_PSDin = LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch);
            % Compute PSD
            [PSD, freq] = calcPSD(DATA_PSDin,tm_coher); %input is tm x ch, output is freq x ch 
            % Store PSD across trials
            PSD_trls(:,:,trl) = PSD;  %(freq x ch x trl)
        end

        % average power across. trials
        powerAvg = mean(PSD_trls,3,"omitnan");       
        
        % normalize power. At each frequency, we divide the power of each channel
        % by that of the channel with the highest power.
        PSD_norm = nan(size(powerAvg)); 
        for freq = 1:size(powerAvg,1)
            powerOfAllChAtThisFreq = powerAvg(freq,:);
            [maxPower,chWithMaxPower] = max(powerOfAllChAtThisFreq);
            PSD_norm(freq,:) = powerOfAllChAtThisFreq ./ maxPower; 
            % output of power_norm is (freq,ch)
        end
        % Store output for condition and penetration
        if cond == 1
            PSD_normCond.dioptic(:,:,penetration) = PSD_norm;
        elseif cond == 3
            PSD_normCond.dichoptic(:,:,penetration) = PSD_norm;
        elseif cond == overallPref
            PSD_normCond.BRFSps(:,:,penetration) = PSD_norm;
        elseif cond == overallNull
            PSD_normCond.BRFSns(:,:,penetration) = PSD_norm;
        end

    end
        disp(strcat('Done with file number: ',string(penetration)))

end %% END OF penetration = 1:size(LFP_trials,1).





%% PSD 
%%%%%%%%%%%%%
%%% POWER %%%
%%%%%%%%%%%%%

%% Average variables across penetration
PSD_grandAvg.dioptic = median(PSD_normCond.dioptic,3,"omitmissing"); % average across penetration
PSD_grandAvg.dichoptic = median(PSD_normCond.dichoptic,3,"omitmissing"); % average across penetration
PSD_grandAvg.BRFSps = median(PSD_normCond.BRFSps,3,"omitmissing"); % average across penetration
PSD_grandAvg.BRFSns = median(PSD_normCond.BRFSns,3,"omitmissing"); % average across penetration

%% Plot PSD

% Visualize the coherence matrix
f = figure;
% set(f,"Position",[-1238 -55 1126 641])
ax(1) = subplot(2,3,1);
imagesc(PSD_grandAvg.dioptic(1:50,:)'); % expected input is ch x freq
colormap(ax(1),'parula');
colorbar;
% clim([.75 1])
xlabel('Frequency');
ylabel('Channel');
title('Dioptic');

ax(2) = subplot(2,3,2);
imagesc(PSD_grandAvg.dichoptic(1:50,:)');
colormap(ax(2),'parula');
colorbar;
% clim([.75 1])
xlabel('Frequency');
ylabel('Channel');
title('Dichoptic');


% Difference plot 
diff_1 = PSD_grandAvg.dioptic(1:50,:)'-PSD_grandAvg.dichoptic(1:50,:)';
ax(3) = subplot(2,3,3);
imagesc(diff_1);
colormap(ax(3),'bone');
clim([0 .3])
e = colorbar;
e.Label.String = "PSD Difference";
e.Label.Rotation = 270;
xlabel('Frequency');
ylabel('Channel');
title('difference');






%% plot BRFS
% Visualize the coherence matrix
% set(f,"Position",[-1238 -55 1126 641])
ax(1) = subplot(2,3,4);
imagesc(PSD_grandAvg.BRFSps(1:50,:)'); % expected input is ch x freq
colormap(ax(1),'parula');
colorbar;
% clim([.75 1])
xlabel('Frequency');
ylabel('Channel');
title('BRFS preferred stimulus flash');

ax(2) = subplot(2,3,5);
imagesc(PSD_grandAvg.BRFSns(1:50,:)');
colormap(ax(2),'parula');
colorbar;
% clim([.75 1])
xlabel('Frequency');
ylabel('Channel');
title('BRFS null stimulus flash');


% Difference plot 
diff_2 = PSD_grandAvg.BRFSps(1:50,:)'-PSD_grandAvg.BRFSns(1:50,:)';
ax(3) = subplot(2,3,6);
imagesc(diff_2);
colormap(ax(3),'bone');
clim([0 .3])
e = colorbar;
e.Label.String = "PSD Difference";
e.Label.Rotation = 270;
xlabel('Frequency');
ylabel('Channel');
title('difference');



%% Save output
save fig
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
        disp('alright, saving figure to plotdir')
        sgtitle('Coherence Penetration Average')
        cd(plotDir)
        saveName = 'fig5_noPowerChange.png';
        saveas(f,saveName) 
        saveName = 'fig5_noPowerChange.svg';
        saveas(f,saveName) 
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end




