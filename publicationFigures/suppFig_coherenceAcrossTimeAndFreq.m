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
    plotDir = 'C:\Users\Brock Carlson\Documents\Coherence manuscript - home pc\officialFigureOutputsFromMATLAB';
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
    % Parameters for mscohere
    fs = 1000;  % Sampling frequency in Hz
    windowSize = 256;  % Window size for computing the coherence
    overlap = windowSize/2;  % Overlap between windows
    % tm_full = 1:2001;
    tm_bl = 1:200;
    % tm_1stOnset = 201:1000; %800ms adapter
    % tm_2ndOnset = 1001:1800;
    % tm_offset = 1801:2001;
    % tm_coher = 1289:1800; % Time window of data. Last 512ms of trial. 

    params = struct( ...
        'name', {'early_0_30','alpha','beta','gamma'}, ...
        'tm_coher', { ...
            1001:1513, ...     % early window
            1289:1800, ...     % original window
            1289:1800, ...
            1289:1800 }, ...
        'freq_idx', { ...
            2:9, ...           % 0–30 Hz
            3:4, ...           % alpha
            5:8, ...           % beta
            9:16 });           % gamma

    % BRFS pref vs null
    for p = 1:length(params)
        tm_coher = params(p).tm_coher;
        f_idx    = params(p).freq_idx;
        count = 0;

        for conditionNumber = [overallPref overallNull]        
            % Get the number of trials for the chosen condition
            numTrials = size(LFP_trials{penetration,1}{conditionNumber,1},1);
            % Initialize coherence matrix
            coherenceMatrix = nan(15,15, numTrials);
            % Loop through all trials and compute coherence for each channel pair
            for trialNumber = 1:numTrials
                for channel1 = v1Ch
                    for channel2 = v1Ch(1):channel1
                        % Extract the LFP_bb data for the chosen channels and trial
                        % format is LFP_trials{penetration,1}{cond,1}{trial,flash}
                        lfpGammaData1 = LFP_trials{penetration,1}{conditionNumber,1}{trialNumber,1}(:,channel1);
                        lfpGammaData2 = LFP_trials{penetration,1}{conditionNumber,1}{trialNumber,1}(:,channel2);
                        %baseline subtract
                        bl1 = median(lfpGammaData1(tm_bl));
                        lfp_blsub_1 = lfpGammaData1 - bl1;
                        bl2 = median(lfpGammaData2(tm_bl));
                        lfp_blsub_2 = lfpGammaData2 - bl2;
        
                        % Compute coherence
                        [coherence, freq] = mscohere(lfp_blsub_1(tm_coher), lfp_blsub_2(tm_coher), windowSize, overlap, [], fs);
            
                        % Store coherence in the matrix. Indexing at the
                        % front allows v1channels of different number to be
                        % set from 1-15
                        coherenceMatrix(channel1-v1Ch(1)+1, channel2-v1Ch(1)+1, trialNumber) = median(coherence(f_idx));  % coherenceMatrix is (ch1 x ch2 x trialNum)
                    end
                end
            end
            % Average across trials and save output
            count = count + 1; % for pref vs null
            averageCoherenceMatrix_BRFS(:,:,count,p,penetration) = median(coherenceMatrix,3); 
            % dims: ch x ch x cond(2) x param(4) x penetration
        end
        disp(strcat('Done with file number: ',string(penetration)))
    end 
end

%% Grand Average over penetration
grandAverageCoherence_BRFS = median(averageCoherenceMatrix_BRFS,5,"omitmissing");
% dims: ch x ch x condition(2) x param(4)


%% plot BRFS
for p = 1:length(params)

    % Extract this param’s data
    grandMat = squeeze(grandAverageCoherence_BRFS(:,:,:,p));
    % dims now: ch x ch x condition(2)

    % --- FIGURE SETUP ---
    f = figure;
    set(f,"Position",[113 458 1699 420])

    %% --- Subplot 1: Preferred ---
    ax(1) = subplot(1,4,1);
    imagesc(grandMat(:,:,1));
    colormap(ax(1),'jet');
    colorbar;
    xlabel('Channel');
    ylabel('Channel');
    title('Preferred Stimulus BRFS flash');

    %% --- Subplot 2: Non-preferred ---
    ax(2) = subplot(1,4,2);
    imagesc(grandMat(:,:,2));
    colormap(ax(2),'jet');
    colorbar;
    xlabel('Channel');
    ylabel('Channel');
    title('Non-preferred stimulus BRFS flash');

    %% --- Subplot 3: Difference ---
    % Pull penetration-level data for stats
    coherenceMatrix1 = squeeze(averageCoherenceMatrix_BRFS(:,:,1,p,:));
    coherenceMatrix2 = squeeze(averageCoherenceMatrix_BRFS(:,:,2,p,:));

    diff_2 = grandMat(:,:,1) - grandMat(:,:,2);

    ax(3) = subplot(1,4,3);
    imagesc(diff_2);
    hline(5.5); hline(10.5);
    vline(5.5); vline(10.5);
    colormap(ax(3),'bone');
    clim([-.04 .04])
    e = colorbar;
    e.Label.String = "Coherence Difference";
    e.Label.Rotation = 270;
    xlabel('Channel');
    ylabel('Channel');
    title('difference');

    %% --- STATISTICS (unchanged logic) ---
    SxS_2 = diff_2(1:5,1:5);
    GxS_2 = diff_2(6:10,1:5);
    IxS_2 = diff_2(11:15,1:5);
    GxG_2 = diff_2(6:10,6:10);
    IxG_2 = diff_2(11:15,6:10);
    IxI_2 = diff_2(11:15,11:15);

    holder_cross_2(:,1) = reshape(GxS_2,[25,1]);
    holder_cross_2(:,2) = reshape(IxS_2,[25,1]);
    holder_cross_2(:,3) = reshape(IxG_2,[25,1]);

    [aov_cross_2{p}, tbl{p}, stats{p}] = anova1(holder_cross_2,[],"off");

    % same-compartment
    SxS_2(SxS_2 == 0) = NaN;
    GxG_2(GxG_2 == 0) = NaN;
    IxI_2(IxI_2 == 0) = NaN;

    holder_same_2(:,1) = reshape(SxS_2,[25,1]);
    holder_same_2(:,2) = reshape(GxG_2,[25,1]);
    holder_same_2(:,3) = reshape(IxI_2,[25,1]);

    aov_same_2{p} = anova1(holder_same_2,[],"off");

    %% --- Subplot 4: Bar plot ---
    means_cross = nanmean(holder_cross_2);
    sems_cross = nanstd(holder_cross_2) ./ sqrt(sum(~isnan(holder_cross_2)));

    ax(4) = subplot(1,4,4);
    b = bar(means_cross, 'FaceColor', [0.2 0.2 0.8]);
    hold on;

    errorbar(1:3, means_cross, sems_cross, 'k', ...
        'LineStyle','none','LineWidth',1.5);

    x1 = 1; x2 = 3;
    y = max(means_cross + sems_cross) * 1.1;
    plot([x1, x2], [y, y], 'k-', 'LineWidth', 1.5);
    text(mean([x1, x2]), y * 1.05, '*', ...
        'HorizontalAlignment','center','FontSize',16);

    title('Cross-Compartment Comparisons');
    xticks(1:3);
    xticklabels({'GxS','IxS','IxG'});
    ylabel('Mean Value');
    xlabel('Group');
    grid on;
    hold off;

    %% --- Title ---
    sgtitle(sprintf('Coherence Penetration Average | %s', params(p).name))

    %% --- Save (parameterized) ---
    answer = questdlg('Would you like to save this figure?', 'Y', 'N');
    switch answer
        case 'Yes'
            cd(plotDir)
            baseName = strcat('suppFig_coherenceAcrossTimeAndFreq_', params(p).name);
            saveas(f, strcat(baseName,'.png'))
            saveas(f, strcat(baseName,'.svg'))
        case 'No'
            cd(plotDir)
    end

end

toc