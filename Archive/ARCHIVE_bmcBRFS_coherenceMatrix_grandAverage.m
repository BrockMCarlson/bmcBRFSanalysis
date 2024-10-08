%% bmcBRFS_coherenceMatrix
% Are there observable differences between trial-types with LFP coherence?
% initialize variables
clear
close all
codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures\grand average');
cd(codeDir)
dataDir = 'D:\sortedData_240229';
% dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\sortedData_240229';
% % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
if ~exist('LFP_trials','var')
    load('LFP_trials.mat') % format is LFP_trials{1,penetration}{cond,1}{trial,flash}
end
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

averageCoherenceMatrix = nan(15,15,2,size(officLamAssign,1));

%% For loop
for penetration = 1:size(LFP_trials,2)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    fileName = strcat('sortedData_',probeName(1:19));
    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',fileName))
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
            disp(fileName)
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('skipping session without full column for now')
            disp(fileName)
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

    %% Coherence matrix across channels
    % Choose the condition
    close all
    count = 0;
    % for conditionNumber = [overallPref overallNull]
    for conditionNumber = [1 3] % Dioptic vs dichoptic
        
        % Get the number of trials for the chosen condition
        %LFP_trials{1 x 31penetrations}{20conditions x 1}{trialNum x 2}{1201tm x 32Ch};
        numTrials = size(LFP_trials{1,penetration}{conditionNumber,1}, 1);
        
        % Parameters for mscohere
        fs = 1000;  % Sampling frequency in Hz
        windowSize = 256;  % Window size for computing the coherence
        overlap = windowSize/2;  % Overlap between windows
        tm_coher = 490:1001; % Time window of data. Last 512ms of trial. 
        tm = 1:1001; % Time window of data. Omitting any data after 1001ms to avoid offset transient
        
        % Initialize coherence matrix
        coherenceMatrix = nan(15,15, numTrials);
        % Loop through all trials and compute coherence for each channel pair
        for trialNumber = 1:numTrials
            for channel1 = v1Ch
                for channel2 = v1Ch(1):channel1
                    % Extract the LFP_bb data for the chosen channels and trial
                    lfpGammaData1 = LFP_trials{1,penetration}{conditionNumber,1}{trialNumber,2}(tm,channel1);
                    lfpGammaData2 = LFP_trials{1,penetration}{conditionNumber,1}{trialNumber,2}(tm,channel2);
                    %baseline subtract
                    bl1 = median(lfpGammaData1(1:200));
                    lfp_blsub_1 = lfpGammaData1 - bl1;
                    bl2 = median(lfpGammaData2(1:200));
                    lfp_blsub_2 = lfpGammaData2 - bl2;
    
                    % Compute coherence
                    [coherence, freq] = mscohere(lfp_blsub_1(tm_coher), lfp_blsub_2(tm_coher), windowSize, overlap, [], fs);
        
                    % Store coherence in the matrix
                    coherenceMatrix(channel1-v1Ch(1)+1, channel2-v1Ch(1)+1, trialNumber) = median(coherence(2:8));  % You can use median or any other aggregation method
                end
            end
        end
        
        % Average across trials and save output
        count = count + 1; % for pref vs null
        averageCoherenceMatrix(:,:,count,penetration) = median(coherenceMatrix,3);
        
    end
    disp(strcat('Done with file number: ',string(penetration)))
end
%% save LFP_trials output


%% plot
grandAverageCoherence = median(averageCoherenceMatrix,4,"omitmissing");

disp('analyzing freq')
disp(freq(2:8))


% Visualize the coherence matrix
figure;
imagesc(grandAverageCoherence(:,:,1));
colormap('jet');
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Preferred Stimulus BRFS flash Grand Average');

figure;
imagesc(grandAverageCoherence(:,:,2));
colormap('jet');
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Non-preferred stimulus BRFS flash Grand Average');

% Difference plot of coherence matrix
% abs(PS - NS)
% 9 NPO --> PO; 10 PO-->NPO
% abs(cond9-cond10)
% abs(cond12-cond11)
differenceMatrix = nan(15,15);
differenceMatrix(:,:) = grandAverageCoherence(:,:,1) - grandAverageCoherence(:,:,2);

% Visualize the diff  matrix
figure;
imagesc(differenceMatrix(:,:));
colormap('bone');
% clim([.5 .75])
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Diff of Pref vs Null');


%% Difference plot
coherenceMatrix1 = squeeze(averageCoherenceMatrix(:,:,1,:));
coherenceMatrix2 = squeeze(averageCoherenceMatrix(:,:,2,:));
for ch1 = 1:15
    for ch2 = 1:ch1
        [h(ch1,ch2),p(ch1,ch2),ci(ch1,ch2,:),stats(ch1,ch2)]...
            = ttest2(...
            coherenceMatrix1(ch1,ch2,:),...
            coherenceMatrix2(ch1,ch2,:));
        tStat(ch1,ch2) = stats(ch1,ch2).tstat;
    end
end

% Visualize the coherence matrix
figure
ax(1) = subplot(1,2,1);
imagesc(tStat);
colormap(ax(1),'bone');
e = colorbar;
e.Label.String = "tStat";
e.Label.Rotation = 270;
xlabel('Channel');
ylabel('Channel');
title('Dioptic vs Dichoptic tScoreMap');

% Visualize the coherence matrix
ax(2) = subplot(1,2,2);
imagesc(h);
colormap(ax(2),'bone');
e = colorbar;
e.Label.String = "h = 1 or 0";
e.Label.Rotation = 270;
xlabel('Channel');
ylabel('Channel');
title('Dioptic vs Dichoptic Hypothesis test');
% sgtitle(probeName(1:end-1),"Interpreter","none")
% % cd(plotDir)
% % saveName = strcat('coherence_',probeName(1:end-1),'.png');
% % saveas(f,saveName) 
% % close all





%% Compare laminar compartments
% % % 
% % % clearvars -except IDX
% % % % Default parameter values
% % % windowSize = 256;  % Window size for computing the coherence
% % % overlap = windowSize/2;  % Overlap between windows
% % % nfft = 512;  % Number of points to use in the FFT
% % % fs = 1000;  % Sampling frequency in Hz
% % % 
% % % % Define channel ranges for different layers
% % % supragranularChannels = 3:7;
% % % granularChannels = 8:12;
% % % infragranularChannels = 13:17;
% % % 
% % % % Frequency range of interest
% % % freqRange = [1 100];
% % % 
% % % % Loop through condition numbers
% % % for conditionNumber = [12, 11]
% % %     % Get the number of trials for the chosen condition
% % %     numTrials = size(IDX(conditionNumber).LFP_bb, 1);
% % % 
% % %     % Initialize freqLimValues outside the loop
% % %     [~, freqLimValues] = mscohere(zeros(windowSize, 1), zeros(windowSize, 1), windowSize, overlap, nfft, fs);
% % % 
% % %     % Initialize arrays to store coherence data for each layer pair
% % %     coherenceGranularSupragranular = zeros(numTrials, length(freqLimValues));
% % %     coherenceGranularInfragranular = zeros(numTrials, length(freqLimValues));
% % %     coherenceSupragranularInfragranular = zeros(numTrials, length(freqLimValues));
% % % 
% % %     % Loop through all trials and compute coherence for each layer pair
% % %     for trialNumber = 1:numTrials
% % %         % Extract the LFP_bb data for each layer
% % %         lfpSupragranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, supragranularChannels);
% % %         lfpGranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, granularChannels);
% % %         lfpInfragranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, infragranularChannels);
% % % 
% % %         % Compute coherence for each layer pair
% % %         [coherenceGranularSupragranularTrial, freq] = mscohere(lfpGranular, lfpSupragranular, windowSize, overlap, nfft, fs);
% % %         [coherenceGranularInfragranularTrial, ~] = mscohere(lfpGranular, lfpInfragranular, windowSize, overlap, nfft, fs);
% % %         [coherenceSupragranularInfragranularTrial, ~] = mscohere(lfpSupragranular, lfpInfragranular, windowSize, overlap, nfft, fs);
% % % 
% % %         % Extract frequency values within the specified range
% % %         freqIndices = find(freq >= freqRange(1) & freq <= freqRange(2));
% % %         freqLimValues = freq(freq >= freqRange(1) & freq <= freqRange(2));
% % % 
% % %         % Store median coherence data for each layer pair
% % %         coherenceGranularSupragranular(trialNumber, :) = median(coherenceGranularSupragranularTrial, 2);
% % %         coherenceGranularInfragranular(trialNumber, :) = median(coherenceGranularInfragranularTrial, 2);
% % %         coherenceSupragranularInfragranular(trialNumber, :) = median(coherenceSupragranularInfragranularTrial, 2);
% % % 
% % %     end
% % % 
% % %     % Plot the average coherence across all trials for each layer pair on a log-scaled x-axis
% % %     figure;
% % %     semilogx(freq, median(coherenceGranularSupragranular, 1), 'DisplayName', 'Granular to Supragranular');
% % %     hold on;
% % %     semilogx(freq, median(coherenceGranularInfragranular, 1), 'DisplayName', 'Granular to Infragranular');
% % %     semilogx(freq, median(coherenceSupragranularInfragranular, 1), 'DisplayName', 'Supragranular to Infragranular');
% % %     xlim([2 100]);
% % %     hold off;
% % %     title(['Coherence Analysis for Condition ', num2str(conditionNumber), ': ', IDX(conditionNumber).conditionString]);
% % %     xlabel('Frequency (Hz)');
% % %     ylabel('Average Coherence');
% % %     legend;
% % % 
% % % end
% % % 



