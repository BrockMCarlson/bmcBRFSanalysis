%% bmcBRFS_coherenceMatrix
% Are there observable differences between trial-types with LFP coherence?
% initialize variables
% % clear
% % close all


%% For loop
codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis';
cd(codeDir)
% % dataDir = 'D:\sortedData_240229';
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
% officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

averageCoherenceMatrix = nan(15,15,2,size(officLamAssign,1));


for file = 1:size(officLamAssign,1)
    if isnan(officLamAssign.stFold4c(file))
        continue
    end
    % load data
    cd(dataDir)
    probeName = char(officLamAssign.Session_probe_(file,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');
    load(fileToLoad)


    granBtm = officLamAssign.stFold4c(file); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',fileToLoad))
        continue
    end
    
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
    if any(v1Ch > 32) || any(v1Ch < 1)
        warning('skipping session without full column for now')
        disp(fileToLoad)
        continue
        % % columnNames = columnNames((v1Ch >= 1) & (v1Ch<=32));
        % % v1Ch        = v1Ch((v1Ch >= 1) & (v1Ch<=32));
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
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(IDX(cond).LFP_bb{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(IDX(cond).LFP_bb{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(IDX(cond).LFP_bb{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(IDX(cond).LFP_bb{trl,1}(:,v1Ch)); 
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
    
    % Skip this file if no tuned units are found
    % % DATAOUT(file).numberOFUnits = sum(tuned);
    % % if sum(tuned) == 0
    % %     warning(strcat('No tuned channels on',fileToLoad))
    % %     continue
    % % end
    % % 
    % % % Creat channel array to use in subsequent steps - only plot tuned units.
    % % chTuned = ch(logical(tuned));
    % % 
    % % 


    % Create array of preference-based data
    % concatenate two timecourses into array
    tm_full = -200:1600; % 1801 total timepoints
    tm1 = 1:801;
    tm2 = 1:1001;
    tm2_concat = 801:1801;

    % pre allocate
    maxTrlLength = max([length(IDX(9).correctTrialIndex),...
        length(IDX(10).correctTrialIndex),...
        length(IDX(11).correctTrialIndex),...
        length(IDX(12).correctTrialIndex)]);
    array_dichopticAdapted_pref = nan(1801,length(v1Ch),maxTrlLength);
    array_dichopticAdapted_null = nan(1801,length(v1Ch),maxTrlLength);
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
    overallNull = mode(nullCondOnFlash);
    %% Coherence matrix across channels
    % Choose the condition
    close all
    count = 0;
    for conditionNumber = [overallPref overallNull]
        
        % Get the number of trials for the chosen condition
        numTrials = size(IDX(conditionNumber).LFP_bb, 1);
        
        % Parameters for mscohere
        fs = 1000;  % Sampling frequency in Hz
        windowSize = 256;  % Window size for computing the coherence
        overlap = windowSize/2;  % Overlap between windows
        tm_coher = 490:1001; % Time window of data. Last 512ms of trial. 
        tm = 1:1001; % Time window of data. Last 512ms of trial. 
        
        % Initialize coherence matrix
        coherenceMatrix = nan(15,15, numTrials);
        % Loop through all trials and compute coherence for each channel pair
        for trialNumber = 1:numTrials
            for channel1 = v1Ch
                for channel2 = v1Ch(1):channel1
                    % Extract the LFP_bb data for the chosen channels and trial
                    lfpGammaData1 = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, channel1);
                    lfpGammaData2 = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, channel2);
                    %baseline subtract
                    bl1 = median(lfpGammaData1(1:200));
                    lfp_blsub_1 = lfpGammaData1 - bl1;
                    bl2 = median(lfpGammaData2(1:200));
                    lfp_blsub_2 = lfpGammaData2 - bl2;
    
                    % Compute coherence
                    [coherence, freq] = mscohere(lfp_blsub_1(tm_coher), lfp_blsub_2(tm_coher), windowSize, overlap, [], fs);
        
                    % Store coherence in the matrix
                    coherenceMatrix(channel1-v1Ch(1)+1, channel2-v1Ch(1)+1, trialNumber) = median(coherence(2:15));  % You can use median or any other aggregation method
                end
            end
        end
        
        % Average across trials and save output
        count = count + 1;
        averageCoherenceMatrix(:,:,count,file) = median(coherenceMatrix,3);
        
    end
    dat_LFP{file} = IDX.LFP_bb;
    disp(strcat('Done with file number: ',string(file)))
end

grandAverageCoherence = median(averageCoherenceMatrix,4,"omitmissing");

disp('analyzing freq')
disp(freq(2:15))


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
differenceMatrix(:,:) = abs(grandAverageCoherence(:,:,1) - grandAverageCoherence(:,:,2));

% Visualize the diff  matrix
figure;
imagesc(differenceMatrix(:,:));
colormap('bone');
% clim([.5 .75])
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Diff of Pref vs Null');





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



