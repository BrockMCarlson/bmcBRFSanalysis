%% bmcBRFS_coherenceMatrix
% Are there observable differences between trial-types with LFP coherence?
% initialize variables
clear
close all


%% For loop
dataDir = 'D:\sortedData_240229';
% dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
allDataFiles = dir('**/*sortedData*.mat');
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\coherenceFigs';
% % averageCoherenceMatrix = nan(15,15,2,length(allDataFiles));


for file = 1:size(officLamAssign,1)
    
    % load data
    cd(dataDir)
    probeName = char(officLamAssign.SessionProbe(file,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');
    load(fileToLoad)
    sessionLabelCell = table2cell(officLamAssign(file,1));
    sessionLabel = sessionLabelCell{1};
    chStr = officLamAssign.ChtoUse(file,1);
    if strcmp(chStr,"1:32")
        ch = 1:32;
    elseif strcmp(chStr,"33:64")
        ch = 33:64;
    end
    %% Parameters for mscohere
    fs = 1000;  % Sampling frequency in Hz
    windowSize = 256;  % Window size for computing the coherence
    overlap = windowSize/2;  % Overlap between windows
    % tm_coher = 1:512; % Time window of data. Last 512ms of trial. 
    tm_coher = 490:1001; % Time window of data. Last 512ms of trial. 
    tm = 1:1001; % Time window of data. Last 512ms of trial. 
    
    %% Dioptic        
    % Get the number of trials for the chosen condition
    totalTrials = size(IDX(1).LFP_bb, 1) + size(IDX(2).LFP_bb, 1);
    
    
    % Initialize coherence matrix
    coherenceMatrix1 = nan(32,32, totalTrials);
    % Loop through all trials and compute coherence for each channel pair
    count = 0;
    for conditionNumber = 1:2
        for trialNumber = 1:size(IDX(conditionNumber).LFP_bb, 1)
            count = count + 1;
            for channel1 = ch
                for channel2 = ch(1):channel1
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
                    if strcmp(chStr,"1:32")
                        coherenceMatrix1(channel1, channel2, count) = median(coherence(2:15));  % You can use median or any other aggregation method
                    elseif strcmp(chStr,"33:64")
                        coherenceMatrix1(channel1-32, channel2-32, count) = median(coherence(2:15));  % You can use median or any other aggregation method
                    end
                end
            end
        end
    end
    
    
    % Average across trials and save output
    averageCoherenceMatrix = median(coherenceMatrix1,3);
    
    % Visualize the coherence matrix
    f = figure;
    set(f,"Position",[7.6667 55.6667 2.5373e+03 1.2953e+03])
    ax(1) = subplot(1,4,1);
    imagesc(averageCoherenceMatrix);
    colormap(ax(1),'jet');
    e = colorbar;
    e.Label.String = "Inter-contact coherence";
    e.Label.Rotation = 270;
    xlabel('Channel');
    ylabel('Channel');
    title('Coherence, Dioptic');
    
    
    %% Dichoptic
    % Get the number of trials for the chosen condition
    totalTrials = size(IDX(3).LFP_bb, 1) + size(IDX(4).LFP_bb, 1);
    
    
    
    % Initialize coherence matrix
    coherenceMatrix2 = nan(32,32, totalTrials);
    % Loop through all trials and compute coherence for each channel pair
    count = 0;
    for conditionNumber = 3:4
        for trialNumber = 1:size(IDX(conditionNumber).LFP_bb, 1)
            count = count + 1;
            for channel1 = ch
                for channel2 = ch(1):channel1
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
                    if strcmp(chStr,"1:32")
                        coherenceMatrix2(channel1, channel2, count) = median(coherence(2:15));  % You can use median or any other aggregation method
                    elseif strcmp(chStr,"33:64")
                        coherenceMatrix2(channel1-32, channel2-32, count) = median(coherence(2:15));  % You can use median or any other aggregation method
                    end
                end
            end
        end
    end
    
    
    
    
    % Average across trials and save output
    averageCoherenceMatrix = median(coherenceMatrix2,3);
    
    % Visualize the coherence matrix
    ax(2) = subplot(1,4,2);
    imagesc(averageCoherenceMatrix);
    colormap(ax(2),'jet');
    e = colorbar;
    e.Label.String = "Inter-contact coherence";
    e.Label.Rotation = 270;
    xlabel('Channel');
    ylabel('Channel');
    title('Coherence, Dichoptic');
    
    
    %% Difference plot
    for ch1 = 1:32
        for ch2 = 1:ch1
            [h(ch1,ch2),p(ch1,ch2),ci(ch1,ch2,:),stats(ch1,ch2)]...
                = ttest2(...
                coherenceMatrix1(ch1,ch2,:),...
                coherenceMatrix2(ch1,ch2,:));
            tStat(ch1,ch2) = stats(ch1,ch2).tstat;
        end
    end
    
    % Visualize the coherence matrix
    ax(3) = subplot(1,4,3);
    imagesc(tStat);
    colormap(ax(3),'bone');
    e = colorbar;
    e.Label.String = "tStat";
    e.Label.Rotation = 270;
    xlabel('Channel');
    ylabel('Channel');
    title('Dioptic vs Dichoptic tScoreMap');
    
    % Visualize the coherence matrix
    ax(4) = subplot(1,4,4);
    imagesc(h);
    colormap(ax(4),'bone');
    e = colorbar;
    e.Label.String = "h = 1 or 0";
    e.Label.Rotation = 270;
    xlabel('Channel');
    ylabel('Channel');
    title('Dioptic vs Dichoptic Hypothesis test');
    sgtitle(probeName(1:end-1),"Interpreter","none")
    cd(plotDir)
    saveName = strcat('coherence_',probeName(1:end-1),'.png');
    saveas(f,saveName) 
    close all
      
end

% % grandAverageCoherence = median(averageCoherenceMatrix,4,"omitmissing");
% % 
% % disp('analyzing freq')
% % disp(freq(2:15))
% % 
% % 
% % % Visualize the coherence matrix
% % figure;
% % imagesc(grandAverageCoherence(:,:,1));
% % colormap('jet');
% % colorbar;
% % xlabel('Channel');
% % ylabel('Channel');
% % title('Preferred Stimulus BRFS flash Grand Average');
% % 
% % figure;
% % imagesc(grandAverageCoherence(:,:,2));
% % colormap('jet');
% % colorbar;
% % xlabel('Channel');
% % ylabel('Channel');
% % title('Non-preferred stimulus BRFS flash Grand Average');
% % 
% % % Difference plot of coherence matrix
% % % abs(PS - NS)
% % % 9 NPO --> PO; 10 PO-->NPO
% % % abs(cond9-cond10)
% % % abs(cond12-cond11)
% % differenceMatrix = nan(15,15);
% % differenceMatrix(:,:) = abs(grandAverageCoherence(:,:,1) - grandAverageCoherence(:,:,2));
% % 
% % % Visualize the diff  matrix
% % figure;
% % imagesc(differenceMatrix(:,:));
% % colormap('bone');
% % % clim([.5 .75])
% % colorbar;
% % xlabel('Channel');
% % ylabel('Channel');
% % title('Diff of Pref vs Null');





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



