%% bmcBRFS_coherenceMatrix
% Are there observable differences between trial-types with LFP coherence?
% initialize variables
fileToLoad = 'S:\bmcBRFS_sortedData_Nov23\sortedData_221206_J_bmcBRFS002.mat';
load(fileToLoad)


%% Coherence matrix across channels
% Choose the condition
close all
averageCoherenceMatrix = nan(32,32,4);
count = 0;
for conditionNumber = [9 10 11 12]
    % conditionNumber = 11;
    
    % Get the number of trials for the chosen condition
    numTrials = size(IDX(conditionNumber).LFP_bb, 1);
    
    % Parameters for mscohere
    fs = 1000;  % Sampling frequency in Hz
    windowSize = 256;  % Window size for computing the coherence
    overlap = windowSize/2;  % Overlap between windows
    tm = 489:1001; % Time window of data. Last 512ms of trial. 
    
    % Initialize coherence matrix
    coherenceMatrix = zeros(32, 32, numTrials);
    % Loop through all trials and compute coherence for each channel pair
    for trialNumber = 1:numTrials
        for channel1 = 1:32
            for channel2 = 1:32
                % Extract the LFP_bb data for the chosen channels and trial
                lfpGammaData1 = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, channel1);
                lfpGammaData2 = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, channel2);
    
                % Compute coherence
                [coherence, freq] = mscohere(lfpGammaData1, lfpGammaData2, windowSize, overlap, [], fs);
    
                % Store coherence in the matrix
                % warning('try maximum coherence values across frequency')
                coherenceMatrix(channel1, channel2, trialNumber) = max(coherence);  % You can use median or any other aggregation method
            end
        end
    end
    
    % Average across trials and save output
    count = count + 1;
    averageCoherenceMatrix(:,:,count) = median(coherenceMatrix,3);
    
end


% Visualize the coherence matrix
count = 0;
for conditionNumber = [9 10 11 12]
    count = count + 1;
    figure;
    imagesc(averageCoherenceMatrix(:,:,count));
    colormap('jet');
    % clim([.5 .75])
    colorbar;
    xlabel('Channel');
    ylabel('Channel');
    title(strcat(string(conditionNumber),'--',IDX(conditionNumber).conditionString));
end

% Difference plot of coherence matrix
% abs(PS - NS)
% 9 NPO --> PO; 10 PO-->NPO
% abs(cond9-cond10)
% abs(cond12-cond11)
differenceMatrix(:,:,1) = abs(averageCoherenceMatrix(:,:,1) - averageCoherenceMatrix(:,:,2));
differenceMatrix(:,:,2) = abs(averageCoherenceMatrix(:,:,4) - averageCoherenceMatrix(:,:,3));

% Visualize the diff  matrix
figure;
imagesc(differenceMatrix(:,:,1));
colormap('bone');
% clim([.5 .75])
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Diff 1');

% Visualize the diff  matrix
figure;
imagesc(differenceMatrix(:,:,2));
colormap('bone');
% clim([.5 .75])
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Diff 2');


%% Compare laminar compartments

clearvars -except IDX
% Default parameter values
windowSize = 256;  % Window size for computing the coherence
overlap = windowSize/2;  % Overlap between windows
nfft = 512;  % Number of points to use in the FFT
fs = 1000;  % Sampling frequency in Hz

% Define channel ranges for different layers
supragranularChannels = 3:7;
granularChannels = 8:12;
infragranularChannels = 13:17;

% Frequency range of interest
freqRange = [1 100];

% Loop through condition numbers
for conditionNumber = [12, 11]
    % Get the number of trials for the chosen condition
    numTrials = size(IDX(conditionNumber).LFP_bb, 1);

    % Initialize freqLimValues outside the loop
    [~, freqLimValues] = mscohere(zeros(windowSize, 1), zeros(windowSize, 1), windowSize, overlap, nfft, fs);
    
    % Initialize arrays to store coherence data for each layer pair
    coherenceGranularSupragranular = zeros(numTrials, length(freqLimValues));
    coherenceGranularInfragranular = zeros(numTrials, length(freqLimValues));
    coherenceSupragranularInfragranular = zeros(numTrials, length(freqLimValues));

    % Loop through all trials and compute coherence for each layer pair
    for trialNumber = 1:numTrials
        % Extract the LFP_bb data for each layer
        lfpSupragranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, supragranularChannels);
        lfpGranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, granularChannels);
        lfpInfragranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, infragranularChannels);

        % Compute coherence for each layer pair
        [coherenceGranularSupragranularTrial, freq] = mscohere(lfpGranular, lfpSupragranular, windowSize, overlap, nfft, fs);
        [coherenceGranularInfragranularTrial, ~] = mscohere(lfpGranular, lfpInfragranular, windowSize, overlap, nfft, fs);
        [coherenceSupragranularInfragranularTrial, ~] = mscohere(lfpSupragranular, lfpInfragranular, windowSize, overlap, nfft, fs);
        
        % Extract frequency values within the specified range
        freqIndices = find(freq >= freqRange(1) & freq <= freqRange(2));
        freqLimValues = freq(freq >= freqRange(1) & freq <= freqRange(2));
        
        % Store median coherence data for each layer pair
        coherenceGranularSupragranular(trialNumber, :) = median(coherenceGranularSupragranularTrial, 2);
        coherenceGranularInfragranular(trialNumber, :) = median(coherenceGranularInfragranularTrial, 2);
        coherenceSupragranularInfragranular(trialNumber, :) = median(coherenceSupragranularInfragranularTrial, 2);

    end

    % Plot the average coherence across all trials for each layer pair on a log-scaled x-axis
    figure;
    semilogx(freq, median(coherenceGranularSupragranular, 1), 'DisplayName', 'Granular to Supragranular');
    hold on;
    semilogx(freq, median(coherenceGranularInfragranular, 1), 'DisplayName', 'Granular to Infragranular');
    semilogx(freq, median(coherenceSupragranularInfragranular, 1), 'DisplayName', 'Supragranular to Infragranular');
    xlim([2 100]);
    hold off;
    title(['Coherence Analysis for Condition ', num2str(conditionNumber), ': ', IDX(conditionNumber).conditionString]);
    xlabel('Frequency (Hz)');
    ylabel('Average Coherence');
    legend;

end




