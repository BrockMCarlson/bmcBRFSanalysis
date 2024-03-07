% Example chatGPT mscohere code
%% Coherence matrix across channels
% Choose the condition
% close all
titleString = {'Binoc PS','Binoc NS', 'Monoc PS'};
count = 0;
for conditionNumber = [1 2 11]
    % conditionNumber = 11;
    
    % Get the number of trials for the chosen condition
    numTrials = size(IDX(conditionNumber).LFP_bb, 1);
    
    % Parameters for mscohere
    fs = 1000;  % Sampling frequency in Hz
    windowSize = 256;  % Window size for computing the coherence
    overlap = 10;  % Overlap between windows
    
    % Initialize coherence matrix
    coherenceMatrix = zeros(32, 32, numTrials);
    % Loop through all trials and compute coherence for each channel pair
    for trialNumber = 1:numTrials
        for channel1 = 3:30
            for channel2 = 3:30
                % Extract the LFP_bb data for the chosen channels and trial
                lfpGammaData1 = IDX(conditionNumber).LFP_bb{trialNumber, 1}(350:1001, channel1);
                lfpGammaData2 = IDX(conditionNumber).LFP_bb{trialNumber, 1}(350:1001, channel2);
    
                % Compute coherence
                [coherence, freq] = mscohere(lfpGammaData1, lfpGammaData2, windowSize, overlap, [], fs);
    
                % Store coherence in the matrix
                coherenceMatrix(channel1, channel2, trialNumber) = mean(coherence);  % You can use mean or any other aggregation method
            end
        end
    end
    
    % Average across trials (optional)
    averageCoherenceMatrix = mean(coherenceMatrix(3:30,3:30,:),3);
    
    % Visualize the coherence matrix
    figure;
    imagesc(averageCoherenceMatrix);
    colormap('jet');
    clim([.5 .65])
    colorbar;
    title('Coherence Matrix across Channels');
    xlabel('Channel');
    ylabel('Channel');
    count = count+1;
    title(titleString{count});
end


%% Compare conditions

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

% Time window of interest
tm = 350:1000;

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
        lfpSupragranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, supragranularChannels);
        lfpGranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, granularChannels);
        lfpInfragranular = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, infragranularChannels);

        % Compute coherence for each layer pair
        [coherenceGranularSupragranularTrial, freq] = mscohere(lfpGranular, lfpSupragranular, windowSize, overlap, nfft, fs);
        [coherenceGranularInfragranularTrial, ~] = mscohere(lfpGranular, lfpInfragranular, windowSize, overlap, nfft, fs);
        [coherenceSupragranularInfragranularTrial, ~] = mscohere(lfpSupragranular, lfpInfragranular, windowSize, overlap, nfft, fs);
        
        % Extract frequency values within the specified range
        freqIndices = find(freq >= freqRange(1) & freq <= freqRange(2));
        freqLimValues = freq(freq >= freqRange(1) & freq <= freqRange(2));
        
        % Store mean coherence data for each layer pair
        coherenceGranularSupragranular(trialNumber, :) = mean(coherenceGranularSupragranularTrial, 2);
        coherenceGranularInfragranular(trialNumber, :) = mean(coherenceGranularInfragranularTrial, 2);
        coherenceSupragranularInfragranular(trialNumber, :) = mean(coherenceSupragranularInfragranularTrial, 2);

    end

    % Plot the average coherence across all trials for each layer pair on a log-scaled x-axis
    figure;
    semilogx(freq, mean(coherenceGranularSupragranular, 1), 'DisplayName', 'Granular to Supragranular');
    hold on;
    semilogx(freq, mean(coherenceGranularInfragranular, 1), 'DisplayName', 'Granular to Infragranular');
    semilogx(freq, mean(coherenceSupragranularInfragranular, 1), 'DisplayName', 'Supragranular to Infragranular');
    xlim([2 100]);
    hold off;
    title(['Coherence Analysis for Condition ', num2str(conditionNumber), ': ', IDX(conditionNumber).conditionString]);
    xlabel('Frequency (Hz)');
    ylabel('Average Coherence');
    legend;

end




