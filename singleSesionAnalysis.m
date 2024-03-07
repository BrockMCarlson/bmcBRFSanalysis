% load('S:\bmcBRFS_sortedData_Nov23\sortedData_221123_J_bmcBRFS001.mat')

% Ideal dimension for a variable are ch x tm x trial
LEcond = 9;
LEdat = nan(32,1201,length(IDX(LEcond).LFP_alpha));
for ch = 1:32
    for trl = 1:length(IDX(LEcond).LFP_alpha)
        LEdat(ch,:,trl) = IDX(LEcond).LFP_alpha{trl,2}(:,ch);
    end
end
LEavgBRFS = mean(LEdat,3);
blAvg = mean(LEavgBRFS(:,1:200),2);
LEblSubBRFS = LEavgBRFS - blAvg;
figure
plot(-200:1000,LEblSubBRFS)
title(IDX(LEcond).conditionString)
% ylim([-200 25])
vline(0)
vline(250)

REcond = 10;
REdat = nan(32,1201,length(IDX(REcond).LFP_alpha));
for ch = 1:32
    for trl = 1:length(IDX(REcond).LFP_alpha)
        REdat(ch,:,trl) = IDX(REcond).LFP_alpha{trl,2}(:,ch);
    end
end
REavgBRFS = mean(REdat,3);
blAvg = mean(REavgBRFS(:,1:200),2);
REblSubBRFS = REavgBRFS - blAvg;
figure
plot(-200:1000,REblSubBRFS)
title(IDX(REcond).conditionString)
% ylim([-200 25])
vline(0)
vline(250)


% Plot CSD
x = [-200:1000];
y = 1:32;
% trial-types
LEcond = 1;
LEdat = nan(32,1201,length(IDX(LEcond).CSD_bb));
for ch = 1:32
    for trl = 1:length(IDX(LEcond).CSD_bb)
        LEdat(ch,:,trl) = IDX(LEcond).CSD_bb{trl,1}(:,ch);
    end
end
LEavgBRFS = mean(LEdat,3);
blAvg = mean(LEavgBRFS(:,1:200),2);
LEblSubBRFS = LEavgBRFS - blAvg;
figure
imagesc(x,y,LEblSubBRFS)
colormap(flipud(jet))
c2 = colorbar;
title(IDX(LEcond).conditionString)
vline(0)

REcond = 2;
REdat = nan(32,1201,length(IDX(REcond).CSD_bb));
for ch = 1:32
    for trl = 1:length(IDX(REcond).CSD_bb)
        REdat(ch,:,trl) = IDX(REcond).CSD_bb{trl,1}(:,ch);
    end
end
REavgBRFS = mean(REdat,3);
blAvg = mean(REavgBRFS(:,1:200),2);
REblSubBRFS = REavgBRFS - blAvg;
figure
imagesc(x,y,REblSubBRFS)
colormap(flipud(jet))
c1 = colorbar;
title(IDX(REcond).conditionString)
vline(0)



%% Plot all condition averages in LFP
for cond = 1:20
    trialNumber = size(IDX(cond).LFP_bb,1);
    fullLineToPlot = nan(1801,32,trialNumber);
    for trl = 1:size(IDX(cond).LFP_bb,1)
        fullLineToPlot(1:1001,:,trl) = IDX(cond).LFP_bb{trl,1};
        fullLineToPlot(901:1801,:,trl) = IDX(cond).LFP_bb{trl,2}(101:1001,:);
    end
    avgDat = mean(fullLineToPlot,3);
    xVals = -200:1600;
    figure
    plot(xVals,avgDat)
    vline(0)
    vline(800)
    ylabel('MicroVolts (uV)')
    title(IDX(cond).conditionString)
    pause(2)
end

trl = 3;
cond = 1;
xVals = -200:1600;
fullLineToPlot = nan(1801,32);
fullLineToPlot(1:1001,:) = IDX(cond).LFP_bb{trl,1}(:,1:32);
fullLineToPlot(901:1801,:) = IDX(cond).LFP_bb{trl,2}(101:1001,1:32);
figure
plot(xVals,fullLineToPlot)
vline(0)
% vline(800)
ylabel('MicroVolts (uV)')
title(IDX(cond).conditionString)

trl = 2;
cond = 11;
xVals = -200:1600;
fullLineToPlot = nan(1801,32);
fullLineToPlot(1:1001,:) = IDX(cond).LFP_bb{trl,1}(:,1:32);
fullLineToPlot(901:1801,:) = IDX(cond).LFP_bb{trl,2}(101:1001,1:32);
figure
plot(xVals,fullLineToPlot)
vline(0)
vline(800)
ylabel('MicroVolts (uV)')
title(IDX(cond).conditionString)





%% LFP_bb is the beautifully time locked across all electrodes
figure
plot(IDX(1).LFP_bb{1,1})
%MUAe is a mess in contrast
figure
plot(IDX(1).MUAe{1,1})
%of course, this is across all electrodes, which - if you think about it -
%is a weird way to look at the data. Normally we look at a single electrode
%averaged across trials.

% So what this actually shows is something quite shocking. That low
% frequency fluctuations are the largest and most impactful on the overall
% signal, and may be totaly lost in spiking data. 
% (of course... it would be helpful if you actually looked at this KLS
% data.........)
test = IDX(1).LFP_bb(:,1);
timeElTrial = nan(1001,32,length(IDX(1).LFP_bb));
for trl = 1:length(IDX(1).LFP_bb)
    timeElTrial(1:1001,1:32,trl) = test{trl};
end
avgDat = mean(timeElTrial,3);
figure
plot(avgDat)




%% Coherence matrix across channels
% Choose the condition
% close all

for conditionNumber = [9 10 11 12]
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
        for channel1 = 1:32
            for channel2 = 1:32
                % Extract the LFP_bb data for the chosen channels and trial
                lfpGammaData1 = IDX(conditionNumber).MUAe{trialNumber, 2}(:, channel1);
                lfpGammaData2 = IDX(conditionNumber).MUAe{trialNumber, 2}(:, channel2);
    
                % Compute coherence
                [coherence, freq] = mscohere(lfpGammaData1, lfpGammaData2, windowSize, overlap, [], fs);
    
                % Store coherence in the matrix
                coherenceMatrix(channel1, channel2, trialNumber) = median(coherence);  % You can use mean or any other aggregation method
            end
        end
        % Visualize individual trials
        % Single contact
        % % error('remov')
        % % figure
        % % set(gcf,"Position",[-1850 -10 1759 610])
        % % subplot(1,2,1)
        % % plot(-200:1000,IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, 10))
        % % hold on
        % % plot(-200:1000,IDX(conditionNumber).LFP_bb{trialNumber, 2}(:, 11))
        % % vline(0)
        % % title('contact 10 and contact 11')
        % % subplot(1,2,2)
        % % [coherence, freq] = mscohere(...
        % %     IDX(conditionNumber).LFP_bb{trialNumber, 2}(400:1001, 10),...
        % %     IDX(conditionNumber).LFP_bb{trialNumber, 2}(400:1001, 11),...
        % %     windowSize, overlap, [], fs);
        % % plot(freq,coherence)
        % % hline(median(coherence))
        % % title('Coherence between signals, with median as output')
        % % 
        % % %Whole electrode
        % % 
        % % figure
        % % set(gcf,"Position",[-1850 -10 1759 610])
        % % % Data 1
        % % subplot(1,2,1)
        % % plot(-200:1000,IDX(conditionNumber).LFP_bb{trialNumber, 2}(:,:))
        % % vline(0)
        % % % imagesc Coherence
        % % subplot(1,2,2)
        % % imagesc(coherenceMatrix(:,:,trialNumber));
        % % colormap('jet');
        % % % clim([.5 .65])
        % % colorbar;
        % % sgtitle('Single trial example of raw LFP');
        % % xlabel('Channel');
        % % ylabel('Channel');
        % % pause(1)
    end
    
    % Average across trials (optional)
    averageCoherenceMatrix = mean(coherenceMatrix,3);
    
    % Visualize the coherence matrix
    figure;
    imagesc(averageCoherenceMatrix);
    colormap('jet');
    % clim([.5 .65])
    colorbar;
    title('MeanAcrossCondition');
    xlabel('Channel');
    ylabel('Channel');
    title(IDX(conditionNumber).conditionString);
end

% Visualize individual trials

    figure;

    imagesc(averageCoherenceMatrix);
    colormap('jet');
    % clim([.5 .65])
    colorbar;
    title('Coherence Matrix across Channels');
    xlabel('Channel');
    ylabel('Channel');


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




