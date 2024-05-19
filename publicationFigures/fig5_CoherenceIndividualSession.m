%% Figure 5: Coherence on individual sessions
% ex: load('sortedData_221214_J_bmcBRFS001.mat')
close all

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
        for channel1 = 1:32
            for channel2 = 1:channel1
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
                coherenceMatrix1(channel1, channel2, count) = median(coherence(2:15));  % You can use median or any other aggregation method
            end
        end
    end
end


% Average across trials and save output
averageCoherenceMatrix = median(coherenceMatrix1,3);

% Visualize the coherence matrix
f = figure;
% set(f,"Position",[-1882 -65 1860 894])
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
        for channel1 = 1:32
            for channel2 = 1:channel1
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
                coherenceMatrix2(channel1, channel2, count) = median(coherence(2:15));  % You can use median or any other aggregation method
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
sgtitle({'221202_J_brfs002','Probe #1'},"Interpreter","none")
cd('S:\formattedDataOutputs')
saveas(f,'coherence_221202_J_brfs002_probe1.png')

