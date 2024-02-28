% Examining ringing in the grand average signal.

codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis';
cd(codeDir)

dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(dataDir)
load('sortedData_221103_J_bmcBRFS001.mat')

%% 1. IDX plots
close all
% plot single trial
figure
subplot(2,1,1)
holder = IDX(1).LFP_beta1{5,1};
plot(-200:800,holder)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Single trial example #1, all 32 electrodes')
subplot(2,1,2)
holder = IDX(1).LFP_beta1{6,1};
plot(-200:800,holder)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Single trial example #2, all 32 electrodes')


% Trial average for a single unit
figure
ch = 10;
singleChDat = nan(1001,length(IDX(1).LFP_beta1));
for i = 1:length(IDX(1).LFP_beta1)
    singleChDat(:,i) = IDX(1).LFP_beta1{i,1}(:,ch);
end
singleChAvg = mean(singleChDat,2);
plot(-200:800,singleChAvg)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Single electrode example, average over all dioptic onsets')


% Grand average for all units on a single session over all trials
figure
allChDat = nan(1001,32,length(IDX(1).LFP_beta1));
for i = 1:length(IDX(1).LFP_beta1)
    allChDat(:,:,i) = IDX(1).LFP_beta1{i,1};
end
allChAvg = mean(allChDat,3);
plot(-200:800,allChAvg)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Single session example of all electrodes, averaged over all dioptic onsets')


%% filteredDataPlots
cd(dataDir)
load('filteredData_221103_J_bmcBRFS001.mat')

filteredData = LFP_alphaBeta(1100001:1108000,:);
plot(filteredData)
xlabel('Time from arbritrary start')
ylabel('Band-limited-Power (uV)')
title('8 seconds of 10-30Hz')