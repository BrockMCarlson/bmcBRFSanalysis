% Examining ringing in the grand average signal.

codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis';
cd(codeDir)

dataDir = 'D:\sortedData_240229';
cd(dataDir)

%% For loop
cd(dataDir)
allDataFiles = dir('**/*sortedData*.mat');
fullDataset = nan(1001,length(allDataFiles));
for file = 1:length(allDataFiles)
    % load data
    cd(dataDir)
    fileToLoad = allDataFiles(file).name;
    load(fileToLoad)
    sessionLabel = allDataFiles(file).name(12:end-4);

    % Trial average across trials for a single unit
    ch = 15;
    singleChDat = nan(1001,length(IDX(1).LFP_beta1));
    for i = 1:length(IDX(1).LFP_beta1)
        singleChDat(:,i) = IDX(1).LFP_beta1{i,1}(:,ch);
    end
    singleChAvg_mean    = mean(singleChDat,2);
    singleChAvg_median  = median(singleChDat,2);

    blAvg_mean = mean(singleChAvg_mean(1:200,1),1);
    singleChBlSub_mean = singleChAvg_mean - blAvg_mean;
    fullDataset_mean(:,file) = singleChBlSub_mean;

    blAvg_median = median(singleChAvg_median(1:200,1),1);
    singleChBlSub_median = singleChAvg_median - blAvg_median;
    fullDataset_median(:,file) = singleChBlSub_median;
    disp('finished file number...')
    disp(file)
    clear IDX

end

%% IDX plot
% % plot(-200:800,singleChBlSub)
% % xlabel('Time from dioptic onset(ms)')
% % ylabel('Band-limited-Power (uV)')
% % vline(0)
% % title('Single electrode example, average over all dioptic onsets')

figure
plot(-200:800,fullDataset_median)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Ch 15 - each line is trial-averaged data from 1 penetration')
    

fullDatMean = mean(fullDataset_mean,2);
figure
plot(-200:800,fullDatMean)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Grand average (mean) of Beta1 LFP recorded at ch15 - average across sessions')

fullDatMedian = median(fullDataset_median,2);
figure
plot(-200:800,fullDatMedian)
xlabel('Time from dioptic onset(ms)')
ylabel('Band-limited-Power (uV)')
vline(0)
title('Grand average (median) of Beta1 LFP recorded at ch15 - average across sessions')