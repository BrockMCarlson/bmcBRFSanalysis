
%% Setup
clear
close all
codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
sortedDataDir = 'D:\bmcBRFS_sortedData';



%% Create dir loop
cd(sortedDataDir)
sortedDataFiles = dir('*sortedData*');
close all
dataFile = 1;
fileToLoad = sortedDataFiles(dataFile).name;
load(fileToLoad)

% Variables
probeLength = size(IDX(1).LFP_bb{1,1},2);
sdftm = STIM(1).sdftm;
timeLength = length(sdftm);
baselineTimeIndex = find(sdftm == -100):find(sdftm == 0); 
xAxisTime = sdftm;

%% Demonstrate monocular preference
% BRFS-like Congruent Adapted Flash is my condition of interest
% IDX [5:8] 
timeWindow = 1; % first 800 ms
ms = 1:1000;
ch = 3:17;
for cond = 5:8
    for trl = size(IDX(cond).MUAe,1)
        monocular(cond,trl,ms,ch) = IDX(cond).MUAe{trl,timeWindow}(ms,ch);
    end
end

You are here: figure out how to simply this system for not and in the future.
4-D variables are untennable. size of monocular returns [1,9] ??

monocularAvg = mean(thingToAverage,"all") ;   % end result should be 1001 x 32

%% Plot response timecourse
% Preferred MUA response is PO_RE for ch [3:17]

prefFlash = IDX(12).MUAe;
nullEye = IDX(10).MUAe;

for i = 1:15
    varForFFT_pref(:,i) = prefFlash{i,2}(:,3:17);
    varForFFT_null(:,i) = nullEye{i,2}(:,3:17);
end


%% Calculate FFT
testFFT = fft(varForFFT_null);
avgFFT = mean(testFFT,2);



%% plot FFT

% Figure settings
set(0,'DefaultFigureWindowStyle','docked')
f = figure;
set(f,'Position',[1000 74.3333 1.1197e+03 1.1633e+03])
plot(abs(avgFFT))
xlim([0 50])
xlabel('Freq (Hz)')
ylim([0 100])
ylabel('|FFT(X)|')
title('FFT on BRFS preferred stimulus flash')