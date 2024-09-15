
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
% Desired output is (1001ms X 32ch X __trls)
timeWindow = 1; % first 800 ms
ms = 1:500;
ch = 10:13;
for trl = size(IDX(5).MUAe,1)
    monocular.Ori1_RightEye(ms,ch,trl) = IDX(5).MUAe{trl,timeWindow}(ms,ch);
end
for trl = size(IDX(6).MUAe,1)
    monocular.Ori2_LeftEye(ms,ch,trl) = IDX(6).MUAe{trl,timeWindow}(ms,ch);
end
for trl = size(IDX(7).MUAe,1)
    monocular.Ori2_RightEye(ms,ch,trl) = IDX(7).MUAe{trl,timeWindow}(ms,ch);
end
for trl = size(IDX(8).MUAe,1)
    monocular.Ori1_LeftEye(ms,ch,trl) = IDX(8).MUAe{trl,timeWindow}(ms,ch);
end

monocularAvg.Ori1_RightEye = mean(mean(monocular.Ori1_RightEye,3),2) ;   % end result should be 1001 x ch
monocularAvg.Ori2_LeftEye = mean(mean(monocular.Ori2_LeftEye,3),2) ;   % end result should be 1001 x ch
monocularAvg.Ori2_RightEye = mean(mean(monocular.Ori2_RightEye,3),2) ;   % end result should be 1001 x ch
monocularAvg.Ori1_LeftEye = mean(mean(monocular.Ori1_LeftEye,3),2) ;   % end result should be 1001 x ch

figure
plot(ms,monocularAvg.Ori1_RightEye); hold on
plot(ms,monocularAvg.Ori2_LeftEye); hold on
plot(ms,monocularAvg.Ori2_RightEye); hold on
plot(ms,monocularAvg.Ori1_LeftEye); hold on

hmm.... This is pretty noisy. It does not look like the MUA tells a clear preference averaged over the units

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