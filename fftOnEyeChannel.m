
%% Setup
clear
close all
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';



%% Create dir loop
cd(outputDir)
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

%% Calculate FFT
% In her original email, Ayelet wrote:
% do you ever have a condition in which one stimulus to one eye onsets ...
% earlier, before another stimulus to the other eye is presented ...
% (so at the end both are on)?

% Conditions of interest: 12 vs 10
% 12: 'BRFS IC Adapted Flash. NPO LeftEye adapting - PO RightEye flashed'
% 10: 'BRFS IC Adapted Flash. PO LeftEye adapting - NPO RightEye flashed'

% Preferred MUA response is PO_RE for ch [3:17]

prefFlash = IDX(12).MUAe;
nullEye = IDX(10).MUAe;

for i = 1:15
    varForFFT_pref(:,i) = prefFlash{i,2}(:,3);
    varForFFT_null(:,i) = nullEye{i,2}(:,3);
end
testFFT = fft(varForFFT_pref);
avgFFT = mean(testFFT,2);



%% Continuous line plots for individual units

% Figure settings
set(0,'DefaultFigureWindowStyle','normal')
f = figure;
set(f,'Position',[1000 74.3333 1.1197e+03 1.1633e+03])
plot(abs(avgFFT))
xlim([0 50])
xlabel('Freq (Hz)')
ylim([0 100])
ylabel('|FFT(X)|')
title('FFT on BRFS preferred stimulus flash')