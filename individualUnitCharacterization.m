%% individualUnitCharacterization.m
% Brock M. Carlson
% July 20th, 2023,
% This is a draft analysis for bmcBRFSanalysis. The questions we are
% interested in at an individual unit level are:
% 1. Visually responseive?
% 2. Feature Selective (Tuned to eye? Ori? Eye And Ori?)
% 3. Supressed by dichoptic stimuli?
% 4. Perceptually modulated? 
% These questions are all based on statistic run across presentation
% (trial) from a binned response window, for transient and sustained.

% Additionally, we want to find the laminar position of each unit.

% For continuous data, we want the following plots
    % 1. 2x2 eye v ori (potentially on the same ordinate) 
    % 2. dichoptic vs dichoptic
    % 3. Same stimulus, different history, Preferred vs null
    % 4. Same stimulus, different history, mixed vs mixed.

%% Setup
clear
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';

cd(outputDir)
load('sortedData_211008_B_bmcBRFS001.mat')

% Variables
cond = 1; % Using most excitatory stimulus, condition 1, 'Simult. Dioptic. PO'
probeLength = size(IDX(cond).LFP_bb{1,1},2);
sdftm = STIM(1).sdftm;
timeLength = length(sdftm);
trlLength = size(IDX(cond).LFP_bb,1);
baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
xAxisTime = sdftm;
xAxisLim = [-100 400];
yAxisChannels = 1:32;

% Figure settings
set(0,'DefaultFigureWindowStyle','docked')
f = figure;

%% Continuous line plots for individual units
% Using unit on contact 15 as example plot
chIdx = 15;

% Feature Selectivity - Eye vs ori
% First 800 ms of condition 5-8,9-12, 13-16, 17-20
clear cond
cond.monoc_PO_LE = [8 10 16 18];
cond.monoc_PO_RE = [5 11 13 19];
cond.monoc_NOP_LE = [6 12 14 20];
cond.monoc_NPO_RE = [7 9 15 17];

fields = fieldnames(cond);
for i = 1:length(fields)
    conditions = cond.(fields{i});
    counter = 0;
    for j = 1:length(conditions)
        trlLength = size(IDX(conditions(j)).MUAe,1);
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls.(fields{i})(:,counter) = IDX(conditions(j)).MUAe{trl,1}(:,chIdx);
        end
    end
    muaConditionmedian = median(muaAllTrls.(fields{i}),2);
    mua_baselinemedian = median(muaConditionmedian(baselineTimeIndex,1));
    mua_blSubAvg.(fields{i}) = muaConditionmedian - mua_baselinemedian;
    maxVals(i) = max(mua_blSubAvg.(fields{i}));
end
[maxVal,maxValLoc] = max(maxVals);
titleText = {'Preferred Response';fields{maxValLoc}};

figure
for i = 1: length(fields)
    subplot(2,2,i)
    plot(sdftm,mua_blSubAvg.(fields{i}))
    ylim([-0.5 maxVal])
    xlim([sdftm(1) sdftm(end)])
    vline(0)
    title(fields{i},'interpreter','none')
end
sgtitle(titleText,'interpreter','none')


% Dioptic vs dichoptic
diopticCond = 1;
trlLength = size(IDX(diopticCond).MUAe,1);
counter = 0;
for trl = 1:trlLength
    counter = counter + 1;
    muaAllTrls_dioptic(:,counter) = IDX(diopticCond).MUAe{trl,1}(:,chIdx);
end
muaConditionmedian_dioptic = median(muaAllTrls_dioptic,2);
mua_baselinemedian_dioptic = median(muaConditionmedian_dioptic(baselineTimeIndex,1));
mua_blSubAvg_dioptic = muaConditionmedian_dioptic - mua_baselinemedian_dioptic;

dichopticCond = 4; % We choose this based on the ocular preferences in previous step
trlLength = size(IDX(dichopticCond).MUAe,1);
counter = 0;
for trl = 1:trlLength
    counter = counter + 1;
    muaAllTrls_dichoptic(:,counter) = IDX(dichopticCond).MUAe{trl,1}(:,chIdx);
end
muaConditionmedian_dichoptic = median(muaAllTrls_dichoptic,2);
mua_baselinemedian_dichoptic = median(muaConditionmedian_dichoptic(baselineTimeIndex,1));
mua_blSubAvg_dichoptic = muaConditionmedian_dichoptic - mua_baselinemedian_dioptic;

figure
plot(sdftm,mua_blSubAvg_dioptic)
hold on
plot(sdftm,mua_blSubAvg_dichoptic)

% Physical Alternation


% Same Stim, different history, preferred vs null


