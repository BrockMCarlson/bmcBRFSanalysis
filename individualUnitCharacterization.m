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
    % 2. Dioptic vs dichoptic
    % 3. Same stimulus, different history, Preferred vs null
    % 4. Same stimulus, different history, mixed vs mixed.

%% Setup
clear
codeDir = 'C:\Users\neuropixel\Documents\GitHub\preProcessEphysData'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';

cd(outputDir)
load('filteredData_211008_B_bmcBRFS001.mat')
load('sortedData_211008_B_bmcBRFS001.mat')

%% Laminar plots, on condition 1, 'Simult. Dioptic. PO'
% CSD


% PSD


% MUAe


%% Continuous line plots, 

% Eye vs ori
% First 800 ms of condition 5-8,9-12, 13-16, 17-20


% Dioptic vs dichoptic


% Same Stim, different history, preferred vs null


