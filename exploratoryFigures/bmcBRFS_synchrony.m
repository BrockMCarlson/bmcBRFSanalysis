% bmcBRFS_coherence.m
% Brock Carlson - 1/9/24
% The goal of this code is to analyze the bmcBRFS data in the same manner
% as Fries, Roelfsema, Engel, Konig, & Singer (1997) PNAS, which suggests
% that differential changes in oscilitory patterning at the early stages of
% visual processing determines which signals are perceived, as opposed to
% neural discharge rates.
% CCH - cross-correlation histogram
% RMA - relative modulation amplitude
% STA - spike triggered averages
% SFC - spike-field coherence.
% We aim to first perform this analysis on individual units, individual
% sessions, laminar compartments, and then full data averages.


%% Setup
clear
% Directories
path1 = strcat(getenv('HOMEDRIVE'),getenv("HOMEPATH"),filesep);
path2 = 'Documents\GitHub\bmcBRFSanalysis';
codeDir = strcat(path1,path2);
cd(codeDir)
path3 = 'Documents\MATLAB\formattedDataOutputs';
plotDir = strcat(path1,filesep,path3);
%sortedDataDir -- need to run new sorted data on home pc
sortedDataDir = 'S:\bmcBRFS_sortedData_Nov23';

%% Choose inidividual file to explore
cd(sortedDataDir)
load('sortedData_211008_B_bmcBRFS001.mat')

%% Import laminar assignments
officLamAssign = importLaminarAssignments( path1 + ...
    "Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx",...
    "Sheet1", [2, Inf]);

granBtm = officLamAssign.Probe11stFold4c(19); % channel corresponding to the bottom of layer 4c


% Calculate V1 ch boundaries
v1Top = granBtm - 9;
v1Btm = granBtm + 5;
v1Ch = v1Top:v1Btm;
% limit ch to cortex only
columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
if any(v1Ch > 32) || any(v1Ch < 1)
    columnNames = columnNames((v1Ch >= 1) & (v1Ch<=32));
    v1Ch        = v1Ch((v1Ch >= 1) & (v1Ch<=32));
end
    
    
%% test for monocular preference
ch = v1Ch;
clear array_ofMonoc1 array_ofMonoc2 array_ofMonoc3 array_ofMonoc4

% Monocular
monoc_1 = [5, 11, 13, 19]; % PO RightEye
monoc_2 = [8, 10, 16, 18]; % PO LeftEye
monoc_3 = [7, 9, 15, 17];  % NPO RightEye
monoc_4 = [6, 12, 14, 20]; % NPO LeftEye

% convert from cell to double and combine monocular conditions
warning('taking abs() of IDX value')
count = 0;
for cond = monoc_1
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc1(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); % 1000 x 32
    end
end
clear cond count trl 
count = 0;
for cond = monoc_2
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc2(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); 
    end
end
clear cond count trl 
count = 0;
for cond = monoc_3
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc3(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); 
    end
end
clear cond count trl 
count = 0;
for cond = monoc_4
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc4(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); 
    end
end
clear cond count trl 


% Test for monocular preference with anova
% The anova function ignores NaN values, <undefined> values, empty 
% characters, and empty strings in y. If factors or tbl contains NaN or 
% <undefined> values, or empty characters or strings, the function ignores 
% the corresponding observations in y. The ANOVA is balanced if each factor 
% value has the same number of observations after the function disregards 
% empty or NaN values. Otherwise, the function performs an unbalanced 
% ANOVA. (We will be doing an unbalanced ANOVA).

% preallocate array with nan - running an unbalanced ANOVA)
A = [size(array_ofMonoc1,3),...
    size(array_ofMonoc2,3),...
    size(array_ofMonoc3,3),...
    size(array_ofMonoc4,3)];
maxTrls = max(A);
minTrls = min(A);
tuned = nan(length(ch),1);
prefMonoc = nan(length(ch),1);

% FOR each individual unit 
for i = 1:length(ch)
    % preallocate y based on trial count
    y = nan(maxTrls,4);
    % create an array of each trial's median response for each monoc 
    % condition (trl x 4 monoc)
    y(1:size(array_ofMonoc1,3),1) = median(squeeze(array_ofMonoc1(200:800,i,:)),1)'; % median of each trial after stim onset
    y(1:size(array_ofMonoc2,3),2) = median(squeeze(array_ofMonoc2(200:800,i,:)),1)'; % median of each trial after stim onset
    y(1:size(array_ofMonoc3,3),3) = median(squeeze(array_ofMonoc3(200:800,i,:)),1)'; % median of each trial after stim onset
    y(1:size(array_ofMonoc4,3),4) = median(squeeze(array_ofMonoc4(200:800,i,:)),1)'; % median of each trial after stim onset
    % Now we perform baseline subtraction
    % First we get the blAvg for this contact across all trials
    y_bl(1,1) = median(median(squeeze(array_ofMonoc1(100:200,i,:)),1)');
    y_bl(1,2) = median(median(squeeze(array_ofMonoc2(100:200,i,:)),1)');
    y_bl(1,3) = median(median(squeeze(array_ofMonoc3(100:200,i,:)),1)');
    y_bl(1,4) = median(median(squeeze(array_ofMonoc4(100:200,i,:)),1)');
    y_blSub = abs(y - median(y_bl));
    p = anova1(y_blSub,[],'off');
    if p < .05
        tuned(i,1) = true;
    else
        tuned(i,1) = false;
    end
    % Now we find the maximum response
    [M,maxRespIdx] = max(median(y_blSub,1,"omitnan"));
    prefMonoc(i,1) = maxRespIdx;
    % And assign the null condition
    if maxRespIdx == 1
        nullRespIdx = 4;
    elseif maxRespIdx == 2
        nullRespIdx = 3;
    elseif maxRespIdx == 3
        nullRespIdx = 2;
    elseif maxRespIdx == 4
        nullRespIdx = 1;
    end
    nullMonoc(i,1) = nullRespIdx;
end

% Plot results
%trl avg
ps_avg = mean(array_ofMonoc1,3);
ps_table = array2table(ps_avg);
ns_avg = mean(array_ofMonoc4,3);
ns_table = array2table(ns_avg);

%Now convert ps_avg (a double array) into a table (input to timetable must
%be a table)
% The goal is to have 32 variables, each as a column, representing a
% different depth. 
Time = milliseconds(-200):milliseconds(1):milliseconds(1000);
ps_TT = table2timetable(ps_table,'RowTimes',Time);
ps_TT2 = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
ns_TT = table2timetable(ns_table,'RowTimes',Time);
ns_TT2 = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);

% stackedplot()
f1=figure;
stackedplot(ps_TT2,ns_TT2)
title('monocularPreference')

%% test auto-correlation and cross-correlation
% units on ch 9 and ch10 look good... lets look at those.
% Preferred stimulus is % PO RightEye: BRFS condition 12
% Null stimulus is NPO LeftEye: BRFS condition11
f2 = figure;
cond = 12;
trl = 1;
ch = 10;
x  = IDX(cond).LFP_gamma1{trl,2}(:,ch);
[c,lags] = xcorr(x,80,'normalized');
stem(lags,c)
xlabel('lag in ms')
ylabel('Auto-Correlation')
title('LFP gamma-band auto-correlation for ch 10')


% Fields to loop through
fields_to_plot = {'LFP_gamma1', 'LFP_alpha', 'LFP_bb', 'CSD_gamma1', 'MUAe'};

% Parameters for auto-correlation plot
max_lag = 80;

% Loop through each field
close all
for field_idx = 1:length(fields_to_plot)
    current_field = fields_to_plot{field_idx};
    
    % Create a new figure for each field
    figure
    x  = IDX(cond).(current_field){trl,1}(:,ch);
    [c,lags] = xcorr(x,80,'normalized');
    stem(lags,c)
    xlabel('lag in ms')
    ylabel('Auto-Correlation')
    title('LFP gamma-band auto-correlation for ch 10')
    % Adjust layout for better visualization
    title(['Auto-correlation for ' current_field]);
end

% There does not to be a clear linear effect accross depth, althoug there
% are differences. Perhaps similar to differences seen trial-to-trial.
% Since those are semi-stable, lets look across conditions to see if we can
% see anything else. 
close all
for cond = 1:size(IDX,1)
    figure; 
    for trl = 1:length(IDX(cond).LFP_gamma1)
        x  = IDX(cond).LFP_gamma1{trl,2}(:,10);
        y = IDX(cond).LFP_gamma1{trl,2}(:,15);
        [c(:,trl),lags] = xcorr(x,y,80,'normalized');
    end
    c_avg = mean(c,2);
    stem(lags,c_avg)
    xlabel('lag in ms')
    ylabel('Cross-correlation')
    title(IDX(cond).conditionString)
end


% % Calculate CCH (cross-correlation histogram) from MUA
% Plotted as +/- 80ms vs RMA
% Detailed methods available in detail at Konig (1994) J Neurosci. Methods
% 1. Calculate correlogram
% 2. Fit a damped cosine wave (Gabor function) to the correlogram
% 3. Ensure that the function accounts for...
%   3a. at least 15% of the data variance 
%   3b. the z-scores of significant peaks must be >2
% 
% Calculate RMA (Response modulation amplitude - a percentage) from MUA
% The strength of synchronization and the regularity of oscillations
% RMA are calculated for the central and the first satellite peak,
% respectively.
% RMA = amplitude of respective peak (measured from the offset caused by
% accidntial coincidences) divided by the offset and multiplied by 100.
% 
% % Calculate STA (Spike-triggered averages) from LFP
% Plotted as +/- 128ms vs LFP uV
% 1. band-pass filter raw data between 1-100 Hz
% 2. average LFPs within a window of +/- 128ms centered on each triggered
% spike.
% Methods detailed in Gray & Singer (1989) PNAS
% 
% 
% % Calculate SFC (spike-field coherence) 
% Plotted as freq from 0 to 102 Hz in 3.9Hz bins vs 0-0.1 SFC
% You should see peak between 39 and 63 Hz
% For each of the LFP segments used for computation of STAs...
% 1. calculate the power spectrum
% 2. Average these spectra (obtaining the spike-triggered power spectrum)
% 3. SFC (computed as a ratio) = power spectrum of STA / spike-triggered
% power spectrum.
