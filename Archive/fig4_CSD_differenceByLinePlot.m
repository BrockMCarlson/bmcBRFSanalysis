%% bmcBRFS_CSD_grandAverage_loadMAT
% Are there observable differences between trial-types with LFP CSD?
% initialize variables

%% Setup
disp('start time')
datetime
clearvars -except LFP_trials
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir\CSDFigs';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\CSDFigs';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)
cd(dataDir)

if ~exist('LFP_trials','var')
    tic
    load('LFP_trials.mat') % format is LFP_trials{penetration,1}{cond,1}{trial,flash}
    toc
end



%% For loop
averageCSDMatrix_diopDichop = nan(2001,15,2,31);
averageCSDMatrix_BRFS       = nan(2001,15,2,31);
for penetration = 1:size(LFP_trials,1)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');

    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',fileToLoad))
        continue
    end
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    if strcmp(string(officLamAssign.ChToUse(penetration)),"1:32")
        if any(v1Ch > 32) || any(v1Ch < 1)
            warning('skipping session without full column for now')
            disp(fileToLoad)
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('skipping session without full column for now')
            disp(fileToLoad)
            continue
        end
    end
    
    %% Obtain Monocular preference
    % bl Sub at average level (better for plotting)
    clear array_ofMonoc1 array_ofMonoc2 array_ofMonoc3 array_ofMonoc4
    
    % Monocular
    monoc_1 = [5, 11, 13, 19]; % PO RightEye
    monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    
    % convert from cell to double and combine monocular conditions
    count = 0;
    for cond = monoc_1
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:size(LFP_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(LFP_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
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
    tuned = nan(length(v1Ch),1);
    prefMonoc = nan(length(v1Ch),1);
    
    % FOR each individual unit 
    for i = 1:length(v1Ch)
        % preallocate y based on trial count
        y = nan(maxTrls,4);
        % create an array of each trial's median response for each monoc 
        % condition (trl x 4 monoc)
        y(1:size(array_ofMonoc1,3),1) = median(squeeze(array_ofMonoc1(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc2,3),2) = median(squeeze(array_ofMonoc2(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc3,3),3) = median(squeeze(array_ofMonoc3(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc4,3),4) = median(squeeze(array_ofMonoc4(200:450,i,:)),1)'; % median of each trial after stim onset
        % Now we perform baseline subtraction
        % First we get the blAvg for this contact across all trials
        y_bl(1,1) = median(median(squeeze(array_ofMonoc1(100:200,i,:)),1)');
        y_bl(1,2) = median(median(squeeze(array_ofMonoc2(100:200,i,:)),1)');
        y_bl(1,3) = median(median(squeeze(array_ofMonoc3(100:200,i,:)),1)');
        y_bl(1,4) = median(median(squeeze(array_ofMonoc4(100:200,i,:)),1)');
        y_blSub = y - median(y_bl);
        p = anova1(y_blSub,[],'off');
        if p < .05
            tuned(i,1) = true;
        else
            tuned(i,1) = false;
        end
        % Now we find the maximum response
        [M,maxRespIdx] = max(median(y_blSub,1,"omitmissing"));
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
   

    for i = 1:length(v1Ch)
        % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
        % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
        % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
        % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
        if prefMonoc(i,1) == 1
            prefCondOnFlash(i,1) = 12;
            nullCondOnFlash(i,1) = 11;
        elseif prefMonoc(i,1) == 2
            prefCondOnFlash(i,1) = 9;
            nullCondOnFlash(i,1) = 10;
        elseif prefMonoc(i,1) == 3
            prefCondOnFlash(i,1) = 10;
            nullCondOnFlash(i,1) = 9;
        elseif prefMonoc(i,1) == 4
            prefCondOnFlash(i,1) = 11;
            nullCondOnFlash(i,1) = 12;
        end

    end
        
    overallPref = mode(prefCondOnFlash);
    if overallPref == 9
        overallNull = 10;
    elseif overallPref == 10
        overallNull = 9;
    elseif overallPref == 11
        overallNull = 12;
    elseif overallPref == 12
        overallNull = 11;
    end

    %% CSD matrix across channels
    chForCSDcalc = v1Ch(1)-1:v1Ch(end)+1;
    % Dioptic vs dichoptic
    count = 0;
    for conditionNumber = [1 3]        
        % Get the number of trials for the chosen condition
        numTrials = size(LFP_trials{penetration,1}{conditionNumber,1},1);
        CSDbinocOut = nan(2001,length(v1Ch),numTrials);
        for trl = 1:numTrials
            LFPbinocTrl = LFP_trials{penetration,1}{conditionNumber,1}{trl,1}(:,chForCSDcalc); % adding ch above and below for CSD calc (to use as vaknin pad)
            blLFPbinoc = mean(LFPbinocTrl(100:200,:),1);
            LFPbinocBlSub = LFPbinocTrl-blLFPbinoc;
            CSDbinocTrl = calcCSD_classic(LFPbinocBlSub);
            CSDbinocOut(:,:,trl) = CSDbinocTrl(:,2:16); % limit to origional V1 ch lim
        end
        % Average across trials and save output
        count = count + 1; % for pref vs null
        averageCSDMatrix_diopDichop(:,:,count,penetration) = median(CSDbinocOut,3); % Average across trl. averageCSDMatrix is (ch1 x ch2 x cond x penetration)
    end

    % BRFS pref vs null
    count = 0;
    for conditionNumber = [overallPref overallNull]        
        % Get the number of trials for the chosen condition
        numTrials = size(LFP_trials{penetration,1}{conditionNumber,1},1);
        CSDflashOut = nan(2001,length(v1Ch),numTrials);
        for trl = 1:numTrials
            LFPflashTrl = LFP_trials{penetration,1}{conditionNumber,1}{trl,1}(:,chForCSDcalc);
            blLFPflash = mean(LFPflashTrl(100:200,:),1);
            LFPflashBlSub = LFPflashTrl-blLFPflash;
            CSDflashTrl = calcCSD_classic(LFPflashBlSub);
            CSDflashOut(:,:,trl) = CSDflashTrl(:,2:16);
        end
        % Average across trials and save output
        count = count + 1; % for pref vs null
        averageCSDMatrix_BRFS(:,:,count,penetration) = median(CSDflashOut,3); % Average across trl. averageCSDMatrix is (ch1 x ch2 x cond x penetration)
    end



    disp(strcat('Done with file number: ',string(penetration)))
end


%% plot dioptic vs dichoptic: averageCSDMatrix_diopDichop

% Organize the data
data_diop = squeeze(averageCSDMatrix_diopDichop(:, :, 1, :));  % Dioptic condition
data_dichop = squeeze(averageCSDMatrix_diopDichop(:, :, 2, :));  % Dichoptic condition

% half-wave rectify the data
halfWave_diop = data_diop;
halfWave_diop(halfWave_diop >= 0) = NaN;
halfWaveRectifiedCSD_diop = abs(halfWave_diop);

halfWave_dichop = data_dichop;
halfWave_dichop(halfWave_dichop >= 0) = NaN;
halfWaveRectifiedCSD_dichop = abs(halfWave_dichop);

% Compute averages for each laminar compartment (across channels)
supragranular_avg_diop = squeeze(mean(halfWaveRectifiedCSD_diop(:, 1:5, :), 2,"omitmissing"));   % Channels 1:5 for dioptic
supragranular_avg_dicho = squeeze(mean(halfWaveRectifiedCSD_dichop(:, 1:5, :), 2,"omitmissing"));  % Channels 1:5 for dichoptic

granular_avg_diop = squeeze(mean(halfWaveRectifiedCSD_diop(:, 6:10, :), 2,"omitmissing"));  % Channels 6:10 for dioptic
granular_avg_dicho = squeeze(mean(halfWaveRectifiedCSD_dichop(:, 6:10, :), 2,"omitmissing"));  % Channels 6:10 for dichoptic

infragranular_avg_diop = squeeze(mean(halfWaveRectifiedCSD_diop(:, 11:15, :), 2,"omitmissing"));  % Channels 11:15 for dioptic
infragranular_avg_dicho = squeeze(mean(halfWaveRectifiedCSD_dichop(:, 11:15, :), 2,"omitmissing"));  % Channels 11:15 for dichoptic

% Percent change from baseline
bl_diop_s = mean(supragranular_avg_diop(1:200,:),1,"omitmissing");
percentChangeFromBl_diop_s = ((supragranular_avg_diop-bl_diop_s)./bl_diop_s)*100; % percent change from baseline
bl_dichop_s = mean(supragranular_avg_dicho(1:200,:),1,"omitmissing");
percentChangeFromBl_dichop_s = ((supragranular_avg_dicho-bl_dichop_s)./bl_dichop_s)*100; % percent change from baseline

bl_diop_g = mean(granular_avg_diop(1:200,:),1,"omitmissing");
percentChangeFromBl_diop_g = ((granular_avg_diop-bl_diop_g)./bl_diop_g)*100; % percent change from baseline
bl_dichop_g = mean(granular_avg_dicho(1:200,:),1,"omitmissing");
percentChangeFromBl_dichop_g = ((granular_avg_dicho-bl_dichop_g)./bl_dichop_g)*100; % percent change from baseline

bl_diop_i = mean(infragranular_avg_diop(1:200,:),1,"omitmissing");
percentChangeFromBl_diop_i = ((infragranular_avg_diop-bl_diop_i)./bl_diop_i)*100; % percent change from baseline
bl_dichop_i = mean(infragranular_avg_dicho(1:200,:),1,"omitmissing");
percentChangeFromBl_dichop_i = ((infragranular_avg_dicho-bl_dichop_i)./bl_dichop_i)*100; % percent change from baseline

% Average over the recording sessions
grandAverageCSD_diop_s = mean(percentChangeFromBl_diop_s, 2, "omitmissing"); % Average across sessions
grandAverageCSD_dichop_s = mean(percentChangeFromBl_dichop_s, 2, "omitmissing"); % Average across sessions

grandAverageCSD_diop_g = mean(percentChangeFromBl_diop_g, 2, "omitmissing"); % Average across sessions
grandAverageCSD_dichop_g = mean(percentChangeFromBl_dichop_g, 2, "omitmissing"); % Average across sessions

grandAverageCSD_diop_i = mean(percentChangeFromBl_diop_i, 2, "omitmissing"); % Average across sessions
grandAverageCSD_dichop_i = mean(percentChangeFromBl_dichop_i, 2, "omitmissing"); % Average across sessions

% Difference between conditions 
supragranular_diff = grandAverageCSD_diop_s - grandAverageCSD_dichop_s;
granular_diff = grandAverageCSD_diop_g - grandAverageCSD_dichop_g;
infragranular_diff = grandAverageCSD_diop_i - grandAverageCSD_dichop_i;

% Low pass filter the data
fs = 1000; % Assuming 1000 Hz sampling rate
cutoff_freq = 45; % 20 Hz cutoff frequency
[b, a] = butter(4, cutoff_freq / (fs / 2), 'low'); % 4th-order Butterworth filter

% Apply low-pass filter to the averaged data
supragranular_avg_diop_filtered = filtfilt(b, a, grandAverageCSD_diop_s);
supragranular_avg_dicho_filtered = filtfilt(b, a, grandAverageCSD_dichop_s);

granular_avg_diop_filtered = filtfilt(b, a, grandAverageCSD_diop_g);
granular_avg_dicho_filtered = filtfilt(b, a, grandAverageCSD_dichop_g);

infragranular_avg_diop_filtered = filtfilt(b, a, grandAverageCSD_diop_i);
infragranular_avg_dicho_filtered = filtfilt(b, a, grandAverageCSD_dichop_i);

% Apply low-pass filter to the difference between conditions
cutoff_freq = 10; % 20 Hz cutoff frequency
[b, a] = butter(4, cutoff_freq / (fs / 2), 'low'); % 4th-order Butterworth filter
supragranular_diff_filtered = filtfilt(b, a, supragranular_diff);
granular_diff_filtered = filtfilt(b, a, granular_diff);
infragranular_diff_filtered = filtfilt(b, a, infragranular_diff);

% Create a modern muted color palette that is easier to differentiate
% Whole number RGB values (0 to 255)
colors_whole = [
    179, 203, 205; % Sky Dioptic/PS
    148, 110, 36;  % Oak Dichoptic/NS
    236, 183, 72;  % Highlight -- Gran
    139, 161, 142; % Sage -- SupraGran
    119, 119, 119; % Dark Gray Diff - Infragran
];

% Convert to decimal format (0 to 1)
colors = colors_whole / 255;

% Create figure for line plots
f = figure;
set(f,"Position",[-1861 -48 1734 866])
set(f,'defaultLegendAutoUpdate','off')

tm = -200:1800;
% First plot: Supragranular
subplot(3, 2, 1);
plot(tm, supragranular_avg_diop_filtered, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'Dioptic'); hold on;
plot(tm, supragranular_avg_dicho_filtered, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'Dichoptic');
xlabel('Time (ms)');
ylabel('Supragranular');
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
box('off')

% Second plot: Granular
subplot(3, 2, 3);
plot(tm, granular_avg_diop_filtered, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'Dioptic'); hold on;
plot(tm, granular_avg_dicho_filtered, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'Dichoptic');
xlabel('Time (ms)');
ylabel('Granular');
legend({'Dioptic','Dichoptic'});
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
box('off')


% Third plot: Infragranular
subplot(3, 2, 5);
plot(tm, infragranular_avg_diop_filtered, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'Dioptic'); hold on;
plot(tm, infragranular_avg_dicho_filtered, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'Dichoptic');
xlabel('Time (ms)');
ylabel({'Infragranular','% Change'});
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
box('off')


% Fourth plot: Difference between conditions
subplot(3, 2, [2 4 6]);
plot(tm, supragranular_diff_filtered, 'Color', colors(4,:), 'LineWidth', 3, 'DisplayName', 'Supragranular Diff'); hold on;
plot(tm, granular_diff_filtered, 'Color', colors(3,:), 'LineWidth', 3, 'DisplayName', 'Granular Diff');
plot(tm, infragranular_diff_filtered, 'Color', colors(5,:), 'LineWidth', 3, 'DisplayName', 'Infragranular Diff');
xlabel('Time (ms)');
ylabel('Difference in % Change');
title('Difference between Conditions');
legend({'Supragranular','Granular','Infragranular'});
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
yline(0, '--', 'LineWidth', 2);
xlim([-200 1800])
box('off')


answer = questdlg('Would you like to save this figure? _diopticDichoptic', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
        disp('alright, saving figure to plotdir')
        cd(plotDir)
        saveName = strcat('CSDlinePlot_diopticDichoptic.png');
        saveas(f,saveName) 
        saveName = strcat('CSDlinePlot_diopticDichoptic.svg');
        saveas(f,saveName)
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end




%% plot BRFS: averageCSDMatrix_BRFS

% Organize the data
data_PS = squeeze(averageCSDMatrix_BRFS(:, :, 1, :));  % Dioptic condition
data_NS = squeeze(averageCSDMatrix_BRFS(:, :, 2, :));  % Dichoptic condition

% half-wave rectify the data
halfWave_PS = data_PS;
halfWave_PS(halfWave_PS >= 0) = NaN;
halfWaveRectifiedCSD_PS = abs(halfWave_PS);

halfWave_NS = data_NS;
halfWave_NS(halfWave_NS >= 0) = NaN;
halfWaveRectifiedCSD_NS = abs(halfWave_NS);

% Compute averages for each laminar compartment (across channels)
supragranular_avg_PS = squeeze(mean(halfWaveRectifiedCSD_PS(:, 1:5, :), 2,"omitmissing"));   % Channels 1:5 for dioptic
supragranular_avg_NS = squeeze(mean(halfWaveRectifiedCSD_NS(:, 1:5, :), 2,"omitmissing"));  % Channels 1:5 for dichoptic

granular_avg_PS = squeeze(mean(halfWaveRectifiedCSD_PS(:, 6:10, :), 2,"omitmissing"));  % Channels 6:10 for dioptic
granular_avg_NS = squeeze(mean(halfWaveRectifiedCSD_NS(:, 6:10, :), 2,"omitmissing"));  % Channels 6:10 for dichoptic

infragranular_avg_PS = squeeze(mean(halfWaveRectifiedCSD_PS(:, 11:15, :), 2,"omitmissing"));  % Channels 11:15 for dioptic
infragranular_avg_NS = squeeze(mean(halfWaveRectifiedCSD_NS(:, 11:15, :), 2,"omitmissing"));  % Channels 11:15 for dichoptic

% Percent change from baseline
bl_PS_s = mean(supragranular_avg_PS(1:200,:),1,"omitmissing");
percentChangeFromBl_PS_s = ((supragranular_avg_PS-bl_PS_s)./bl_PS_s)*100; % percent change from baseline
bl_NS_s = mean(supragranular_avg_NS(1:200,:),1,"omitmissing");
percentChangeFromBl_NS_s = ((supragranular_avg_NS-bl_NS_s)./bl_NS_s)*100; % percent change from baseline

bl_PS_g = mean(granular_avg_PS(1:200,:),1,"omitmissing");
percentChangeFromBl_PS_g = ((granular_avg_PS-bl_PS_g)./bl_PS_g)*100; % percent change from baseline
bl_NS_g = mean(granular_avg_NS(1:200,:),1,"omitmissing");
percentChangeFromBl_NS_g = ((granular_avg_NS-bl_NS_g)./bl_NS_g)*100; % percent change from baseline

bl_PS_i = mean(infragranular_avg_PS(1:200,:),1,"omitmissing");
percentChangeFromBl_PS_i = ((infragranular_avg_PS-bl_PS_i)./bl_PS_i)*100; % percent change from baseline
bl_NS_i = mean(infragranular_avg_NS(1:200,:),1,"omitmissing");
percentChangeFromBl_NS_i = ((infragranular_avg_NS-bl_NS_i)./bl_NS_i)*100; % percent change from baseline

% Average over the recording sessions
grandAverageCSD_PS_s = mean(percentChangeFromBl_PS_s, 2, "omitmissing"); % Average across sessions
grandAverageCSD_NS_s = mean(percentChangeFromBl_NS_s, 2, "omitmissing"); % Average across sessions

grandAverageCSD_PS_g = mean(percentChangeFromBl_PS_g, 2, "omitmissing"); % Average across sessions
grandAverageCSD_NS_g = mean(percentChangeFromBl_NS_g, 2, "omitmissing"); % Average across sessions

grandAverageCSD_PS_i = mean(percentChangeFromBl_PS_i, 2, "omitmissing"); % Average across sessions
grandAverageCSD_NS_i = mean(percentChangeFromBl_NS_i, 2, "omitmissing"); % Average across sessions

% Difference between conditions 
supragranular_diff = grandAverageCSD_PS_s - grandAverageCSD_NS_s;
granular_diff = grandAverageCSD_PS_g - grandAverageCSD_NS_g;
infragranular_diff = grandAverageCSD_PS_i - grandAverageCSD_NS_i;

% Low pass filter the data
fs = 1000; % Assuming 1000 Hz sampling rate
cutoff_freq = 45; % 20 Hz cutoff frequency
[b, a] = butter(4, cutoff_freq / (fs / 2), 'low'); % 4th-order Butterworth filter

% Apply low-pass filter to the averaged data
supragranular_avg_PS_filtered = filtfilt(b, a, grandAverageCSD_PS_s);
supragranular_avg_NS_filtered = filtfilt(b, a, grandAverageCSD_NS_s);

granular_avg_PS_filtered = filtfilt(b, a, grandAverageCSD_PS_g);
granular_avg_NS_filtered = filtfilt(b, a, grandAverageCSD_NS_g);

infragranular_avg_PS_filtered = filtfilt(b, a, grandAverageCSD_PS_i);
infragranular_avg_NS_filtered = filtfilt(b, a, grandAverageCSD_NS_i);

% Apply low-pass filter to the difference between conditions
cutoff_freq = 10; % 20 Hz cutoff frequency
[b, a] = butter(4, cutoff_freq / (fs / 2), 'low'); % 4th-order Butterworth filter
supragranular_diff_filtered = filtfilt(b, a, supragranular_diff);
granular_diff_filtered = filtfilt(b, a, granular_diff);
infragranular_diff_filtered = filtfilt(b, a, infragranular_diff);

% Create a modern muted color palette that is easier to differentiate
% Whole number RGB values (0 to 255)
colors_whole = [
    179, 203, 205; % Sky Dioptic/PS
    148, 110, 36;  % Oak Dichoptic/NS
    236, 183, 72;  % Highlight -- Gran
    139, 161, 142; % Sage -- SupraGran
    119, 119, 119; % Dark Gray Diff - Infragran
];

% Convert to decimal format (0 to 1)
colors = colors_whole / 255;

% Create figure for line plots
f = figure;
set(f,"Position",[-1861 -48 1734 866])
set(f,'defaultLegendAutoUpdate','off')

tm = -200:1800;
% First plot: Supragranular
subplot(3, 2, 1);
plot(tm, supragranular_avg_PS_filtered, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'PS'); hold on;
plot(tm, supragranular_avg_NS_filtered, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'NS');
xlabel('Time (ms)');
ylabel('Supragranular');
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(800, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
box('off')

% Second plot: Granular
subplot(3, 2, 3);
plot(tm, granular_avg_PS_filtered, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'PS'); hold on;
plot(tm, granular_avg_NS_filtered, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'NS');
xlabel('Time (ms)');
ylabel('Granular');
legend({'PS','NS'});
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(800, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
box('off')


% Third plot: Infragranular
subplot(3, 2, 5);
plot(tm, infragranular_avg_PS_filtered, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'PS'); hold on;
plot(tm, infragranular_avg_NS_filtered, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'NS');
xlabel('Time (ms)');
ylabel({'Infragranular','% Change'});
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(800, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
box('off')


% Fourth plot: Difference between conditions
subplot(3, 2, [2 4 6]);
plot(tm, supragranular_diff_filtered, 'Color', colors(4,:), 'LineWidth', 3, 'DisplayName', 'Supragranular Diff'); hold on;
plot(tm, granular_diff_filtered, 'Color', colors(3,:), 'LineWidth', 3, 'DisplayName', 'Granular Diff');
plot(tm, infragranular_diff_filtered, 'Color', colors(5,:), 'LineWidth', 3, 'DisplayName', 'Infragranular Diff');
xlabel('Time (ms)');
ylabel('Difference in % Change');
title('Difference between Conditions');
legend({'Supragranular','Granular','Infragranular'});
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(800, '--', 'Stimulus onset', 'LineWidth', 2);
xline(1600, '--', 'Stimulus offset', 'LineWidth', 2);
yline(0, '--', 'LineWidth', 2);
xlim([-200 1800])
box('off')


answer = questdlg('Would you like to save this figure? _PSticNS', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
        disp('alright, saving figure to plotdir')
        cd(plotDir)
        saveName = strcat('fig4_CSDlinePlot_PSticNS.png');
        saveas(f,saveName) 
        saveName = strcat('fig4_CSDlinePlot_PSticNS.svg');
        saveas(f,saveName)
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end
