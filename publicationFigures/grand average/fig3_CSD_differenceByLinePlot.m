%% bmcBRFS_CSD_grandAverage_loadMAT
% Are there observable differences between trial-types with LFP CSD?
% initialize variables

%% Setup
disp('start time')
datetime
clearvars -except LFP_trials
workingPC = 'office'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures\grand average';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir\CSDFigs';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures\grand average';
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


%% plot dioptic vs dichoptic
%% plot dioptic vs dichoptic
% Set smoothing factor (adjust if needed)
smoothing_factor = 20;  % You can tweak this value for more/less smoothing

% half-wave rectify the data
% halfWaveRectify = averageCSDMatrix_BRFS;
halfWaveRectify = averageCSDMatrix_diopDichop;
halfWaveRectify(halfWaveRectify > 0) = 0;
rectifiedCSD = abs(halfWaveRectify);

% Average over the recording sessions
grandAverageCSD = median(rectifiedCSD, 4, "omitmissing"); % Average across sessions

% Extract the data for dioptic and dichoptic conditions
dioptic_data = grandAverageCSD(:, :, 1);  % Dioptic condition
dichoptic_data = grandAverageCSD(:, :, 2);  % Dichoptic condition


% Smoothing the data before averaging
dioptic_data_smooth = smoothdata(dioptic_data, 'gaussian', smoothing_factor);
dichoptic_data_smooth = smoothdata(dichoptic_data, 'gaussian', smoothing_factor);

% Compute averages for each laminar compartment (across channels)
supragranular_avg_diop = mean(dioptic_data_smooth(:, 1:5), 2);   % Channels 1:5 for dioptic
supragranular_avg_dicho = mean(dichoptic_data_smooth(:, 1:5), 2);  % Channels 1:5 for dichoptic

granular_avg_diop = mean(dioptic_data_smooth(:, 6:10), 2);  % Channels 6:10 for dioptic
granular_avg_dicho = mean(dichoptic_data_smooth(:, 6:10), 2);  % Channels 6:10 for dichoptic

infragranular_avg_diop = mean(dioptic_data_smooth(:, 11:15), 2);  % Channels 11:15 for dioptic
infragranular_avg_dicho = mean(dichoptic_data_smooth(:, 11:15), 2);  % Channels 11:15 for dichoptic

% Difference between conditions (before filtering)
supragranular_diff = supragranular_avg_diop - supragranular_avg_dicho;
granular_diff = granular_avg_diop - granular_avg_dicho;
infragranular_diff = infragranular_avg_diop - infragranular_avg_dicho;

% Define the cutoff frequency and sample rate for the low-pass filter
fs = 1000; % Assuming 1000 Hz sampling rate
cutoff_freq = 30; % 20 Hz cutoff frequency
[b, a] = butter(4, cutoff_freq / (fs / 2), 'low'); % 4th-order Butterworth filter

% Apply low-pass filter to the averaged data
supragranular_avg_diop_filtered = filtfilt(b, a, supragranular_avg_diop);
supragranular_avg_dicho_filtered = filtfilt(b, a, supragranular_avg_dicho);

granular_avg_diop_filtered = filtfilt(b, a, granular_avg_diop);
granular_avg_dicho_filtered = filtfilt(b, a, granular_avg_dicho);

infragranular_avg_diop_filtered = filtfilt(b, a, infragranular_avg_diop);
infragranular_avg_dicho_filtered = filtfilt(b, a, infragranular_avg_dicho);

% Apply low-pass filter to the difference between conditions
supragranular_diff_filtered = filtfilt(b, a, supragranular_diff);
granular_diff_filtered = filtfilt(b, a, granular_diff);
infragranular_diff_filtered = filtfilt(b, a, infragranular_diff);

% Create a modern muted color palette that is easier to differentiate
colors = [
    0.121, 0.466, 0.705; % Blue
    0.682, 0.780, 0.902; % Light Blue
    0.890, 0.102, 0.110; % Red
    0.600, 0.400, 0.000; % Olive Green
    0.659, 0.661, 0.661; % Gray
];

% Create figure for line plots
f = figure;
tm = -200:1800;
% First plot: Supragranular
subplot(3, 2, 1);
plot(tm, supragranular_avg_diop_filtered, 'Color', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Dioptic'); hold on;
plot(tm, supragranular_avg_dicho_filtered, 'Color', colors(2,:), 'LineWidth', 1.5, 'DisplayName', 'Dichoptic');
xlabel('Time (ms)');
ylabel('Supragranular');
legend;
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
box('off')

% Second plot: Granular
subplot(3, 2, 3);
plot(tm, granular_avg_diop_filtered, 'Color', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Dioptic'); hold on;
plot(tm, granular_avg_dicho_filtered, 'Color', colors(2,:), 'LineWidth', 1.5, 'DisplayName', 'Dichoptic');
xlabel('Time (ms)');
ylabel('Granular');
legend;
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
box('off')


% Third plot: Infragranular
subplot(3, 2, 5);
plot(tm, infragranular_avg_diop_filtered, 'Color', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Dioptic'); hold on;
plot(tm, infragranular_avg_dicho_filtered, 'Color', colors(2,:), 'LineWidth', 1.5, 'DisplayName', 'Dichoptic');
xlabel('Time (ms)');
ylabel({'Infragranular','nA/mm^3'});
legend;
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
box('off')


% Fourth plot: Difference between conditions
subplot(3, 2, [2 4 6]);
plot(tm, supragranular_diff_filtered, 'Color', colors(4,:), 'LineWidth', 1.5, 'DisplayName', 'Supragranular Diff'); hold on;
plot(tm, granular_diff_filtered, 'Color', colors(3,:), 'LineWidth', 1.5, 'DisplayName', 'Granular Diff');
plot(tm, infragranular_diff_filtered, 'Color', colors(5,:), 'LineWidth', 1.5, 'DisplayName', 'Infragranular Diff');
xlabel('Time (ms)');
ylabel('Difference (nA/mm^3)');
title('Difference between Conditions');
legend;
xline(0, '--', 'Stimulus onset', 'LineWidth', 2);
xline(800, '--', 'Stimulus offset', 'LineWidth', 2);
hline(0)
box('off')




%% plot BRFS
% % % grandAverageCSD_BRFS = median(averageCSDMatrix_BRFS,4,"omitmissing"); % average across penetration
% % % filterCSD_BRFSps = filterCSD(grandAverageCSD_BRFS(801:2001,:,1)',.05);
% % % filterCSD_BRFSns = filterCSD(grandAverageCSD_BRFS(801:2001,:,2)',.05);
% % % 
% % % % Visualize the CSD matrix
% % % ax(6) = subplot(2,3,4);
% % % imagesc(-200:1000,1:15,filterCSD_BRFSps);
% % % colormap(ax(6),newcmap)
% % % clim([-1500 1500]);
% % % % xlabel('Time (ms)');
% % % ylabel('channel');
% % % cb = colorbar(); 
% % % % ylabel(cb,'(nA/mm)^3','FontSize',12)
% % % xl = xline(0,'--','Stimulus onset','LineWidth',3);
% % % xl = xline(800,'--','Stimulus offset','LineWidth',3);
% % % title('Preferred stimulus BRFS flash');
% % % box('off')
% % % set(gca,'XTick',[])
% % % 
% % % 
% % % 
% % % 
% % % ax(7) = subplot(2,3,5);
% % % imagesc(-200:1000,1:15,filterCSD_BRFSns);
% % % colormap(ax(7),newcmap)
% % % clim([-1500 1500]);
% % % xlabel('Time (ms)');
% % % % ylabel('channel');
% % % cb = colorbar(); 
% % % % ylabel(cb,'(nA/mm)^3','FontSize',12)
% % % xl = xline(0,'--','Stimulus onset','LineWidth',3);
% % % xl = xline(800,'--','Stimulus offset','LineWidth',3);
% % % title('Non-preferred stimulus BRFS flash');
% % % box('off')
% % % 
% % % 
% % % % Visualize the raw diff
% % % differenceMatrix_BRFS = ...
% % %     abs(grandAverageCSD_BRFS(801:2001,:,1)')-abs(grandAverageCSD_BRFS(801:2001,:,2)');
% % % filterCSD_BRFSDiff =  filterCSD(differenceMatrix_BRFS,.05);
% % % 
% % % ax(8) = subplot(2,3,6);
% % % imagesc(-200:1000,1:15,filterCSD_BRFSDiff);
% % % colormap(ax(8),'bone');
% % % clim([-500 500]);
% % % % xlabel('Time (ms)');
% % % cb = colorbar(); 
% % % ylabel(cb,'(nA/mm)^3','FontSize',12)
% % % xl = xline(0,'--','Stimulus onset','LineWidth',3);
% % % xl = xline(800,'--','Stimulus offset','LineWidth',3);
% % % title('|Difference|: |PS-NS|');
% % % box('off')
% % % set(gca,'XTick',[])
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % sgtitle('CSD Penetration Average')


%% Save output
%save fig
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
        disp('alright, saving figure to plotdir')
        cd(plotDir)
        saveName = strcat('CSDlinePlot.png');
        saveas(f,saveName) 
        saveName = strcat('CSDlinePlot.svg');
        saveas(f,saveName)
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end
