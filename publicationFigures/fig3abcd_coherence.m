%% fig2_coherence
% Are there observable differences between trial-types with LFP coherence?
% initialize variables
clearvars -except LFP_trials
tic
workingPC = 'home'; % options: 'home', 'office'

%% Setup
disp('start time')
datetime
if strcmp(workingPC,'home')
    codeDir     = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)


cd(dataDir)
if ~exist('LFP_trials','var')
    load('LFP_trials.mat') % format is LFP_trials{penetration,1}{cond,1}{trial,1}(2001tm x 65Ch)
end

averageCoherenceMatrix_diopDichop = nan(15,15,2,size(officLamAssign,1));

%% for loop
for penetration = 1:size(LFP_trials,1)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    penetrationFileName = strcat('trialTriggeredData_',probeName(1:19),'.mat');
    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',penetrationFileName))
        continue
    end
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
    if strcmp(string(officLamAssign.ChToUse(penetration)),"1:32")
        if any(v1Ch > 32) || any(v1Ch < 1)
            warning('skipping session without full column for now')
            disp(penetrationFileName)
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('skipping session without full column for now')
            disp(penetrationFileName)
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
    




        


    %% Coherence matrix across channels
    % Parameters for mscohere
    fs = 1000;  % Sampling frequency in Hz
    windowSize = 256;  % Window size for computing the coherence
    overlap = windowSize/2;  % Overlap between windows
    tm_full = 1:2001;
    tm_bl = 1:200;
    tm_1stOnset = 201:1000; %800ms adapter
    tm_2ndOnset = 1001:1800;
    tm_offset = 1801:2001;
    tm_coher = 1289:1800; % Time window of data. Last 512ms of trial. 
    % Dioptic vs dichoptic
    count = 0;
    for conditionNumber = [1 3]        % dioptic and dichoptic
        % Get the number of trials for the chosen condition
        numTrials = size(LFP_trials{penetration,1}{conditionNumber,1},1);
        % Initialize coherence matrix
        coherenceMatrix = nan(15,15, numTrials);
        % Loop through all trials and compute coherence for each channel pair
        for trialNumber = 1:numTrials
            for channel1 = v1Ch
                for channel2 = v1Ch(1):channel1
                    % Extract the LFP_bb data for the chosen channels and trial
                    % format is LFP_trials{penetration,1}{cond,1}{trial,flash}
                    lfpGammaData1 = LFP_trials{penetration,1}{conditionNumber,1}{trialNumber,1}(:,channel1);
                    lfpGammaData2 = LFP_trials{penetration,1}{conditionNumber,1}{trialNumber,1}(:,channel2);
                    %baseline subtract
                    bl1 = median(lfpGammaData1(tm_bl));
                    lfp_blsub_1 = lfpGammaData1 - bl1;
                    bl2 = median(lfpGammaData2(tm_bl));
                    lfp_blsub_2 = lfpGammaData2 - bl2;
    
                    % Compute coherence
                    [coherence, freq] = mscohere(lfp_blsub_1(tm_coher), lfp_blsub_2(tm_coher), windowSize, overlap, [], fs);
        
                    % Store coherence in the matrix
                    coherenceMatrix(channel1-v1Ch(1)+1, channel2-v1Ch(1)+1, trialNumber) = median(coherence(2:9));  % coherenceMatrix is (ch1 x ch2 x trialNum)
                end
            end
        end
        % Average across trials and save output
        count = count + 1; % for pref vs null
        averageCoherenceMatrix_diopDichop(:,:,count,penetration) = median(coherenceMatrix,3); % Average across trl. averageCoherenceMatrix is (ch1 x ch2 x cond x penetration)
    end

    disp(strcat('Done with file number: ',string(penetration)))
end


%% plot dioptic vs dichoptic
grandAverageCoherence_diopDichop = median(averageCoherenceMatrix_diopDichop,4,"omitmissing"); % average across penetration
  
disp('analyzing freq')
disp(freq(2:9))


% Visualize the coherence matrix
f = figure;
set(f,"Position",[-1611 318 1499 268])
ax(1) = subplot(1,4,1);
imagesc(grandAverageCoherence_diopDichop(:,:,1));
colormap(ax(1),'jet');
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Dioptic');

ax(2) = subplot(1,4,2);
imagesc(grandAverageCoherence_diopDichop(:,:,2));
colormap(ax(2),'jet');
colorbar;
xlabel('Channel');
ylabel('Channel');
title('Dichoptic');


% Difference plot 
coherenceMatrix1 = squeeze(averageCoherenceMatrix_diopDichop(:,:,1,:)); 
coherenceMatrix2 = squeeze(averageCoherenceMatrix_diopDichop(:,:,2,:));
diff_1 = grandAverageCoherence_diopDichop(:,:,1)-grandAverageCoherence_diopDichop(:,:,2);

ax(3) = subplot(1,4,3);
imagesc(diff_1);
hline(5.5)
hline(10.5)
vline(5.5)
vline(10.5)
colormap(ax(3),'bone');
clim([-.04 .04])
e = colorbar;
e.Label.String = "Coherence Difference";
e.Label.Rotation = 270;
xlabel('Channel');
ylabel('Channel');
title('difference');




%% Statistical test - ANOVA between compartment comparisons
SxS_1 = diff_1(1:5,1:5); %half block
GxS_1 = diff_1(6:10,1:5);
IxS_1 = diff_1(11:15,1:5);
GxG_1 = diff_1(6:10,6:10); % half block
IxG_1 = diff_1(11:15,6:10);
IxI_1 = diff_1(11:15,11:15); % half block

% aov = anova(y) performs a one-way ANOVA and returns the anova object...
% aov for the response data in the matrix y. Each column of y is treated...
% as a different factor value.
% construct y for ANOVA fro cross- compartment comparisons
holder_cross_1(:,1) = reshape(GxS_1,[25,1]);
holder_cross_1(:,2) = reshape(IxS_1,[25,1]);
holder_cross_1(:,3) = reshape(IxG_1,[25,1]);
[aov_cross_1,tbl,stats] = anova1(holder_cross_1,[],"off");
disp(aov_cross_1)

%ANOVA on within-compartment comparisons
SxS_1(SxS_1 == 0) = NaN;
GxG_1(GxG_1 == 0) = NaN;
IxI_1(IxI_1 == 0) = NaN;
holder_same_1(:,1) = reshape(SxS_1,[25,1]);
holder_same_1(:,2) = reshape(GxG_1,[25,1]);
holder_same_1(:,3) = reshape(IxI_1,[25,1]);
[aov_same_1]= anova1(holder_same_1,[],"off");
%  disp(aov_same_1)

%% Plotting the results of aov_cross_1 as bar plots

% Calculate means and standard errors for cross comparisons
means_cross = nanmean(holder_cross_1);
sems_cross = nanstd(holder_cross_1) ./ sqrt(sum(~isnan(holder_cross_1))); % standard error of the mean (SEM)

% Plot bar plot for cross comparisons
ax(4) = subplot(1,4,4);
b = bar(means_cross, 'FaceColor', [0.2 0.2 0.8]); % Blue color for bars
hold on;

% Add error bars
num_groups = length(means_cross);
errorbar(1:num_groups, means_cross, sems_cross, 'k', 'LineStyle', 'none', 'LineWidth', 1.5); % Add error bars

% Add significance indicators (for example, between bars 1 and 2)
x1 = 1; x2 = 3; % Positions of the bars to connect
y = max(means_cross + sems_cross) * 1.1; % Height of the line above the highest bar
plot([x1, x2], [y, y], 'k-', 'LineWidth', 1.5); % Horizontal line
text(mean([x1, x2]), y * 1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 16); % Star to indicate significance

% Customize the plot
title('Cross-Compartment Comparisons');
xticks(1:3);
xticklabels({'GxS', 'IxS', 'IxG'});
ylabel('Mean Value');
xlabel('Group');
grid on;
hold off;


%% Save output
%save fig
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
        disp('alright, saving figure to plotdir')
        sgtitle('Coherence Penetration Average')
        cd(plotDir)
        saveName = strcat('fig3_coherence_diopDichop.png');
        saveas(f,saveName) 
        saveName = strcat('fig3_coherence_diopDichop.svg');
        saveas(f,saveName) 
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end


