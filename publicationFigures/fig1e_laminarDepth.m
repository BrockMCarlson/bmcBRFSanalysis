%% fig1e_laminarDepth
% MUA Grand Average by Laminar Depth

%% Setup
disp('start time')
datetime
close all
clearvars -except MUA_trials
tm_full = -200:1800; % Time vector

workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)
cd(dataDir)

if ~exist('MUA_trials','var')
    tic
    load('MUA_trials.mat') % format is MUA_trials{penetration,1}{condition,1}{trial,flash}
    toc
end

%% Preallocate data storage for 15 channels across all penetrations
averageMUAMatrix_Diop = nan(2001, 15, 31);
averageMUAMatrix_Dichop = nan(2001, 15, 31);

for penetration = 1:size(MUA_trials,1)
    
    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('No sink found for penetration ', num2str(penetration)))
        continue
    end
    
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    
    if any(v1Ch < 1) || any(v1Ch > 32)
        warning(strcat('Skipping penetration due to incomplete channel range: ', num2str(penetration)))
        continue
    end
    
    for ch_idx = 1:length(v1Ch)
        % Dioptic (condition 1) processing
        numTrials = size(MUA_trials{penetration,1}{1,1}, 1);
        MUAflashOut_diop = nan(2001, numTrials);
        for trl = 1:numTrials
            MUAflashOut_diop(:, trl) = MUA_trials{penetration,1}{1,1}{trl,1}(:, v1Ch(ch_idx));
        end
        
        avg_diop = median(MUAflashOut_diop, 2, 'omitnan');
        bl_diop = median(avg_diop(1:200));
        avg_diop_percent = 100 * ((avg_diop - bl_diop) / bl_diop);
        
        averageMUAMatrix_Diop(:, ch_idx, penetration) = avg_diop_percent;
        
        % Dichoptic (condition 3) processing
        numTrials = size(MUA_trials{penetration,1}{3,1}, 1);
        MUAflashOut_dichop = nan(2001, numTrials);
        for trl = 1:numTrials
            MUAflashOut_dichop(:, trl) = MUA_trials{penetration,1}{3,1}{trl,1}(:, v1Ch(ch_idx));
        end
        
        avg_dichop = median(MUAflashOut_dichop, 2, 'omitnan');
        bl_dichop = median(avg_dichop(1:200));
        avg_dichop_percent = 100 * ((avg_dichop - bl_dichop) / bl_dichop);
        
        averageMUAMatrix_Dichop(:, ch_idx, penetration) = avg_dichop_percent;
    end
end

%% Calculate grand averages
psAvg = mean(averageMUAMatrix_Diop, 3, 'omitnan');
nsAvg = mean(averageMUAMatrix_Dichop, 3, 'omitnan');

%% Plot each depth channel stacked on top of each other
figure;
hold on;
for ch_idx = 1:15
    offset = (ch_idx - 1) * 50; % Offset each channel to separate visually
    plot(tm_full, psAvg(:, ch_idx) + offset, 'b', 'LineWidth', 1.5);
    plot(tm_full, nsAvg(:, ch_idx) + offset, 'r', 'LineWidth', 1.5);
end
xlabel('Time (ms)');
ylabel('Depth Channels');
title('MUA Grand Average by Laminar Depth');

%% Save figure
answer = questdlg('Would you like to save this figure?', 'Save Figures', 'Yes', 'No', 'No');
switch answer
    case 'Yes'
        disp('Saving figure to plotDir...');
        cd(plotDir);
        figName = strcat('MUA_fig1e_laminarDepth_', '.png');
        saveas(gcf, figName);
        figName = strcat('MUA_fig1e_laminarDepth_', '.svg');
        saveas(gcf, figName);
    case 'No'
        disp('See plotDir for last saved versions.');
end
