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
close all
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';

cd(outputDir)
fileName = 'sortedData_211008_B_bmcBRFS001.mat';
load(fileName)

% Variables
probeLength = size(IDX(1).LFP_bb{1,1},2);
sdftm = STIM(1).sdftm;
timeLength = length(sdftm);
baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
xAxisTime = sdftm;
yAxisChannels = 1:32;

for chIdx = yAxisChannels
    %% Continuous line plots for individual units
    
    % Figure settings
    set(0,'DefaultFigureWindowStyle','normal')
    f = figure;
    set(f,'Position',[-2.0583e+03 7.6667 1186 1.3447e+03])
    
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
            trlLength = size(IDX(conditions(j)).CSD_gamma,1);
            for trl = 1:trlLength
                counter = counter + 1;
                muaAllTrls.(fields{i})(:,counter) = IDX(conditions(j)).CSD_gamma{trl,1}(:,chIdx);
            end
        end
        muaConditionmean = mean(muaAllTrls.(fields{i}),2);
        mua_baselinemean = mean(muaConditionmean(baselineTimeIndex,1));
        mua_blSubAvg.(fields{i}) = muaConditionmean - mua_baselinemean;
        minVals(i) = min(mua_blSubAvg.(fields{i}));
    end
    [sortedResponse,sortedResponseIdx] = sort(minVals,'ascend'); 
    [minVal,minValLoc] = min(minVals); % Location: 1="PO_LE", 2="PO_RE",etc.

    fields = fieldnames(cond);
    preferedStimulus = fields{minValLoc}(7:end);
    
    for i = 1: length(sortedResponseIdx)
        subplot(4,2,i)
        plot(sdftm,smoothdata(mua_blSubAvg.(fields{sortedResponseIdx(i)}),"gaussian",20))
    ylim([minVal -minVal])
        vline(0)
        xlim([sdftm(1) sdftm(end)])
        title(fields{sortedResponseIdx(i)},'interpreter','none')
        ylabel('uV')
    end
    
    titleText = {strcat('Channel Index = ',string(chIdx));strcat('Preferred Response = ',preferedStimulus)};
    sgtitle(titleText,'interpreter','none')
    
    

    %% Physical Alternation
    % Chose condition based on ocular preferences
    % cond.monoc_PO_LE = [8 10 16 18];
    % cond.monoc_PO_RE = [5 11 13 19];
    % cond.monoc_NOP_LE = [6 12 14 20];
    % cond.monoc_NPO_RE = [7 9 15 17];
    if minValLoc == 1 % PO LE
        nullToPreferred = 17;
        preferredToNull = 18;
    elseif minValLoc == 2 % PO RE
        nullToPreferred = 20;
        preferredToNull = 19;
    elseif minValLoc == 3 % NPO_LE
        nullToPreferred = 19;
        preferredToNull = 20;
    elseif minValLoc == 4 % NPO_RE
        nullToPreferred = 18;
        preferredToNull = 17;
    end
    
    % Average data
    % Null to preferred
    trlLength = size(IDX(nullToPreferred).CSD_gamma,1);
    counter = 0;
    for trl = 1:trlLength
        counter = counter + 1;
        muaAllTrls_NtP_adapter(:,counter) = IDX(nullToPreferred).CSD_gamma{trl,1}(:,chIdx); %first 800
        muaAllTrls_NtP_flash(:,counter) = IDX(nullToPreferred).CSD_gamma{trl,2}(:,chIdx); % second 800
    end
    muaConditionmean_NtP(:,1) = mean(muaAllTrls_NtP_adapter,2);
    muaConditionmean_NtP(:,2) = mean(muaAllTrls_NtP_flash,2);
    mua_baselinemean_NtP = mean(muaConditionmean_NtP(baselineTimeIndex,:));
    mua_blSubAvg_NtP = muaConditionmean_NtP - mua_baselinemean_NtP;
    
    % preferredToNull
    trlLength = size(IDX(preferredToNull).CSD_gamma,1);
    counter = 0;
    for trl = 1:trlLength
        counter = counter + 1;
        muaAllTrls_PtN(:,counter,1) = IDX(preferredToNull).CSD_gamma{trl,1}(:,chIdx); %first 800
        muaAllTrls_PtN(:,counter,2) = IDX(preferredToNull).CSD_gamma{trl,2}(:,chIdx); %second 800
    end
    muaConditionmean_PtN = squeeze(mean(muaAllTrls_PtN,2));
    mua_baselinemean_PtN = mean(muaConditionmean_PtN(baselineTimeIndex,:));
    mua_blSubAvg_PtN = muaConditionmean_PtN - mua_baselinemean_PtN(1);
    
    % Plot Data
    smoothed_NtP = smoothdata(mua_blSubAvg_NtP,"gaussian",20);
    smoothed_PtN = smoothdata(mua_blSubAvg_PtN,"gaussian",20);
    subplot(4,2,5)
    plot(sdftm,smoothed_NtP(:,1))
    hold on
    plot(sdftm,smoothed_PtN(:,1))
ylim([minVal -minVal])
    vline(0)
    legend ('Null adapter','Preferred Adapter')
    title('Physical Alternation, monocular adaptation')
    ylabel('uV')
    
    subplot(4,2,6)
    plot(sdftm,smoothed_NtP(:,2))
    hold on
    plot(sdftm,smoothed_PtN(:,2))
ylim([minVal -minVal])
    vline(0)
    legend ('Preferred flash (after null adapter)','Null flash (after preferred adapter)')
    title('Physical Alternation, monocular flash')
    
    %% Same Stim, different history, preferred vs null (BRFS)
    % Chose condition based on ocular preferences
    % cond.monoc_PO_LE = [8 10 16 18];
    % cond.monoc_PO_RE = [5 11 13 19];
    % cond.monoc_NOP_LE = [6 12 14 20];
    % cond.monoc_NPO_RE = [7 9 15 17];
    if minValLoc == 1 % PO LE
        nullToPreferred = 9;
        preferredToNull = 10;
    elseif minValLoc == 2 % PO RE
        nullToPreferred = 12;
        preferredToNull = 11;
    elseif minValLoc == 3 % NPO_LE
        nullToPreferred = 11;
        preferredToNull = 12;
    elseif minValLoc == 4 % NPO_RE
        nullToPreferred = 10;
        preferredToNull = 9;
    end
    
    % Average data
    % Null to preferred
    trlLength = size(IDX(nullToPreferred).CSD_gamma,1);
    counter = 0;
    for trl = 1:trlLength
        counter = counter + 1;
        muaAllTrls_NtP_adapter(:,counter) = IDX(nullToPreferred).CSD_gamma{trl,1}(:,chIdx); %first 800
        muaAllTrls_NtP_flash(:,counter) = IDX(nullToPreferred).CSD_gamma{trl,2}(:,chIdx); % second 800
    end
    muaConditionmean_NtP(:,1) = mean(muaAllTrls_NtP_adapter,2);
    muaConditionmean_NtP(:,2) = mean(muaAllTrls_NtP_flash,2);
    mua_baselinemean_NtP = mean(muaConditionmean_NtP(baselineTimeIndex,:));
    mua_blSubAvg_NtP = muaConditionmean_NtP - mua_baselinemean_NtP;
    
    % preferredToNull
    trlLength = size(IDX(preferredToNull).CSD_gamma,1);
    counter = 0;
    for trl = 1:trlLength
        counter = counter + 1;
        muaAllTrls_PtN(:,counter,1) = IDX(preferredToNull).CSD_gamma{trl,1}(:,chIdx); %first 800
        muaAllTrls_PtN(:,counter,2) = IDX(preferredToNull).CSD_gamma{trl,2}(:,chIdx); %second 800
    end
    muaConditionmean_PtN = squeeze(mean(muaAllTrls_PtN,2));
    mua_baselinemean_PtN = mean(muaConditionmean_PtN(baselineTimeIndex,:));
    mua_blSubAvg_PtN = muaConditionmean_PtN - mua_baselinemean_PtN(1);
    
    % Plot Data
    smoothed_NtP = smoothdata(mua_blSubAvg_NtP,"gaussian",20);
    smoothed_PtN = smoothdata(mua_blSubAvg_PtN,"gaussian",20);
    subplot(4,2,7)
    plot(sdftm,smoothed_NtP(:,1))
    hold on
    plot(sdftm,smoothed_PtN(:,1))
ylim([minVal -minVal])
    vline(0)
    legend ('Null adapter','Preferred Adapter')
    title('BRFS, monocular adaptation')
    ylabel('uV')
    xlabel('Time (ms) from monocular onset')

    subplot(4,2,8)
    plot(sdftm,smoothed_NtP(:,2))
    hold on
    plot(sdftm,smoothed_PtN(:,2))
ylim([minVal -minVal])
    vline(0)
    legend ('Preferred flash (after null adapter)','Null flash (after preferred adapter)')
    title('BRFS, binocular dichoptic flash')
    xlabel('Time (ms) from flash onset')

    %% Save output
    cd(outputDir)
    pdfOutputName = strcat(fileName(12:end-4),'_individualUnitClassification_CSD.pdf');
    exportgraphics(f,pdfOutputName,"Append",true)

end