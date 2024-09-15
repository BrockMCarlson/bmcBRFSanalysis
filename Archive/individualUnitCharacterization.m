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



%% Create dir loop
cd(outputDir)
sortedDataFiles = dir('*sortedData*');
for dataFile = 1:length(sortedDataFiles)
    close all
    fileToLoad = sortedDataFiles(dataFile).name;

    load(fileToLoad)
    
    % Variables
    probeLength = size(IDX(1).LFP_bb{1,1},2);
    sdftm = STIM(1).sdftm;
    timeLength = length(sdftm);
    baselineTimeIndex = find(sdftm == -100):find(sdftm == 0); 
    xAxisTime = sdftm;
    if dataFile <=21
        yAxisChannels = 1:32;
         yAxisLim =      [1 32];
    elseif dataFile >=22
        yAxisChannels = 1:64;
         yAxisLim =      [1 64];
    end

    for chIdx = yAxisChannels
        %% Continuous line plots for individual units
        
        % Figure settings
        set(0,'DefaultFigureWindowStyle','normal')
        f = figure;
        set(f,'Position',[1000 74.3333 1.1197e+03 1.1633e+03])
        
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
            muaConditionmean = mean(muaAllTrls.(fields{i}),2);
            mua_baselinemean = mean(muaConditionmean(baselineTimeIndex,1));
            mua_blSubAvg.(fields{i}) = muaConditionmean - mua_baselinemean;
            maxVals(i) = max(mua_blSubAvg.(fields{i}));
            minVals(i) = min(mua_blSubAvg.(fields{i}));
        end
        [sortedResponse,sortedResponseIdx] = sort(maxVals,'descend'); 
        [maxVal,maxValLoc] = max(maxVals); % Location: 1="PO_LE", 2="PO_RE",etc.
        [minVal,inValLoc] = min(maxVals); % Location: 1="PO_LE", 2="PO_RE",etc.
    
        fields = fieldnames(cond);
        preferedStimulus = fields{maxValLoc}(7:end);
        
        for i = 1: length(sortedResponseIdx)
            subplot(4,2,i)
            plot(sdftm,smoothdata(mua_blSubAvg.(fields{sortedResponseIdx(i)}),"gaussian",10))
            ylim([-1 maxVal])
            vline(0)
            xlim([sdftm(1) sdftm(end)])
            title(fields{sortedResponseIdx(i)},'interpreter','none')
            ylabel('uV')
        end
        
        titleText = {strcat('Channel Index = ',string(chIdx));strcat('Preferred Response = ',fields{maxValLoc}(7:end))};
        sgtitle(titleText,'interpreter','none')
        
        
        %% Dioptic vs dichoptic
        % % % diopticCond = 1;
        % % % trlLength = size(IDX(diopticCond).MUAe,1);
        % % % counter = 0;
        % % % for trl = 1:trlLength
        % % %     counter = counter + 1;
        % % %     muaAllTrls_dioptic(:,counter) = IDX(diopticCond).MUAe{trl,1}(:,chIdx);
        % % % end
        % % % muaConditionmean_dioptic = mean(muaAllTrls_dioptic,2);
        % % % mua_baselinemean_dioptic = mean(muaConditionmean_dioptic(baselineTimeIndex,1));
        % % % mua_blSubAvg_dioptic = muaConditionmean_dioptic - mua_baselinemean_dioptic;
        % % % 
        % % % if maxValLoc == 1 %Left Eye preferred
        % % %     dichopticCond = 3;
        % % % elseif maxValLoc == 2 % Right Eye Preferred
        % % %     dichopticCond = 4;
        % % % else
        % % %     warning('unit preferres "null" orientation')
        % % % end
        % % % trlLength = size(IDX(dichopticCond).MUAe,1);
        % % % counter = 0;
        % % % for trl = 1:trlLength
        % % %     counter = counter + 1;
        % % %     muaAllTrls_dichoptic(:,counter) = IDX(dichopticCond).MUAe{trl,1}(:,chIdx);
        % % % end
        % % % muaConditionmean_dichoptic = mean(muaAllTrls_dichoptic,2);
        % % % mua_baselinemean_dichoptic = mean(muaConditionmean_dichoptic(baselineTimeIndex,1));
        % % % mua_blSubAvg_dichoptic = muaConditionmean_dichoptic - mua_baselinemean_dichoptic;
        % % % 
        % % % %Plot
        % % % smoothed_dioptic = smoothdata(mua_blSubAvg_dioptic,"gaussian",10);
        % % % smoothed_dichoptic = smoothdata(mua_blSubAvg_dichoptic,"gaussian",10);
        % % % 
        % % % figure
        % % % plot(sdftm,smoothed_dioptic)
        % % % hold on
        % % % plot(sdftm,smoothed_dichoptic)
        
        %% Physical Alternation
        % Chose condition based on ocular preferences
        % cond.monoc_PO_LE = [8 10 16 18];
        % cond.monoc_PO_RE = [5 11 13 19];
        % cond.monoc_NOP_LE = [6 12 14 20];
        % cond.monoc_NPO_RE = [7 9 15 17];
        if maxValLoc == 1 % PO LE
            nullToPreferred = 17;
            preferredToNull = 18;
        elseif maxValLoc == 2 % PO RE
            nullToPreferred = 20;
            preferredToNull = 19;
        elseif maxValLoc == 3 % NPO_LE
            nullToPreferred = 19;
            preferredToNull = 20;
        elseif maxValLoc == 4 % NPO_RE
            nullToPreferred = 18;
            preferredToNull = 17;
        end
        
        % Average data
        % Null to preferred
        trlLength = size(IDX(nullToPreferred).MUAe,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_NtP_adapter(:,counter) = IDX(nullToPreferred).MUAe{trl,1}(:,chIdx); %first 800
            muaAllTrls_NtP_flash(:,counter) = IDX(nullToPreferred).MUAe{trl,2}(:,chIdx); % second 800
        end
        muaConditionmean_NtP(:,1) = mean(muaAllTrls_NtP_adapter,2);
        muaConditionmean_NtP(:,2) = mean(muaAllTrls_NtP_flash,2);
        mua_baselinemean_NtP = mean(muaConditionmean_NtP(baselineTimeIndex,:));
        mua_blSubAvg_NtP = muaConditionmean_NtP - mua_baselinemean_NtP;
        maxValBRFS = max(mua_blSubAvg_NtP(:,2));
        
        % preferredToNull
        trlLength = size(IDX(preferredToNull).MUAe,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_PtN(:,counter,1) = IDX(preferredToNull).MUAe{trl,1}(:,chIdx); %first 800
            muaAllTrls_PtN(:,counter,2) = IDX(preferredToNull).MUAe{trl,2}(:,chIdx); %second 800
        end
        muaConditionmean_PtN = squeeze(mean(muaAllTrls_PtN,2));
        mua_baselinemean_PtN = mean(muaConditionmean_PtN(baselineTimeIndex,:));
        mua_blSubAvg_PtN = muaConditionmean_PtN - mua_baselinemean_PtN(1);

        
        % Plot Data
        smoothed_NtP = smoothdata(mua_blSubAvg_NtP,"gaussian",10);
        smoothed_PtN = smoothdata(mua_blSubAvg_PtN,"gaussian",10);
        subplot(4,2,5)
        plot(sdftm,smoothed_NtP(:,1))
        hold on
        plot(sdftm,smoothed_PtN(:,1))
        ylim([-1 maxValBRFS+.1])
        vline(0)
        legend ('Null adapter','Preferred Adapter')
        title('Physical Alternation, monocular adaptation')
        ylabel('uV')
        
        subplot(4,2,6)
        plot(sdftm,smoothed_NtP(:,2))
        hold on
        plot(sdftm,smoothed_PtN(:,2))
        ylim([-1 maxValBRFS+.1])
        vline(0)
        legend ('Preferred flash (after null adapter)','Null flash (after preferred adapter)')
        title('Physical Alternation, monocular flash')
        
        %% Same Stim, different history, preferred vs null (BRFS)
        % Chose condition based on ocular preferences
        % cond.monoc_PO_LE = [8 10 16 18];
        % cond.monoc_PO_RE = [5 11 13 19];
        % cond.monoc_NOP_LE = [6 12 14 20];
        % cond.monoc_NPO_RE = [7 9 15 17];
        if maxValLoc == 1 % PO LE
            nullToPreferred = 9;
            preferredToNull = 10;
        elseif maxValLoc == 2 % PO RE
            nullToPreferred = 12;
            preferredToNull = 11;
        elseif maxValLoc == 3 % NPO_LE
            nullToPreferred = 11;
            preferredToNull = 12;
        elseif maxValLoc == 4 % NPO_RE
            nullToPreferred = 10;
            preferredToNull = 9;
        end
        
        % Average data
        % Null to preferred
        trlLength = size(IDX(nullToPreferred).MUAe,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_NtP_adapter(:,counter) = IDX(nullToPreferred).MUAe{trl,1}(:,chIdx); %first 800
            muaAllTrls_NtP_flash(:,counter) = IDX(nullToPreferred).MUAe{trl,2}(:,chIdx); % second 800
        end
        muaConditionmean_NtP(:,1) = mean(muaAllTrls_NtP_adapter,2);
        muaConditionmean_NtP(:,2) = mean(muaAllTrls_NtP_flash,2);
        mua_baselinemean_NtP = mean(muaConditionmean_NtP(baselineTimeIndex,:));
        mua_blSubAvg_NtP = muaConditionmean_NtP - mua_baselinemean_NtP;
        
        % preferredToNull
        trlLength = size(IDX(preferredToNull).MUAe,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_PtN(:,counter,1) = IDX(preferredToNull).MUAe{trl,1}(:,chIdx); %first 800
            muaAllTrls_PtN(:,counter,2) = IDX(preferredToNull).MUAe{trl,2}(:,chIdx); %second 800
        end
        muaConditionmean_PtN = squeeze(mean(muaAllTrls_PtN,2));
        mua_baselinemean_PtN = mean(muaConditionmean_PtN(baselineTimeIndex,:));
        mua_blSubAvg_PtN = muaConditionmean_PtN - mua_baselinemean_PtN(1);
        
        % Plot Data
        smoothed_NtP = smoothdata(mua_blSubAvg_NtP,"gaussian",10);
        smoothed_PtN = smoothdata(mua_blSubAvg_PtN,"gaussian",10);
        subplot(4,2,7)
        plot(sdftm,smoothed_NtP(:,1))
        hold on
        plot(sdftm,smoothed_PtN(:,1))
        ylim([-1 maxValBRFS+.1])
        vline(0)
        legend ('Null adapter','Preferred Adapter')
        title('BRFS, monocular adaptation')
        ylabel('uV')
        xlabel('Time (ms) from monocular onset')
    
        subplot(4,2,8)
        plot(sdftm,smoothed_NtP(:,2))
        hold on
        plot(sdftm,smoothed_PtN(:,2))
        ylim([-1 maxValBRFS+.1])
        vline(0)
        legend ('Preferred flash (after null adapter)','Null flash (after preferred adapter)')
        title('BRFS, binocular dichoptic flash')
        xlabel('Time (ms) from flash onset')
    
        %% Save output
        cd(outputDir)
        pdfOutputName = strcat('individualUnitClassification_MUAe_',fileToLoad(12:end-4),'_.pdf');
        if chIdx == 1
            exportgraphics(f,pdfOutputName)
        else
            exportgraphics(f,pdfOutputName,"Append",true)
        end
    
    end
end

% Hello World