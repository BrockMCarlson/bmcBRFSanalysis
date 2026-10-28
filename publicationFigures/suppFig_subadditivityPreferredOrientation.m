%% figS_subadditivityPreferredOrientation.m
% Supplementary analysis: preferred-orientation subadditivity in MUA.
%
% Compares:
%   1) Monocular preferred-orientation response in the preferred eye
%   2) Monocular preferred-orientation response in the null eye
%   3) Sum of those two monocular responses
%   4) Dioptic binocular response (MUA_trials condition 1)
%
% Responses are expressed as percent change from baseline, following the
% existing fig1b_monocPrefGrandAvg.m and fig1d_DiopDichopAvg.m scripts.
%
% Primary test: paired, one-tailed t-test on each contact's mean response
% over samples 488:1000 (final 512 ms by elapsed time).
%
% Secondary test: paired, one-tailed t-tests in 100-ms bins from 0-800 ms,
% with Bonferroni correction across the eight bins.
%
% NOTE: As in the existing scripts, contacts are the observations for the
% t-tests. Contacts from the same penetration may not be independent.
%
% Data format: MUA_trials{penetration,1}{condition,1}{trial,flash}

%% Setup
disp('Starting preferred-orientation subadditivity analysis')
close all

% Keep an existing MUA_trials variable if one is already in the workspace.
if ~exist('MUA_trials','var')
    workingPC = 'home'; % options: 'home', 'office'

    if strcmp(workingPC,'home')
        codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
        dataDir = 'S:\TrialTriggeredLFPandMUA';
        plotDir = 'C:\Users\Brock Carlson\Documents\Coherence manuscript - home pc\officialFigureOutputsFromMATLAB';
        laminarFile = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx';
    elseif strcmp(workingPC,'office')
        codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
        dataDir = 'D:\TrialTriggeredLFPandMUA';
        plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
        laminarFile = 'C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx';
    end

    addpath(codeDir)
    officLamAssign = importLaminarAssignments(laminarFile, "AnalysisList", [2, Inf]);
    cd(dataDir)
    load('MUA_trials.mat')
else
    % MUA_trials already exists. Set these paths for loading assignments
    % and saving figures if they have not been set in the workspace.
    workingPC = 'home'; % options: 'home', 'office'
    if strcmp(workingPC,'home')
        codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
        dataDir = 'S:\TrialTriggeredLFPandMUA';
        plotDir = 'C:\Users\Brock Carlson\Documents\Coherence manuscript - home pc\officialFigureOutputsFromMATLAB';
        laminarFile = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx';
    else
        codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
        dataDir = 'D:\TrialTriggeredLFPandMUA';
        plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
        laminarFile = 'C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx';
    end
    addpath(codeDir)
    if ~exist('officLamAssign','var')
        officLamAssign = importLaminarAssignments(laminarFile, "AnalysisList", [2, Inf]);
    end
end

% Time indices follow the existing scripts and your specified analysis.
tm_full = (-200:1800)';   % 2001 samples; sample 201 corresponds to 0 ms
tm_plot = tm_full(1:1000); % -200 to 799 ms: first 800-ms stimulus only
tm_bl = 1:200;
tm_1stOnset = 201:1000;
analysisWindow = 488:1000; % final 512 ms by elapsed time (513 samples inclusive)

% Monocular condition groups, copied from the existing monocular analysis.
monoc_1 = [5, 11, 13, 19]; % PO RightEye
monoc_2 = [8, 10, 16, 18]; % PO LeftEye
monoc_3 = [7, 9, 15, 17]; % NPO RightEye
monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
monocGroups = {monoc_1, monoc_2, monoc_3, monoc_4};

nPen = size(MUA_trials,1);
nChV1 = 15;

% Store contact-level time courses (time x V1 contact x penetration).
monoPrefEye = nan(2001,nChV1,nPen);
monoNullEye = nan(2001,nChV1,nPen);
diopticResp = nan(2001,nChV1,nPen);
prefMonocAll = nan(nChV1,nPen);
validPen = false(1,nPen);

%% Calculate responses
for penetration = 1:nPen
    fprintf('Processing penetration %d of %d\n', penetration, nPen);

    % Match the V1 contact selection used in the existing figure scripts.
    granBtm = officLamAssign.stFold4c(penetration);
    if isnan(granBtm)
        warning('No layer 4C bottom assignment for penetration %d; skipping.', penetration);
        continue
    end

    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;

    if strcmp(string(officLamAssign.ChToUse(penetration)),"1:32")
        if any(v1Ch > 32) || any(v1Ch < 1)
            warning('Skipping penetration %d because its V1 column is incomplete.', penetration);
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('Skipping penetration %d because its V1 column is incomplete.', penetration);
            continue
        end
    end

    % Gather monocular trials for each of the four eye/orientation groups.
    % Data dimensions are time x contact x trial.
    monocData = cell(1,4);
    groupBaseline = nan(4,nChV1);
    groupResponse = nan(4,nChV1);

    for g = 1:4
        trialMatrices = {};
        for cond = monocGroups{g}
            if cond > numel(MUA_trials{penetration,1}) || ...
                    isempty(MUA_trials{penetration,1}{cond,1})
                continue
            end
            nTrials = size(MUA_trials{penetration,1}{cond,1},1);
            for trl = 1:nTrials
                thisTrial = MUA_trials{penetration,1}{cond,1}{trl,1};
                if isempty(thisTrial) || size(thisTrial,1) < 2001 || size(thisTrial,2) < max(v1Ch)
                    continue
                end
                trialMatrices{end+1} = abs(thisTrial(:,v1Ch)); %#ok<SAGROW>
            end
        end

        if isempty(trialMatrices)
            monocData{g} = [];
            continue
        end

        allTrials = cat(3,trialMatrices{:});
        monocData{g} = allTrials;

        % Preference classification uses the same time range and baseline
        % procedure as the existing monocular figure script.
        trialMedianResponse = squeeze(median(allTrials(200:450,:,:),1,'omitnan'));
        if isvector(trialMedianResponse)
            trialMedianResponse = reshape(trialMedianResponse,size(allTrials,2),[]);
        end
        groupResponse(g,:) = median(trialMedianResponse,2,'omitnan')';

        trialMedianBaseline = squeeze(median(allTrials(tm_bl,:,:),1,'omitnan'));
        if isvector(trialMedianBaseline)
            trialMedianBaseline = reshape(trialMedianBaseline,size(allTrials,2),[]);
        end
        groupBaseline(g,:) = median(trialMedianBaseline,2,'omitnan');
    end

    % Skip if any monocular group is missing; preference cannot be assigned.
    if any(cellfun(@isempty,monocData))
        warning('Missing monocular trials for penetration %d; skipping it.', penetration);
        continue
    end

    % Existing preference assignment subtracts the median baseline across
    % the four groups, then chooses the group with the largest median response.
    for ch = 1:nChV1
        baselineForPreference = median(groupBaseline(:,ch),'omitnan');
        responseForPreference = groupResponse(:,ch) - baselineForPreference;
        if all(isnan(responseForPreference))
            continue
        end
        [~,prefGroup] = max(responseForPreference);
        prefMonocAll(ch,penetration) = prefGroup;

        % Preferred-orientation monocular groups in each eye:
        % group 1 = PO right eye; group 2 = PO left eye;
        % group 3 = NPO right eye; group 4 = NPO left eye.
        % The null eye receives the same orientation as the preferred eye.
        switch prefGroup
            case 1
                prefEyeGroup = 1; nullEyeGroup = 2;
            case 2
                prefEyeGroup = 2; nullEyeGroup = 1;
            case 3
                prefEyeGroup = 3; nullEyeGroup = 4;
            case 4
                prefEyeGroup = 4; nullEyeGroup = 3;
        end

        % Median across trials, then percent change from that condition's
        % own baseline, matching the existing monocular script.
        prefAvg = median(monocData{prefEyeGroup}(:,ch,:),3,'omitnan');
        prefBl = median(prefAvg(tm_bl),'omitnan');
        if prefBl == 0, prefBl = 0.1; end
        monoPrefEye(:,ch,penetration) = 100*((prefAvg-prefBl)./prefBl);

        nullAvg = median(monocData{nullEyeGroup}(:,ch,:),3,'omitnan');
        nullBl = median(nullAvg(tm_bl),'omitnan');
        if nullBl == 0, nullBl = 0.1; end
        monoNullEye(:,ch,penetration) = 100*((nullAvg-nullBl)./nullBl);
    end

    % Dioptic response: condition 1 is the dioptic preferred-orientation
    % condition in fig1d_DiopDichopAvg.m. Use its median across trials and
    % normalize to its own baseline, as in that script.
    if numel(MUA_trials{penetration,1}) >= 1 && ...
            ~isempty(MUA_trials{penetration,1}{1,1})
        nTrials = size(MUA_trials{penetration,1}{1,1},1);
        for ch = 1:nChV1
            trialData = nan(2001,nTrials);
            for trl = 1:nTrials
                thisTrial = MUA_trials{penetration,1}{1,1}{trl,1};
                if ~isempty(thisTrial) && size(thisTrial,1) >= 2001 && size(thisTrial,2) >= max(v1Ch)
                    trialData(:,trl) = thisTrial(:,v1Ch(ch));
                end
            end
            diopticAvg = median(trialData,2,'omitnan');
            diopticBl = median(diopticAvg(tm_bl),'omitnan');
            if diopticBl == 0, diopticBl = 0.1; end
            diopticResp(:,ch,penetration) = 100*((diopticAvg-diopticBl)./diopticBl);
        end
        validPen(penetration) = true;
    end
end

%% Flatten valid contacts across penetrations
% Exclude contacts without all three usable response traces.
prefFlat = reshape(monoPrefEye(1:1000,:,:),1000,[]);
nullFlat = reshape(monoNullEye(1:1000,:,:),1000,[]);
diopticFlat = reshape(diopticResp(1:1000,:,:),1000,[]);
validContacts = all(~isnan(prefFlat(1:200,:)),1) & ...
                all(~isnan(nullFlat(1:200,:)),1) & ...
                all(~isnan(diopticFlat(1:200,:)),1);
prefFlat = prefFlat(:,validContacts);
nullFlat = nullFlat(:,validContacts);
diopticFlat = diopticFlat(:,validContacts);

if isempty(prefFlat)
    error('No contacts had valid preferred-eye, null-eye, and dioptic responses.');
end

sumFlat = prefFlat + nullFlat;
nContacts = size(sumFlat,2);
timeMs = tm_plot;

% Grand-average traces and SEM across included contacts.
prefMean = mean(prefFlat,2,'omitnan');
nullMean = mean(nullFlat,2,'omitnan');
sumMean = mean(sumFlat,2,'omitnan');
diopticMean = mean(diopticFlat,2,'omitnan');

prefSEM = std(prefFlat,0,2,'omitnan') ./ sqrt(sum(~isnan(prefFlat),2));
nullSEM = std(nullFlat,0,2,'omitnan') ./ sqrt(sum(~isnan(nullFlat),2));
sumSEM = std(sumFlat,0,2,'omitnan') ./ sqrt(sum(~isnan(sumFlat),2));
diopticSEM = std(diopticFlat,0,2,'omitnan') ./ sqrt(sum(~isnan(diopticFlat),2));

%% Primary statistical test: mean response over samples 488:1000
% Convert the full-data indices to the truncated (1:1000) arrays.
primaryIdx = analysisWindow;
prefWindow = mean(prefFlat(primaryIdx,:),1,'omitnan');
nullWindow = mean(nullFlat(primaryIdx,:),1,'omitnan');
sumWindow = mean(sumFlat(primaryIdx,:),1,'omitnan');
diopticWindow = mean(diopticFlat(primaryIdx,:),1,'omitnan');

validPrimary = ~isnan(sumWindow) & ~isnan(diopticWindow);
sumWindow = sumWindow(validPrimary);
diopticWindow = diopticWindow(validPrimary);
differenceWindow = sumWindow - diopticWindow;

[~,pPrimary,ciPrimary,statsPrimary] = ttest(sumWindow,diopticWindow,'Tail','right');
fprintf('\nPRIMARY TEST: summed monocular > dioptic, samples %d:%d\n',analysisWindow(1),analysisWindow(end));
fprintf('Contacts included: %d\n',numel(differenceWindow));
fprintf('Mean summed monocular response: %.4f %% change from baseline\n',mean(sumWindow,'omitnan'));
fprintf('Mean dioptic response: %.4f %% change from baseline\n',mean(diopticWindow,'omitnan'));
fprintf('Mean paired difference (sum - dioptic): %.4f\n',mean(differenceWindow,'omitnan'));
fprintf('Paired t(%d) = %.4f, one-tailed p = %.6g; 95%% CI for difference = [%.4f, %.4f]\n', ...
    statsPrimary.df,statsPrimary.tstat,pPrimary,ciPrimary(1),ciPrimary(2));

%% Secondary time-binned tests with Bonferroni correction
runTimeBinnedTests = true;
binWidth = 100;
timeBinEdges = 0:binWidth:800; % stimulus-relative time in ms
nBins = numel(timeBinEdges)-1;
pBins = nan(1,nBins);
tBins = nan(1,nBins);
meanDiffBins = nan(1,nBins);
binCenters = (timeBinEdges(1:end-1)+timeBinEdges(2:end))/2;
bonferroniAlpha = 0.05/nBins;

if runTimeBinnedTests
    for b = 1:nBins
        % tm_plot is -200:799; each bin is defined relative to stimulus onset.
        idx = find(timeMs >= timeBinEdges(b) & timeMs < timeBinEdges(b+1));
        sumBin = mean(sumFlat(idx,:),1,'omitnan');
        diopticBin = mean(diopticFlat(idx,:),1,'omitnan');
        validBin = ~isnan(sumBin) & ~isnan(diopticBin);
        [~,pBins(b),~,statsBin] = ttest(sumBin(validBin),diopticBin(validBin),'Tail','right');
        tBins(b) = statsBin.tstat;
        meanDiffBins(b) = mean(sumBin(validBin)-diopticBin(validBin),'omitnan');
    end
    fprintf('\nTIME-BINNED TESTS (one-tailed; Bonferroni alpha = %.5f for %d bins)\n',bonferroniAlpha,nBins);
    disp(table(timeBinEdges(1:end-1)',timeBinEdges(2:end)',meanDiffBins',tBins',pBins', ...
        (pBins' < bonferroniAlpha), ...
        'VariableNames',{'Start_ms','End_ms','MeanSumMinusDioptic','tStatistic','pValue','BonferroniSignificant'}));
end

%% Plot supplementary figure
figS = figure('Color','w','Position',[100 100 1000 600]);
hold on

% SEM shading helper blocks are written explicitly for compatibility.
fill([timeMs; flipud(timeMs)], ...
    [sumMean-sumSEM; flipud(sumMean+sumSEM)], ...
    [0.75 0.75 0.75],'FaceAlpha',0.35,'EdgeColor','none');
fill([timeMs; flipud(timeMs)], ...
    [diopticMean-diopticSEM; flipud(diopticMean+diopticSEM)], ...
    [0.65 0.78 0.95],'FaceAlpha',0.35,'EdgeColor','none');

% Include the two monocular components as lighter/dashed lines so the sum
% is transparent, while emphasizing the summed prediction and dioptic data.
plot(timeMs,prefMean,'Color',[0.45 0.45 0.45],'LineStyle',':','LineWidth',1.2);
plot(timeMs,nullMean,'Color',[0.25 0.55 0.35],'LineStyle',':','LineWidth',1.2);
plot(timeMs,sumMean,'Color',[0.15 0.15 0.15],'LineWidth',2.2);
plot(timeMs,diopticMean,'Color',[0 0.25 0.85],'LineWidth',2.2);

xline(0,'k-');
xline(800,'k:');
yline(0,'k:');
xlim([-200 800]);
xlabel('Time from stimulus onset (ms)');
ylabel('% change from baseline');
title('Preferred-orientation subadditivity: summed monocular vs. dioptic MUA');
legend({'Summed monocular SEM','Dioptic SEM','Preferred-eye monocular','Null-eye monocular', ...
    'Summed monocular','Dioptic'},'Location','best','Interpreter','none');
set(gca,'FontSize',12,'LineWidth',1);
box off

% Mark bins with Bonferroni-significant evidence that sum > dioptic.
if runTimeBinnedTests
    yLimits = ylim;
    markerY = yLimits(2) - 0.05*(yLimits(2)-yLimits(1));
    for b = 1:nBins
        if pBins(b) < bonferroniAlpha
            text(binCenters(b),markerY,'*','FontSize',15,'FontWeight','bold', ...
                'HorizontalAlignment','center','Color','k');
        end
    end
end

%% Save outputs
if ~exist(plotDir,'dir')
    warning('Output directory does not exist: %s. Figure will remain open but not be saved.',plotDir);
else
    saveas(figS,fullfile(plotDir,'figS_Subadditivity_PreferredOrientation.png'));
    saveas(figS,fullfile(plotDir,'figS_Subadditivity_PreferredOrientation.svg'));
    save(fullfile(plotDir,'figS_Subadditivity_PreferredOrientation.mat'), ...
        'monoPrefEye','monoNullEye','diopticResp','prefMonocAll','validPen', ...
        'validContacts','prefFlat','nullFlat','diopticFlat','sumFlat','timeMs', ...
        'analysisWindow','sumWindow','diopticWindow','differenceWindow', ...
        'pPrimary','ciPrimary','statsPrimary','pBins','tBins','meanDiffBins', ...
        'timeBinEdges','bonferroniAlpha','nContacts');
    fprintf('\nSaved supplementary figure and analysis data to:\n%s\n',plotDir);
end

fprintf('\nAnalysis complete. Number of contacts included in grand average: %d\n',nContacts);

