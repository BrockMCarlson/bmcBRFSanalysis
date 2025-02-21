%% fig1d_DiopDichopAvg.m
% MUA grand average analysis with dioptic vs. dichoptic comparison
% Using a paired one-tailed test (dioptic > dichoptic) in 50 ms bins
% plus one a priori test in the final 500 ms of stimulation (no multiple-comparison correction).

%% Setup
disp('start time')
datetime
close
clearvars -except MUA_trials
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)
cd(dataDir)

if ~exist('MUA_trials','var')
    tic
    load('MUA_trials.mat') % format is MUA_trials{penetration,1}{cond,1}{trial,flash}
    toc
end

%% Preallocate
averageMUAMatrix_BRFSps = nan(2001,15,31);
averageMUAMatrix_BRFSns = nan(2001,15,31);

%% Loop through each penetration
for penetration = 1:size(MUA_trials,1)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');

    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to bottom of layer 4c
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

    %% Obtain Monocular preference via ANOVA (for potential classification)
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
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
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
        % % % if p < .05
            tuned(i,1) = true;
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
            % % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
            % % % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
            % % % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
            % % % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
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
        % % % else
        % % %     tuned(i,1) = false;
        % % %     prefMonoc(i,1) = NaN;
        % % %     nullMonoc(i,1) = NaN;
        % % %     prefCondOnFlash(i,1) = NaN;
        % % %     nullCondOnFlash(i,1) = NaN;
        % % % end

    end
    %% MUA in percent change for the BRFS conditions
    % Condition 1 (index=1): 'ps' (pref-stim?)
    % Condition 3 (index=3): 'ns' (null-stim?), i.e., dichoptic condition
    for i = 1:length(v1Ch)
        if ~tuned(i,1)
            continue
        end
        % Dioptic/Pref
        numTrials_ps = size(MUA_trials{penetration,1}{1,1},1);
        MUAflashOut_ps = nan(2001,numTrials_ps);
        for trl = 1:numTrials_ps
            MUAflashOut_ps(:,trl) = MUA_trials{penetration,1}{1,1}{trl,1}(:,v1Ch(i));
        end
        ps_avg = median(MUAflashOut_ps,2,"omitmissing");    
        psBl = median(ps_avg(1:200,:));
        if psBl == 0, psBl = 0.1; end
        ps_PercentC = 100*((ps_avg-psBl)./psBl);
        averageMUAMatrix_BRFSps(:,i,penetration) = ps_PercentC;

        % Dichoptic/Null
        numTrials_ns = size(MUA_trials{penetration,1}{3,1},1);
        MUAflashOut_ns = nan(2001,numTrials_ns);
        for trl = 1:numTrials_ns
            MUAflashOut_ns(:,trl) = MUA_trials{penetration,1}{3,1}{trl,1}(:,v1Ch(i));
        end
        ns_avg = median(MUAflashOut_ns,2,"omitmissing");
        nsBl = median(ns_avg(1:200,:));
        if nsBl == 0, nsBl = 0.1; end
        ns_PercentC = 100*((ns_avg-nsBl)./nsBl);  
        averageMUAMatrix_BRFSns(:,i,penetration) = ns_PercentC;
    end
    disp(strcat('Done with file number: ',string(penetration)))
end

%% Reshape across penetrations (405 contacts total if all included)
useIdx = squeeze(~isnan(averageMUAMatrix_BRFSps(1,1,:)));
ps_reshaped = reshape(averageMUAMatrix_BRFSps(:,:,useIdx),[2001,405]);
ns_reshaped = reshape(averageMUAMatrix_BRFSns(:,:,useIdx),[2001,405]);

%% Mean & SEM across all 405 "units" (or channels)
psAvg = smoothdata(mean(ps_reshaped,2,"omitmissing"),"gaussian",20);
nsAvg = smoothdata(mean(ns_reshaped,2,"omitmissing"),"gaussian",20);

contactNum = size(ps_reshaped,2);
psSEM = std(ps_reshaped,0,2,"omitmissing")./sqrt(contactNum);
nsSEM = std(ns_reshaped,0,2,"omitmissing")./sqrt(contactNum);

ps_S_avgPlusSEM  = psAvg + psSEM;
ps_S_avgMinusSEM = psAvg - psSEM;
ns_S_avgPlusSEM  = nsAvg + nsSEM;
ns_S_avgMinusSEM = nsAvg - nsSEM;

% Time vector
tm_full = (-200:1800)';  % Adjust if needed

%% Figure generation
lamCom = figure('Name','fig1d_DiopDichopAvg','Color','w');
hold on; box off;

% Shaded region for dioptic
fill([tm_full; flipud(tm_full)], ...
    [ps_S_avgPlusSEM; flipud(ps_S_avgMinusSEM)], ...
    [0.2, 0.4, 1], 'FaceAlpha', 0.4, 'EdgeColor', [0 0 0.5]);

% Shaded region for dichoptic
fill([tm_full; flipud(tm_full)], ...
    [ns_S_avgPlusSEM; flipud(ns_S_avgMinusSEM)], ...
    [1, 0.4, 0.4], 'FaceAlpha', 0.4, 'EdgeColor', [0.5 0 0]);

% Mean lines
plot(tm_full, psAvg, 'color', [0, 0, 1], 'LineWidth', 1.5); 
plot(tm_full, nsAvg, 'color', [1, 0, 0], 'LineWidth', 1.5);

ylim([0 45])
xlim([-200 1800])
plot([0 0], ylim, 'k', 'LineWidth', 2);
plot([1633 1633], ylim, 'k', 'LineWidth', 2);
plot(xlim, [0 0], 'k', 'LineWidth', 2);

title('Grand Average: Dioptic vs. Dichoptic')
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('% change from baseline', 'FontSize', 14, 'FontWeight', 'bold')
legend({'Dioptic','Dichoptic'}, 'Location','best')
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

%% --------------- NEW SECTION FOR SIGNIFICANCE (PAIRED, ONE-TAILED) ---------------

% 1) Bin-based test from 0 to 1600 ms in 50 ms intervals, Bonferroni corrected
bin_width = 100; 
time_bins = 0:bin_width:1600; % Bins from 0 to 1600 in steps of 50
num_bins  = length(time_bins) - 1;

alpha_original = 0.05;
alpha_corrected = alpha_original / num_bins;

pvals = nan(num_bins,1);
sigBins = false(num_bins,1);

% For plotting significance asterisks, pick a Y position above the max of the dioptic trace
maxVal_all = max([psAvg; nsAvg]);
y_pos = maxVal_all - 5;  % a bit below the highest average

for i = 1:num_bins
    % Indices for this 50 ms window
    bin_idx = tm_full >= time_bins(i) & tm_full < time_bins(i+1);
    
    % Average within this bin for each neuron
    diop_bin   = mean(ps_reshaped(bin_idx, :), 1, 'omitnan');  % (1 x 405)
    dichop_bin = mean(ns_reshaped(bin_idx, :), 1, 'omitnan');
    
    % Paired one-tailed test, dioptic > dichoptic
    [~, p] = ttest(diop_bin, dichop_bin, 'Tail','right');
    pvals(i) = p;
    
    if p < alpha_corrected
        sigBins(i) = true;
        % Mark significance at the bin's midpoint
        x_pos = mean(time_bins(i:i+1));
        text(x_pos, y_pos, '*', 'FontSize',14, 'HorizontalAlignment','center');
    end
end

%% Save figure
answer = questdlg('Would you like to save these figures?', ...
    'Save Figures', 'Yes', 'No', 'No');
switch answer
    case 'Yes'
        disp('Saving figures to plotDir...')
        cd(plotDir)
        figName_lamCom_png = strcat('fig1d','_DiopticVsDichoptic_','.png');
        saveas(lamCom, figName_lamCom_png)
        figName_lamCom_svg = strcat('fig1d','_DiopticVsDichoptic_','.svg');
        saveas(lamCom, figName_lamCom_svg)
    case 'No'
        disp(plotDir)
        disp(codeDir)
        disp('Please see plotDir for last saved versions if needed.')
end


%% 2) A priori single test (final 500 ms)

startWindow = 1200;
endWindow   = 1800;
win_idx = tm_full >= startWindow & tm_full < endWindow;

diop_window   = mean(ps_reshaped(win_idx, :), 1, 'omitnan');
dichop_window = mean(ns_reshaped(win_idx, :), 1, 'omitnan');

% One-tailed paired test for dioptic > dichoptic
[H, p_500ms, CI, stats_500ms]  = ttest(diop_window, dichop_window, 'Tail','right');

fprintf('\n==== A Priori 500 ms Window Test ====\n');
fprintf('Paired one-tailed t-test (dioptic > dichoptic) in last 500 ms:\n');
fprintf('p = %.5f (no multiple-comparison correction applied)\n\n', p_500ms);


% Gather your data
cond1Mean = mean(diop_window);      % dioptic
cond1SD   = std(diop_window);
cond2Mean = mean(dichop_window);    % dichoptic
cond2SD   = std(dichop_window);

meanDiff  = cond1Mean - cond2Mean;
n         = length(diop_window);
SE        = std(diop_window - dichop_window)/sqrt(n);

tVal      = stats_500ms.tstat;
dfVal     = stats_500ms.df;
pVal      = p_500ms;

% --- CREATE A TABLE ---
contrastName = "Dioptic vs Dichoptic";
resultsTable = table(contrastName, ...
                     cond1Mean, cond1SD, ...
                     cond2Mean, cond2SD, ...
                     meanDiff, SE, tVal, dfVal, pVal, ...
    'VariableNames', ...
    {'Contrast','MeanCond1','SDCond1','MeanCond2','SDCond2','MeanDiff','SE','tVal','df','pValue'});

% --- DISPLAY TABLE IN COMMAND WINDOW ---
disp(resultsTable);
