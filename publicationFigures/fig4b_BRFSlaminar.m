%% fig3b
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur

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


%%
averageMUAMatrix_BRFSps = nan(2001,15,31);
averageLFPMatrix_BRFSns = nan(2001,15,31);
for penetration = 1:size(MUA_trials,1)
    
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
   

    % MUA in percent change
    % BRFS pref vs null
    % Get the number of trials for the chosen condition
    for i = 1:length(v1Ch)
        numTrials_ps = size(MUA_trials{penetration,1}{prefCondOnFlash(i,1),1},1);
        MUAflashOut_ps = nan(2001,numTrials_ps);
        for trl = 1:numTrials_ps
            MUAflashOut_ps(:,trl) = MUA_trials{penetration,1}{prefCondOnFlash(i,1),1}{trl,1}(:,v1Ch(i));
        end
        
        % median across trials
        ps_avg = median(MUAflashOut_ps,2,"omitmissing"); % input is (tm,trl)    
        % Calculate as Percent Change
        %              X(t) - avgBl
        % %Ch = 100 * -------------
        %                 avgBl
        psBl = median(ps_avg(1:200,:));
        if psBl == 0
            psBl = .1;
        end
        ps_PercentC = 100*((ps_avg-psBl)./psBl);
        averageMUAMatrix_BRFSps(:,i,penetration) = ps_PercentC;

    
        % Get the number of trials for the chosen condition
        numTrials_ns = size(MUA_trials{penetration,1}{nullCondOnFlash(i,1) ,1},1);
        MUAflashOut_ns = nan(2001,numTrials_ns);
        for trl = 1:numTrials_ns
            MUAflashOut_ns(:,trl) = MUA_trials{penetration,1}{nullCondOnFlash(i,1),1}{trl,1}(:,v1Ch(i));
        end
        % Average across trials and save output
        ns_avg = median(MUAflashOut_ns,2,"omitmissing"); % input is (tm,trl)    
        nsBl = median(ns_avg(1:200,:));
        if nsBl == 0
            nsBl = .1;
        end
        ns_PercentC = 100*((ns_avg-nsBl)./nsBl);  
        averageMUAMatrix_BRFSns(:,i,penetration) = ns_PercentC; % Average across trl. averageMUAMatrix is (tm x ch x x penetration)
    end
    disp(strcat('Done with file number: ',string(penetration)))
end



%% Organize into compartments, median and std across contacts + penetration
useIdx = squeeze(~isnan(averageMUAMatrix_BRFSps(1,1,:))); 
ps_S = reshape(averageMUAMatrix_BRFSps(:,1:5,useIdx), [2001, 135]);
ps_G = reshape(averageMUAMatrix_BRFSps(:,6:10,useIdx), [2001, 135]);
ps_I = reshape(averageMUAMatrix_BRFSps(:,11:15,useIdx), [2001, 135]);
ns_S = reshape(averageMUAMatrix_BRFSns(:,1:5,useIdx), [2001, 135]);
ns_G = reshape(averageMUAMatrix_BRFSns(:,6:10,useIdx), [2001, 135]);
ns_I = reshape(averageMUAMatrix_BRFSns(:,11:15,useIdx), [2001, 135]);

% Smoothed averages for plotting
ps_S_avg = smoothdata(mean(ps_S, 2, "omitmissing"), "gaussian", 50);
ps_G_avg = smoothdata(mean(ps_G, 2, "omitmissing"), "gaussian", 50);
ps_I_avg = smoothdata(mean(ps_I, 2, "omitmissing"), "gaussian", 50);
ns_S_avg = smoothdata(mean(ns_S, 2, "omitmissing"), "gaussian", 50);
ns_G_avg = smoothdata(mean(ns_G, 2, "omitmissing"), "gaussian", 50);
ns_I_avg = smoothdata(mean(ns_I, 2, "omitmissing"), "gaussian", 50);

ps_data_all = {ps_S, ps_G, ps_I}; % This takes the raw (unsmoothed) data)
ns_data_all = {ns_S, ns_G, ns_I};
ps_avg_all = {ps_S_avg, ps_G_avg, ps_I_avg};
ns_avg_all = {ns_S_avg, ns_G_avg, ns_I_avg};
%% Statistical Analysis for Significance Testing
tm_full = -200:1800;
bin_width = 50; % ms
time_bins = 0:bin_width:max(tm_full);
num_bins = length(time_bins) - 1;
original_threshold = 0.05;
bonferroni_threshold = (original_threshold / (3*num_bins));

compartments = {'Supragranular', 'Granular', 'Infragranular'};

lamCom = figure;
set(gcf, "Position", [1000 123.6667 757.6667 1.1140e+03]);
t = tiledlayout(3, 1);

%% Plot smoothed data with SEM
for idx = 1:3
    nexttile
    ps_data = ps_data_all{idx};
    ns_data = ns_data_all{idx};
    ps_avg = ps_avg_all{idx};
    ns_avg = ns_avg_all{idx};
    
    % Calculate SEM for each time point
    contactNum = size(ps_data, 2);
    ps_sem = std(ps_data, 0, 2, 'omitnan') ./ sqrt(contactNum);
    ns_sem = std(ns_data, 0, 2, 'omitnan') ./ sqrt(contactNum);
    
    % Plot shaded SEM regions
    fill([tm_full, fliplr(tm_full)], ...
         [ps_avg' + ps_sem', fliplr(ps_avg' - ps_sem')], ...
         [230/255 97/255 1/255], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on;
    fill([tm_full, fliplr(tm_full)], ...
         [ns_avg' + ns_sem', fliplr(ns_avg' - ns_sem')], ...
         [94/255 60/255 153/255], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    % Plot smoothed data
    plot(tm_full, ps_avg, 'color', [230/255 97/255 1/255], 'LineWidth', 1.5);
    plot(tm_full, ns_avg, 'color', [94/255 60/255 153/255], 'LineWidth', 1.5);


    % Plot asterisks for significant bins
    y_pos = max(max(ps_avg), max(ns_avg)) - 0.1 * range([ps_avg(:); ns_avg(:)]); 
    for i = 1:num_bins
        bin_indices = find(tm_full >= time_bins(i) & tm_full < time_bins(i+1));
        
        % Extract raw data for the current bin
        ps_bin_data = mean(ps_data(bin_indices, :), 1, 'omitnan'); % Mean across time bin
        ns_bin_data = mean(ns_data(bin_indices, :), 1, 'omitnan');
        
        % Perform t-test across electrodes
        [~, p] = ttest2(ps_bin_data, ns_bin_data);
        
        % Annotate significance if p < Bonferroni threshold
        if p < bonferroni_threshold
            x_pos = mean(time_bins(i:i+1));
            text(x_pos, y_pos, '*', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'center');
        end
    end
    
    xlabel('Time (ms)');
    ylabel('Percent Change');
    ylim([0 35]);
    title([compartments{idx}, ' Compartment']);
    % Add black lines for stimulus onset and offset times
    xline(0, 'k', 'LineWidth', 2);
    xline(800, 'k', 'LineWidth', 2);
    xline(1600, 'k', 'LineWidth', 2);
end

title(t, 'Laminar Compartmental MUA Responses');


%% Save figure
answer = questdlg('Would you like to save this figure?', 'Y', 'N');
switch answer
    case 'Yes'
        cd(plotDir)
        saveas(lamCom, 'fig4b_BRFSlaminar.png');
        saveas(lamCom, 'fig4b_BRFSlaminar.svg');
    case 'No'
        disp('Figure not saved.')
end


