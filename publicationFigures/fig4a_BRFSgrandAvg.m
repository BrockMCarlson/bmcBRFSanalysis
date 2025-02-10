%% fig3a
% MUA grand average

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
   

   

    %% MUA in percent change
    % BRFS pref vs null
    % Get the number of trials for the chosen condition
    for i = 1:length(v1Ch)
        if ~tuned(i,1)
            continue
        end
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


%% Organize into compartments, median and std across contacts+penetration
% reshape
useIdx = squeeze(~isnan(averageMUAMatrix_BRFSps(1,1,:))); 
ps_reshaped = reshape(averageMUAMatrix_BRFSps(:,:,useIdx),[2001,405]);
ns_reshaped = reshape(averageMUAMatrix_BRFSns(:,:,useIdx),[2001,405]);

% Average across penetrations
psAvg = smoothdata(mean(ps_reshaped,2,"omitmissing"),"gaussian",20);
nsAvg = smoothdata(mean(ns_reshaped,2,"omitmissing"),"gaussian",20);



% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
contactNum = size(ps_reshaped,2);
psSEM = std(ps_reshaped,0,2,"omitmissing")./sqrt(contactNum); 
nsSEM = std(ns_reshaped,0,2,"omitmissing")./sqrt(contactNum); 

% Mean +/- SEM
ps_S_avgPlusSEM = psAvg + psSEM; %pref stim -- median plus sem 1 
ps_S_avgMinusSEM = psAvg - psSEM; %pref stim -- median minus sem 1 


ns_S_avgPlusSEM = nsAvg + nsSEM; 
ns_S_avgMinusSEM = nsAvg - nsSEM;  



%% Figure generation! 

tm_full = (-200:1800)'; % 1801 total timepoints
bin_width = 50; % ms
time_bins = 0:bin_width:max(tm_full);
num_bins = length(time_bins) - 1;
original_threshold = 0.05;
bonferroni_threshold = original_threshold / num_bins;

f = figure;
hold on;

ps_avg_data = psAvg;
ns_avg_data = nsAvg;
ps_plus = ps_S_avgPlusSEM;
ps_minus = ps_S_avgMinusSEM;
ns_plus = ns_S_avgPlusSEM;
ns_minus = ns_S_avgMinusSEM;



% Plot the shaded region for dioptic (psAvg) in light blue
fill([tm_full; flipud(tm_full)], [ps_plus; flipud(ps_minus)], [230/255 97/255 1/255], 'FaceAlpha', 0.4, 'EdgeColor', [0 0 0.5]); 
hold on

% Plot the shaded region for dichoptic (nsAvg) in light red
fill([tm_full; flipud(tm_full)], [ns_plus; flipud(ns_minus)], [94/255 60/255 153/255], 'FaceAlpha', 0.4, 'EdgeColor', [0.5 0 0]); 

% Plot main lines for dioptic and dichoptic
plot(tm_full, ps_avg_data, 'Color', [230/255 97/255 1/255], 'LineWidth', 1.5); 
plot(tm_full, ns_avg_data, 'Color', [94/255 60/255 153/255], 'LineWidth', 1.5); 

% Add black lines for stimulus onset and offset times
xline(0, 'k', 'LineWidth', 2);
xline(800, 'k', 'LineWidth', 2);
xline(1600, 'k', 'LineWidth', 2);

% Bonferroni-adjusted significance testing
y_pos = max(max(ps_avg_data), max(ns_avg_data)) - 0.1 * range([ps_avg_data(:); ns_avg_data(:)]); 
for i = 1:num_bins
    bin_indices = find(tm_full >= time_bins(i) & tm_full < time_bins(i+1));
    ps_data_bin = ps_avg_data(bin_indices);
    ns_data_bin = ns_avg_data(bin_indices);
    [~, p] = ttest2(ps_data_bin, ns_data_bin);
    if p < bonferroni_threshold
        x_pos = mean(time_bins(i:i+1)); 
        text(x_pos, y_pos, '*', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'center');
    end
end

xlabel('Time (ms)');
ylabel('Percent Change');


% Add global title
title('Laminar Compartmental MUA Responses');

%% save fig
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
       disp('alright, saving figure to plotdir')
        cd(plotDir)
        saveas(f, 'fig3a_BRFSgrandAvg.png');
        saveas(f, 'fig3a_BRFSgrandAvg.svg');
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end


