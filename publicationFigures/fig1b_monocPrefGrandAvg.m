%% fig1b
% MonocPrefGrandAverage

%% Setup
disp('start time')
datetime
close
clearvars -except MUA_trials
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir\fig3_MUA';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\fig3_MUA';
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
averageMUAMatrix_POPE = nan(2001,15,31);
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
        % Adjust the trial-type assignments
        if prefMonoc(i, 1) == 1
            prefOri_prefEye{i, 1} = [5, 11, 13, 19]; % PO RightEye
            prefOri_nullEye{i, 1} = [8, 10, 16, 18]; % PO LeftEye
            nullOri_prefEye{i, 1} = [7, 9, 15, 17];  % NPO RightEye
            nullOri_nullEye{i, 1} = [6, 12, 14, 20]; % NPO LeftEye
        elseif prefMonoc(i, 1) == 2
            prefOri_prefEye{i, 1} = [8, 10, 16, 18]; % PO LeftEye
            prefOri_nullEye{i, 1} = [5, 11, 13, 19]; % PO RightEye
            nullOri_prefEye{i, 1} = [6, 12, 14, 20]; % NPO LeftEye
            nullOri_nullEye{i, 1} = [7, 9, 15, 17];  % NPO RightEye
        elseif prefMonoc(i, 1) == 3
            prefOri_prefEye{i, 1} = [7, 9, 15, 17];  % NPO RightEye
            prefOri_nullEye{i, 1} = [6, 12, 14, 20]; % NPO LeftEye
            nullOri_prefEye{i, 1} = [5, 11, 13, 19]; % PO RightEye
            nullOri_nullEye{i, 1} = [8, 10, 16, 18]; % PO LeftEye
        elseif prefMonoc(i, 1) == 4
            prefOri_prefEye{i, 1} = [6, 12, 14, 20]; % NPO LeftEye
            prefOri_nullEye{i, 1} = [7, 9, 15, 17];  % NPO RightEye
            nullOri_prefEye{i, 1} = [8, 10, 16, 18]; % PO LeftEye
            nullOri_nullEye{i, 1} = [5, 11, 13, 19]; % PO RightEye
        end
    end
    %% Average  
    % Preallocate matrices for the four conditions
    averageMUAMatrix_POPE = nan(2001, length(v1Ch), size(MUA_trials, 1));
    averageMUAMatrix_PONE = nan(2001, length(v1Ch), size(MUA_trials, 1));
    averageMUAMatrix_NOPE = nan(2001, length(v1Ch), size(MUA_trials, 1));
    averageMUAMatrix_NONE = nan(2001, length(v1Ch), size(MUA_trials, 1));
    
    for i = 1:length(v1Ch)
        % Define the conditions and their corresponding trial groups
        conditionIndices = {prefOri_prefEye{i, 1}, prefOri_nullEye{i, 1}, ...
                            nullOri_prefEye{i, 1}, nullOri_nullEye{i, 1}};
        
        % Loop through each condition
        for cIdx = 1:numel(conditionIndices)
            % Get the trial types for the condition
            trialTypes = conditionIndices{cIdx};
            MUAflashOut = [];
            
            % Loop through each trial type and accumulate trials
            for trialType = trialTypes
                if isempty(MUA_trials{penetration, 1}{trialType, 1})
                    continue;
                end
                
                numTrials = size(MUA_trials{penetration, 1}{trialType, 1}, 1);
                for trl = 1:numTrials
                    flashData = MUA_trials{penetration, 1}{trialType, 1}{trl, 1}(:, v1Ch(i));
                    MUAflashOut = cat(2, MUAflashOut, flashData); % Concatenate trials
                end
            end
            
            % Calculate the median response
            avg = median(MUAflashOut, 2, "omitmissing");
            bl = median(avg(1:200, :)); % Baseline average
            if bl == 0
                bl = 0.1; % Avoid division by zero
            end
            
            % Calculate percent change
            percentChange = 100 * ((avg - bl) ./ bl);
            
            % Assign the results to the corresponding matrix
            if cIdx == 1
                averageMUAMatrix_POPE(:, i, penetration) = percentChange;
            elseif cIdx == 2
                averageMUAMatrix_PONE(:, i, penetration) = percentChange;
            elseif cIdx == 3
                averageMUAMatrix_NOPE(:, i, penetration) = percentChange;
            elseif cIdx == 4
                averageMUAMatrix_NONE(:, i, penetration) = percentChange;
            end
        end
    end



%% end of the for loop here
    disp(strcat('Done with file number: ',string(penetration)))
end

%% Median and standard error across contacts and penetrations
% Reshape matrices for all four conditions
useIdx = squeeze(~isnan(averageMUAMatrix_POPE(1, 1, :)));

% Reshape dynamically based on non-NaN data
ps_reshaped = reshape(averageMUAMatrix_POPE(:, :, useIdx), 2001, []);
pone_reshaped = reshape(averageMUAMatrix_PONE(:, :, useIdx), 2001, []);
nope_reshaped = reshape(averageMUAMatrix_NOPE(:, :, useIdx), 2001, []);
none_reshaped = reshape(averageMUAMatrix_NONE(:, :, useIdx), 2001, []);

% Average across penetrations for all four conditions
psAvg = smoothdata(mean(ps_reshaped, 2, "omitmissing"), "gaussian", 20);
poneAvg = smoothdata(mean(pone_reshaped, 2, "omitmissing"), "gaussian", 20);
nopeAvg = smoothdata(mean(nope_reshaped, 2, "omitmissing"), "gaussian", 20);
noneAvg = smoothdata(mean(none_reshaped, 2, "omitmissing"), "gaussian", 20);

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
contactNum = size(ps_reshaped, 2);
psSEM = std(ps_reshaped, 0, 2, "omitmissing") ./ sqrt(contactNum);
poneSEM = std(pone_reshaped, 0, 2, "omitmissing") ./ sqrt(contactNum);
nopeSEM = std(nope_reshaped, 0, 2, "omitmissing") ./ sqrt(contactNum);
noneSEM = std(none_reshaped, 0, 2, "omitmissing") ./ sqrt(contactNum);

% Mean +/- SEM for each condition
ps_S_avgPlusSEM = psAvg + psSEM; 
ps_S_avgMinusSEM = psAvg - psSEM;
pone_S_avgPlusSEM = poneAvg + poneSEM; 
pone_S_avgMinusSEM = poneAvg - poneSEM;
nope_S_avgPlusSEM = nopeAvg + nopeSEM; 
nope_S_avgMinusSEM = nopeAvg - nopeSEM;
none_S_avgPlusSEM = noneAvg + noneSEM; 
none_S_avgMinusSEM = noneAvg - noneSEM;

%% Figure generation for all conditions
close all;
tm_full = -200:1800; % 1801 total timepoints
% Trim time and data to the first 800ms
tm_trim = tm_full(1:1001); % Timepoints corresponding to the first 800ms
psAvg_trim = psAvg(1:1001);
poneAvg_trim = poneAvg(1:1001);
nopeAvg_trim = nopeAvg(1:1001);
noneAvg_trim = noneAvg(1:1001);

ps_S_avgPlusSEM_trim = ps_S_avgPlusSEM(1:1001);
ps_S_avgMinusSEM_trim = ps_S_avgMinusSEM(1:1001);

pone_S_avgPlusSEM_trim = pone_S_avgPlusSEM(1:1001);
pone_S_avgMinusSEM_trim = pone_S_avgMinusSEM(1:1001);

nope_S_avgPlusSEM_trim = nope_S_avgPlusSEM(1:1001);
nope_S_avgMinusSEM_trim = nope_S_avgMinusSEM(1:1001);

none_S_avgPlusSEM_trim = none_S_avgPlusSEM(1:1001);
none_S_avgMinusSEM_trim = none_S_avgMinusSEM(1:1001);

% Generate the plot
lamCom = figure;
hold on;

% Plot each condition with different colors
plot(tm_trim, psAvg_trim, 'color', [230/255 97/255 1/255], 'LineWidth', 1.5); % prefOri_prefEye
plot(tm_trim, poneAvg_trim, 'color', [94/255 153/255 60/255], 'LineWidth', 1.5); % prefOri_nullEye
plot(tm_trim, nopeAvg_trim, 'color', [60/255 153/255 230/255], 'LineWidth', 1.5); % nullOri_prefEye
plot(tm_trim, noneAvg_trim, 'color', [153/255 60/255 94/255], 'LineWidth', 1.5); % nullOri_nullEye

% Plot SEM as shaded areas or dotted lines
plot(tm_trim, ps_S_avgPlusSEM_trim, 'color', [230/255 97/255 1/255], 'LineWidth', 1, 'LineStyle', ':');
plot(tm_trim, ps_S_avgMinusSEM_trim, 'color', [230/255 97/255 1/255], 'LineWidth', 1, 'LineStyle', ':');

plot(tm_trim, pone_S_avgPlusSEM_trim, 'color', [94/255 153/255 60/255], 'LineWidth', 1, 'LineStyle', ':');
plot(tm_trim, pone_S_avgMinusSEM_trim, 'color', [94/255 153/255 60/255], 'LineWidth', 1, 'LineStyle', ':');

plot(tm_trim, nope_S_avgPlusSEM_trim, 'color', [60/255 153/255 230/255], 'LineWidth', 1, 'LineStyle', ':');
plot(tm_trim, nope_S_avgMinusSEM_trim, 'color', [60/255 153/255 230/255], 'LineWidth', 1, 'LineStyle', ':');

plot(tm_trim, none_S_avgPlusSEM_trim, 'color', [153/255 60/255 94/255], 'LineWidth', 1, 'LineStyle', ':');
plot(tm_trim, none_S_avgMinusSEM_trim, 'color', [153/255 60/255 94/255], 'LineWidth', 1, 'LineStyle', ':');

vline(0); % Add vertical line at 0 ms
title('Grand Averages for Four Conditions (First 800ms)');
xlabel('Time (ms)');
ylabel('% Change from Baseline');
box("off");
legend('prefOri_prefEye', 'prefOri_nullEye', 'nullOri_prefEye', 'nullOri_nullEye','Interpreter', 'none');


%% ave fig
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
       disp('alright, saving figure to plotdir')
        cd(plotDir)
        figName_lamCom = strcat('MUA_fig1b','_MonocAvg_','.png');
        saveas(lamCom,figName_lamCom)
        figName_lamCom = strcat('MUA_fig1b','_MonocAvg_','.svg');
        saveas(lamCom,figName_lamCom)
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end


