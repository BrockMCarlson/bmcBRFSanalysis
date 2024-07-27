%% assignPrefFromMUA.m
% Are there observable differences between trial-types with LFP CSD?
% initialize variables
clearvars -except MUA_trials
close all

saveName = 'eventPref_211008_B_bmcBRFS001.mat';
outDir = 'D:\KiloSort-ed\211008_B_bmcBRFS001';

%% For loop
plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\CSDFigs';
dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\sortedData_240229';
% % dataDir = 'D:\sortedData_240229';
% % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
% % tic %takes just over 2 min
if ~exist('MUA_trials','var')
    load('MUA_trials.mat') % format is MUA_trials{1,penetration}{cond,1}{trial,flash}
end
% % toc
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);



%% Former for loop
penetration = 1;
    
probeName = char(officLamAssign.Session_probe_(penetration,1));
fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');

granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
if isnan(granBtm)
    warning(strcat('no sink found on _',fileToLoad))
    error()
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
        error()
    end
elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
    if any(v1Ch > 64) || any(v1Ch < 33)
        warning('skipping session without full column for now')
        disp(fileToLoad)
        error()
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
    for trl = 1:size(MUA_trials{1,penetration}{cond,1},1)
        count = count + 1;
        array_ofMonoc1(:,:,count) = abs(MUA_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
    end
end
clear cond count trl 

count = 0;
for cond = monoc_2
    for trl = 1:size(MUA_trials{1,penetration}{cond,1},1)
        count = count + 1;
        array_ofMonoc2(:,:,count) = abs(MUA_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); 
    end
end
clear cond count trl 


count = 0;
for cond = monoc_3
    for trl = 1:size(MUA_trials{1,penetration}{cond,1},1)
        count = count + 1;
        array_ofMonoc3(:,:,count) = abs(MUA_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); 
    end
end
clear cond count trl 


count = 0;
for cond = monoc_4
    for trl = 1:size(MUA_trials{1,penetration}{cond,1},1)
        count = count + 1;
        array_ofMonoc4(:,:,count) = abs(MUA_trials{1,penetration}{cond,1}{trl,1}(:,v1Ch)); 
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



%% save output
chAsStr = string(num2str(v1Ch'));
prefOutTable = table(prefCondOnFlash,nullCondOnFlash,'RowNames',chAsStr);

disp('saving session output')
cd(outDir)
save(saveName,'prefOutTable','overallPref','overallNull','-v7.3');

disp(strcat('Finished processing and saved file...',saveName))