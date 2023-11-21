% Plot average of all conditions


%% Setup
clear
% sessionLabel = '211008_B_bmcBRFS001';
sessionLabel = '221206_J_bmcBRFS003';
% sessionLabel = '221102_J_bmcBRFS001'; % the problem session

% Directories
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(outputDir)

% load data
load(strcat('sortedData_',sessionLabel,'.mat'))


%% bl Sub at average level (better for plotting)
ch = 1:32;

% Monocular
monoc_1 = [5, 11, 13, 19]; % PO RightEye
monoc_2 = [8, 10, 16, 18]; % PO LeftEye
monoc_3 = [7, 9, 15, 17];  % NPO RightEye
monoc_4 = [6, 12, 14, 20]; % NPO LeftEye

% convert from cell to double and combine monocular conditions
count = 0;
for cond = monoc_1
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc1(:,:,count) = IDX(cond).LFP_gamma{trl,1}(:,ch); % 1000 x 32
    end
end
clear cond count trl 

count = 0;
for cond = monoc_2
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc2(:,:,trl) = IDX(cond).LFP_gamma{trl,1}(:,ch); 
    end
end
clear cond count trl 


count = 0;
for cond = monoc_3
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc3(:,:,trl) = IDX(cond).LFP_gamma{trl,1}(:,ch); 
    end
end
clear cond count trl 


count = 0;
for cond = monoc_4
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc4(:,:,trl) = IDX(cond).LFP_gamma{trl,1}(:,ch); 
    end
end
clear cond count trl 





%% Test for monocular preference with anova
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
tuned = nan(length(ch),1);
prefMonoc = nan(length(ch),1);

% FOR each individual unit 
for i = ch
    % preallocate y based on trial count
    y = nan(maxTrls,4);
    % create an array of each trial's median response for each monoc 
    % condition (trl x 4 monoc)
    y(1:size(array_ofMonoc1,3),1) = median(squeeze(array_ofMonoc1(200:300,i,:)),1)'; % median of each trial after stim onset
    y(1:size(array_ofMonoc2,3),2) = median(squeeze(array_ofMonoc2(200:300,i,:)),1)'; % median of each trial after stim onset
    y(1:size(array_ofMonoc3,3),3) = median(squeeze(array_ofMonoc3(200:300,i,:)),1)'; % median of each trial after stim onset
    y(1:size(array_ofMonoc4,3),4) = median(squeeze(array_ofMonoc4(200:300,i,:)),1)'; % median of each trial after stim onset
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

%% Plot monocular based on monocular preference
avg_ofMonoc(:,:,1) = median(array_ofMonoc1,3); 
avg_ofMonoc(:,:,2) = median(array_ofMonoc2,3); 
avg_ofMonoc(:,:,3) = median(array_ofMonoc3,3); 
avg_ofMonoc(:,:,4) = median(array_ofMonoc4,3); 

for i = ch
    sorted_ofMonoc_pref(:,i) = avg_ofMonoc(:,i,prefMonoc(i,1));
    sorted_ofMonoc_null(:,i) = avg_ofMonoc(:,i,nullMonoc(i,1));
end

% Average over all contacts
GrandAvg_Monoc_pref = median(sorted_ofMonoc_pref,2);
GrandAvg_Monoc_null = median(sorted_ofMonoc_null,2);

% Baseline subtractions
bl_pref = median(GrandAvg_Monoc_pref(100:200,1));
bl_null = median(GrandAvg_Monoc_null(100:200,1));
blSub_pref = GrandAvg_Monoc_pref - bl_pref;
blSub_null = GrandAvg_Monoc_null - bl_null;

% Smooth data
smooth_ofMonoc_pref = smoothdata(blSub_pref,"gaussian",20);
smooth_ofMonoc_null = smoothdata(blSub_null,"gaussian",20);


% Plot
figure
set(gcf,"Position",[17.6667 59 2.5313e+03 1.2813e+03])
t = tiledlayout(3,3,'TileSpacing','compact');
set(0, 'DefaultLineLineWidth', 1.5);
nexttile
tm = -200:800;
plot(tm,smooth_ofMonoc_pref); hold on
plot(tm,smooth_ofMonoc_null); hold on
vline(0)
legend('Preferred stimulus','Null stimulus')
title('Average of preferred and null monocular stimulus across units')
box off

%% binocular dioptic
% concatenate two timecourses into array
tm_full = -200:1600; % 1801 total timepoints
tm1 = 1:801;
tm2 = 1:1001;
tm2_concat = 801:1801;

for trl = 1:length(IDX(1).correctTrialIndex)
    array_bi_dioptic_1(tm1,:,trl) = IDX(1).LFP_gamma{trl,1}(tm1,ch); 
    array_bi_dioptic_1(tm2_concat,:,trl) = IDX(1).LFP_gamma{trl,2}(tm2,ch); 
end

for trl = 1:length(IDX(2).correctTrialIndex)
    array_bi_dioptic_2(tm1,:,trl) = IDX(2).LFP_gamma{trl,1}(tm1,ch); 
    array_bi_dioptic_2(tm2_concat,:,trl) = IDX(2).LFP_gamma{trl,2}(tm2,ch); 
end

% Align arrays based on preference
for j = ch
    if prefMonoc(j,1) <= 2
        prefBinoc(j,1) = 1;
        nullBinoc(j,1) = 2;
    elseif prefMonoc(j,1) >= 3
        prefBinoc(j,1) = 2;
        nullBinoc(j,1) = 1;
    end
end

avg_bi_dioptic(:,:,1) = median(array_bi_dioptic_1,3); % average across trials
avg_bi_dioptic(:,:,2) = median(array_bi_dioptic_2,3); % average across trials


for i = ch
    sorted_bi_dioptic_pref(:,i) = avg_bi_dioptic(:,i,prefBinoc(i,1));
    sorted_bi_dioptic_null(:,i) = avg_bi_dioptic(:,i,nullBinoc(i,1));
end

% Average over all contacts
GrandAvg_Binoc_pref = median(sorted_bi_dioptic_pref,2);
GrandAvg_Binoc_null = median(sorted_bi_dioptic_null,2);

% Baseline subtractions
bl_Binoc_pref = median(GrandAvg_Binoc_pref(100:200,1));
bl_Binoc_null = median(GrandAvg_Binoc_null(100:200,1));
blSub_Binoc_pref = GrandAvg_Binoc_pref - bl_Binoc_pref;
blSub_Binoc_null = GrandAvg_Binoc_null - bl_Binoc_null;

% Smooth data
smooth_Binoc_pref = smoothdata(blSub_Binoc_pref,"gaussian",20);
smooth_Binoc_null = smoothdata(blSub_Binoc_null,"gaussian",20);


% plot
nexttile
plot(tm_full,smooth_Binoc_pref); hold on
plot(tm_full,smooth_Binoc_null); hold on
vline(0); vline(800)
legend('binoc dioptic PO','binoc dioptic NPO')
title('Binocular onset - dioptic')
xlim([-200 1600])
box off

%% binocular dichoptic
for trl = 1:length(IDX(3).correctTrialIndex)
    array_bi_dichoptic_1(tm1,:,trl) = IDX(3).LFP_gamma{trl,1}(tm1,ch); 
    array_bi_dichoptic_1(tm2_concat,:,trl) = IDX(3).LFP_gamma{trl,2}(tm2,ch); 
end

for trl = 1:length(IDX(4).correctTrialIndex)
    array_bi_dichoptic_2(tm1,:,trl) = IDX(4).LFP_gamma{trl,1}(tm1,ch); 
    array_bi_dichoptic_2(tm2_concat,:,trl) = IDX(4).LFP_gamma{trl,2}(tm2,ch); 
end


% Align arrays based on preference
for j = ch
    if prefMonoc(j,1) == 1 || prefMonoc(j,1) == 4 
        prefBinoc_dichoptic(j,1) = 2;
        nullBinoc_dichoptic(j,1) = 1;
    elseif prefMonoc(j,1) == 2 || prefMonoc(j,1) == 3
        prefBinoc_dichoptic(j,1) = 1;
        nullBinoc_dichoptic(j,1) = 2;
    end
end

avg_bi_dichoptic(:,:,1) = median(array_bi_dichoptic_1,3); % average across trials
avg_bi_dichoptic(:,:,2) = median(array_bi_dichoptic_2,3); % average across trials


for i = ch
    sorted_bi_dichoptic_pref(:,i) = avg_bi_dichoptic(:,i,prefBinoc_dichoptic(i,1));
    sorted_bi_dichoptic_null(:,i) = avg_bi_dichoptic(:,i,nullBinoc_dichoptic(i,1));
end

% Average over all contacts
GrandAvg_bi_dichoptic_pref = median(sorted_bi_dichoptic_pref,2);
GrandAvg_bi_dichoptic_null = median(sorted_bi_dichoptic_null,2);

% Baseline subtractions
bl_bi_dichoptic_pref = median(GrandAvg_bi_dichoptic_pref(100:200,1));
bl_bi_dichoptic_null = median(GrandAvg_bi_dichoptic_null(100:200,1));
blSub_bi_dichoptic_pref = GrandAvg_bi_dichoptic_pref - bl_bi_dichoptic_pref;
blSub_bi_dichoptic_null = GrandAvg_bi_dichoptic_null - bl_bi_dichoptic_null;

% Smooth data
smooth_bi_dichoptic_pref = smoothdata(blSub_bi_dichoptic_pref,"gaussian",20);
smooth_bi_dichoptic_null = smoothdata(blSub_bi_dichoptic_null,"gaussian",20);


% plot
nexttile
plot(tm_full,smooth_bi_dichoptic_pref); hold on
plot(tm_full,smooth_bi_dichoptic_null); hold on
vline(0); vline(800)
legend('binoc dichoptic 1','binoc dichoptic 2')
title('Binocular onset - dichoptic')
xlim([-200 1600])
box off

%% BRFS-like dioptic
% Ok, you have a good formula. Now lets try to do this programmatically
nexttile(5)
for i = ch
    % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
    % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    if prefMonoc(i,1) == 1
        prefCondOnFlash = 8;
        nullCondOnFlash = 5;
    elseif prefMonoc(i,1) == 2
        prefCondOnFlash = 5;
        nullCondOnFlash = 8;
    elseif prefMonoc(i,1) == 3
        prefCondOnFlash = 6;
        nullCondOnFlash = 7;
    elseif prefMonoc(i,1) == 4
        prefCondOnFlash = 7;
        nullCondOnFlash = 6;
    end

    % convert cell to array
    for trl = 1:length(IDX(prefCondOnFlash).correctTrialIndex)
        array_diopticAdapted_pref(tm1,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_diopticAdapted_pref(tm2_concat,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
    for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
        array_diopticAdapted_null(tm1,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_diopticAdapted_null(tm2_concat,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
end

% and average across whole electrode
avg_diopticAdapted_pref = median(median(array_diopticAdapted_pref,3),2); % average across trials and average across electrode
bl_diopticAdapted_pref(:,1) = avg_diopticAdapted_pref(:,1) - median(avg_diopticAdapted_pref(100:200,1));
smooth_diopticAdapted_pref(:,1) = smoothdata(bl_diopticAdapted_pref(:,1),"gaussian",20);
plot(tm_full,smooth_diopticAdapted_pref); hold on
avg_diopticAdapted_null= median(median(array_diopticAdapted_null,3),2); % average across trials and average across electrode
bl_diopticAdapted_null(:,1) = avg_diopticAdapted_null(:,1) - median(avg_diopticAdapted_null(100:200,1));
smooth_diopticAdapted_null(:,1) = smoothdata(bl_diopticAdapted_null(:,1),"gaussian",20);
plot(tm_full,smooth_diopticAdapted_null); hold on

vline(0); vline(800)
legend('Preferred condition flashed','Same ori other eye flashed')
title('BRFS-like adaptation - dioptic')
xlim([-200 1600])
box off

%% BRFS (dichoptic)
nexttile(6)
% preallocate
maxTrlLength = max([length(IDX(9).correctTrialIndex),...
    length(IDX(10).correctTrialIndex),...
    length(IDX(11).correctTrialIndex),...
    length(IDX(12).correctTrialIndex)]);
array_dichopticAdapted_pref = nan(1801,32,maxTrlLength);
array_dichopticAdapted_null = nan(1801,32,maxTrlLength);
for i = ch
    % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
    % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    if prefMonoc(i,1) == 1
        prefCondOnFlash = 12;
        nullCondOnFlash = 11;
    elseif prefMonoc(i,1) == 2
        prefCondOnFlash = 9;
        nullCondOnFlash = 10;
    elseif prefMonoc(i,1) == 3
        prefCondOnFlash = 10;
        nullCondOnFlash = 9;
    elseif prefMonoc(i,1) == 4
        prefCondOnFlash = 11;
        nullCondOnFlash = 12;
    end

    % convert cell to array
    for trl = 1:length(IDX(prefCondOnFlash).correctTrialIndex)
        array_dichopticAdapted_pref(tm1,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_dichopticAdapted_pref(tm2_concat,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
    for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
        array_dichopticAdapted_null(tm1,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_dichopticAdapted_null(tm2_concat,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
end

% and average across whole electrode
avg_dichopticAdapted_pref = median(median(array_dichopticAdapted_pref,3,"omitmissing"),2,"omitmissing"); % average across trials and average across electrode
bl_dichopticAdapted_pref(:,1) = avg_dichopticAdapted_pref(:,1) - median(avg_dichopticAdapted_pref(100:200,1));
smooth_dichopticAdapted_pref(:,1) = smoothdata(bl_dichopticAdapted_pref(:,1),"gaussian",20);
plot(tm_full,smooth_dichopticAdapted_pref); hold on
avg_dichopticAdapted_null= median(median(array_dichopticAdapted_null,3,"omitmissing"),2,"omitmissing");% average across trials and average across electrode
bl_dichopticAdapted_null(:,1) = avg_dichopticAdapted_null(:,1) - median(avg_dichopticAdapted_null(100:200,1));
smooth_dichopticAdapted_null(:,1) = smoothdata(bl_dichopticAdapted_null(:,1),"gaussian",20);
plot(tm_full,smooth_dichopticAdapted_null); hold on

vline(0); vline(800)
legend('Preferred condition flashed','Same simulus, but null flashed')
title('BRFS')
xlim([-200 1600])
box off


%% Monocular alternation dioptic
nexttile(8)
% preallocate
maxTrlLength = max([length(IDX(13).correctTrialIndex),...
    length(IDX(14).correctTrialIndex),...
    length(IDX(15).correctTrialIndex),...
    length(IDX(16).correctTrialIndex)]);
array_diopticMonocAlt_pref = nan(1801,32,maxTrlLength);
array_diopticMonocAlt_null = nan(1801,32,maxTrlLength);
for i = ch
    % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
    % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    if prefMonoc(i,1) == 1
        prefCondOnFlash = 16;
        nullCondOnFlash = 13;
    elseif prefMonoc(i,1) == 2
        prefCondOnFlash = 13;
        nullCondOnFlash = 16;
    elseif prefMonoc(i,1) == 3
        prefCondOnFlash = 14;
        nullCondOnFlash = 15;
    elseif prefMonoc(i,1) == 4
        prefCondOnFlash = 15;
        nullCondOnFlash = 14;
    end

    % convert cell to array
    for trl = 1:length(IDX(prefCondOnFlash).correctTrialIndex)
        array_diopticMonocAlt_pref(tm1,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_diopticMonocAlt_pref(tm2_concat,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
    for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
        array_diopticMonocAlt_null(tm1,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_diopticMonocAlt_null(tm2_concat,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
end

% and average across whole electrode
avg_diopticMonocAlt_pref = median(median(array_diopticMonocAlt_pref,3,"omitmissing"),2,"omitmissing"); % average across trials and average across electrode
bl_diopticMonocAlt_pref(:,1) = avg_diopticMonocAlt_pref(:,1) - median(avg_diopticMonocAlt_pref(100:200,1));
smooth_diopticMonocAlt_pref(:,1) = smoothdata(bl_diopticMonocAlt_pref(:,1),"gaussian",20);
plot(tm_full,smooth_diopticMonocAlt_pref); hold on
avg_diopticMonocAlt_null= median(median(array_diopticMonocAlt_null,3,"omitmissing"),2,"omitmissing");% average across trials and average across electrode
bl_diopticMonocAlt_null(:,1) = avg_diopticMonocAlt_null(:,1) - median(avg_diopticMonocAlt_null(100:200,1));
smooth_diopticMonocAlt_null(:,1) = smoothdata(bl_diopticMonocAlt_null(:,1),"gaussian",20);
plot(tm_full,smooth_diopticMonocAlt_null); hold on


vline(0); vline(800);
legend('Preferred condition flashed','Same ori, null eye flashed')
title('Monocular alternation - dioptic')
xlim([-200 1600])
box off

%% Monocular alternation dichoptic
nexttile(9)

% preallocate
maxTrlLength = max([length(IDX(17).correctTrialIndex),...
    length(IDX(18).correctTrialIndex),...
    length(IDX(19).correctTrialIndex),...
    length(IDX(20).correctTrialIndex)]);
array_dichopticMonocAlt_pref = nan(1801,32,maxTrlLength);
array_dichopticMonocAlt_null = nan(1801,32,maxTrlLength);
for i = ch
    % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
    % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    if prefMonoc(i,1) == 1
        prefCondOnFlash = 20;
        nullCondOnFlash = 19;
    elseif prefMonoc(i,1) == 2
        prefCondOnFlash = 17;
        nullCondOnFlash = 18;
    elseif prefMonoc(i,1) == 3
        prefCondOnFlash = 18;
        nullCondOnFlash = 17;
    elseif prefMonoc(i,1) == 4
        prefCondOnFlash = 19;
        nullCondOnFlash = 20;
    end

    % convert cell to array
    for trl = 1:length(IDX(prefCondOnFlash).correctTrialIndex)
        array_dichopticMonocAlt_pref(tm1,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_dichopticMonocAlt_pref(tm2_concat,i,trl) = IDX(prefCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
    for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
        array_dichopticMonocAlt_null(tm1,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,1}(tm1,i); % now we index the cell for the second 800ms
        array_dichopticMonocAlt_null(tm2_concat,i,trl) = IDX(nullCondOnFlash).LFP_gamma{trl,2}(tm2,i); % now we index the cell for the second 800ms
    end
end

% and average across whole electrode
avg_dichopticMonocAlt_pref = median(median(array_dichopticMonocAlt_pref,3,"omitmissing"),2,"omitmissing"); % average across trials and average across electrode
bl_dichopticMonocAlt_pref(:,1) = avg_dichopticMonocAlt_pref(:,1) - median(avg_dichopticMonocAlt_pref(100:200,1));
smooth_dichopticMonocAlt_pref(:,1) = smoothdata(bl_dichopticMonocAlt_pref(:,1),"gaussian",20);
plot(tm_full,smooth_dichopticMonocAlt_pref); hold on
avg_dichopticMonocAlt_null= median(median(array_dichopticMonocAlt_null,3,"omitmissing"),2,"omitmissing");% average across trials and average across electrode
bl_dichopticMonocAlt_null(:,1) = avg_dichopticMonocAlt_null(:,1) - median(avg_dichopticMonocAlt_null(100:200,1));
smooth_dichopticMonocAlt_null(:,1) = smoothdata(bl_dichopticMonocAlt_null(:,1),"gaussian",20);
plot(tm_full,smooth_dichopticMonocAlt_null); hold on

vline(0); vline(800);
legend('Preferred condition flashed','Same simulus, but null flashed')
title('Monocular alternation - dichoptic')
xlim([-200 1600])
box off

%% Figure title and axis
title(t,sessionLabel,'Interpreter','none')
xlabel(t,'Time (ms)')
ylabel(t,'Neural reponse (uV)')

