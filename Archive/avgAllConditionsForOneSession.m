% Plot average of all conditions


%% Setup
clear
sessionLabel = '211008_B_bmcBRFS001';
% sessionLabel = '221206_J_bmcBRFS003';
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
        array_ofMonoc1(:,:,count) = IDX(cond).MUAe{trl,1}(:,ch); % 1000 x 32
    end
end
avg_ofMonoc1 = median(median(array_ofMonoc1,3),2); 
clear cond count trl 

count = 0;
for cond = monoc_2
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc2(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); 
    end
end
avg_ofMonoc2 = median(median(array_ofMonoc2,3),2); 
clear cond count trl 


count = 0;
for cond = monoc_3
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc3(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); 
    end
end
avg_ofMonoc3 = median(median(array_ofMonoc3,3),2);
clear cond count trl 


count = 0;
for cond = monoc_4
    for trl = 1:length(IDX(cond).correctTrialIndex)
        count = count + 1;
        array_ofMonoc4(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); 
    end
end
avg_ofMonoc4 = median(median(array_ofMonoc4,3),2);
clear cond count trl 

% Baseline subtractions
bl_1 = median(avg_ofMonoc1(100:200,1));
bl_2 = median(avg_ofMonoc2(100:200,1));
bl_3 = median(avg_ofMonoc3(100:200,1));
bl_4 = median(avg_ofMonoc4(100:200,1));
blSub1 = avg_ofMonoc1 - bl_1;
blSub2 = avg_ofMonoc2 - bl_2;
blSub3 = avg_ofMonoc3 - bl_3;
blSub4 = avg_ofMonoc4 - bl_4;

% Smooth data
smooth_ofMonoc1 = smoothdata(blSub1,"gaussian",20);
smooth_ofMonoc2 = smoothdata(blSub2,"gaussian",20);
smooth_ofMonoc3 = smoothdata(blSub3,"gaussian",20);
smooth_ofMonoc4 = smoothdata(blSub4,"gaussian",20);

% Plot
figure
set(gcf,"Position",[17.6667 59 2.5313e+03 1.2813e+03])
t = tiledlayout(3,3,'TileSpacing','compact');
set(0, 'DefaultLineLineWidth', 1.5);
nexttile
tm = -200:800;
plot(tm,smooth_ofMonoc1); hold on
plot(tm,smooth_ofMonoc2); hold on
plot(tm,smooth_ofMonoc3); hold on
plot(tm,smooth_ofMonoc4); hold on
vline(0)
legend('PO RightEye','PO LeftEye','NPO RightEye','NPO LeftEye')
title('Monocular')
box off

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
    [M,maxRespIdx] = max(median(y_bl,1));
    prefMonoc(i,1) = maxRespIdx;

end

%% binocular dioptic
% concatenate two timecourses into array
tm_full = -200:1600; % 1801 total timepoints
tm1 = 1:801;
tm2 = 1:1001;
tm2_concat = 801:1801;

for trl = 1:length(IDX(1).correctTrialIndex)
    array_bi_dioptic_PO(tm1,:,trl) = IDX(1).MUAe{trl,1}(tm1,ch); 
    array_bi_dioptic_PO(tm2_concat,:,trl) = IDX(1).MUAe{trl,2}(tm2,ch); 
    avg_bi_dioptic_PO = median(median(array_bi_dioptic_PO,3),2);
end
for trl = 1:length(IDX(2).correctTrialIndex)
    array_bi_dioptic_NPO(tm1,:,trl) = IDX(2).MUAe{trl,1}(tm1,ch); 
    array_bi_dioptic_NPO(tm2_concat,:,trl) = IDX(2).MUAe{trl,2}(tm2,ch); 
    avg_bi_dioptic_NPO = median(median(array_bi_dioptic_NPO,3),2);
end

% Baseline subtraction
bl_bi_dioptic_PO = avg_bi_dioptic_PO - median(avg_bi_dioptic_PO(100:200,1));
bl_bi_dioptic_NPO = avg_bi_dioptic_NPO - median(avg_bi_dioptic_NPO(100:200,1));

% Smooth data
smooth_bi_dioptic_PO = smoothdata(bl_bi_dioptic_PO,"gaussian",20);
smooth_bi_dioptic_NPO = smoothdata(bl_bi_dioptic_NPO,"gaussian",20);

% plot
nexttile
plot(tm_full,smooth_bi_dioptic_PO); hold on
plot(tm_full,smooth_bi_dioptic_NPO); hold on
vline(0); vline(800)
legend('binoc dioptic PO','binoc dioptic NPO')
title('Binocular onset - dioptic')
xlim([-200 1600])
box off

%% binocular dichoptic
for trl = 1:length(IDX(3).correctTrialIndex)
    array_bi_dichoptic_1(tm1,:,trl) = IDX(3).MUAe{trl,1}(tm1,ch); 
    array_bi_dichoptic_1(tm2_concat,:,trl) = IDX(3).MUAe{trl,2}(tm2,ch); 
    avg_bi_dichoptic_1 = median(median(array_bi_dichoptic_1,3),2);
end
for trl = 1:length(IDX(4).correctTrialIndex)
    array_bi_dichoptic_2(tm1,:,trl) = IDX(4).MUAe{trl,1}(tm1,ch); 
    array_bi_dichoptic_2(tm2_concat,:,trl) = IDX(4).MUAe{trl,2}(tm2,ch); 
    avg_bi_dichoptic_2 = median(median(array_bi_dichoptic_2,3),2);
end

% Baseline subtraction
bl_bi_dichoptic_1 = avg_bi_dichoptic_1 - median(avg_bi_dichoptic_1(100:200,1));
bl_bi_dichoptic_2 = avg_bi_dichoptic_2 - median(avg_bi_dichoptic_2(100:200,1));

% Smooth data
smooth_bi_dichoptic_1 = smoothdata(bl_bi_dichoptic_1,"gaussian",20);
smooth_bi_dichoptic_2 = smoothdata(bl_bi_dichoptic_2,"gaussian",20);

% plot
nexttile
plot(tm_full,smooth_bi_dichoptic_1); hold on
plot(tm_full,smooth_bi_dichoptic_2); hold on
vline(0); vline(800)
legend('binoc dichoptic 1','binoc dichoptic 2')
title('Binocular onset - dichoptic')
xlim([-200 1600])
box off

%% BRFS-like dioptic
% Ok, you have a good formula. Now lets try to do this progromatically
nexttile(5)
count = 0;
for cond = [5,6,7,8]
    count = count + 1;
    % convert cell to array
    clear array_diopticAdapted
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_diopticAdapted(tm1,:,trl) = IDX(cond).MUAe{trl,1}(tm1,ch); % now we index the cell for the second 800ms
        array_diopticAdapted(tm2_concat,:,trl) = IDX(cond).MUAe{trl,2}(tm2,ch); % now we index the cell for the second 800ms
    end
    % and average across whole electrode
    avg_diopticAdapted(:,count) = median(median(array_diopticAdapted,3),2); % average across trials and average across electrode
    bl_diopticAdapted(:,count) = avg_diopticAdapted(:,count) - median(avg_diopticAdapted(100:200,count));
    smooth_diopticAdapted(:,count) = smoothdata(bl_diopticAdapted(:,count),"gaussian",20);
    plot(tm_full,smooth_diopticAdapted(:,count)); hold on
end
vline(0); vline(800)
legend('PO RightEye -> PO LeftEye','NPO LeftEye -> NPO RightEye','NPO RightEye -> NPO LeftEye','PO LeftEye -> PO RightEye')
title('BRFS-like adaptation - dioptic')
xlim([-200 1600])
box off

%% BRFS (dichoptic)
nexttile(6)
count = 0;
for cond = [9,10,11,12]
    count = count + 1;
    % convert cell to array
    clear array_dichopticAdapted
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_dichopticAdapted(tm1,:,trl) = IDX(cond).MUAe{trl,1}(tm1,ch); % now we index the cell for the second 800ms
        array_dichopticAdapted(tm2_concat,:,trl) = IDX(cond).MUAe{trl,2}(tm2,ch); % now we index the cell for the second 800ms
    end
    % and average across whole electrode
    avg_dichopticAdapted(:,count) = median(median(array_dichopticAdapted,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    bl_dichopticAdapted(:,count) = avg_dichopticAdapted(:,count) - median(avg_dichopticAdapted(100:200,count));
    smooth_dichopticAdapted(:,count) = smoothdata(bl_dichopticAdapted(:,count),"gaussian",20);
    plot(tm_full,smooth_dichopticAdapted(:,count)); hold on
end
vline(0); vline(800)
legend('NPO RightEye -> PO LeftEye','PO LeftEye -> NPO RightEye','PO RightEye -> NPO LeftEye','NPO LeftEye -> PO RightEye')
title('BRFS')
xlim([-200 1600])
box off
%% Monocular alternation dioptic
nexttile(8)
count = 0;
for cond = [13,14,15,16]
    count = count + 1;
    % convert cell to array
    clear array_diopticMonocAlt
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_diopticMonocAlt(tm1,:,trl) = IDX(cond).MUAe{trl,1}(tm1,ch); % now we index the cell for the second 800ms
        array_diopticMonocAlt(tm2_concat,:,trl) = IDX(cond).MUAe{trl,2}(tm2,ch); % now we index the cell for the second 800ms
    end
    % and average across whole electrode
    avg_diopticMonocAlt(:,count) = median(median(array_diopticMonocAlt,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    bl_diopticMonocAlt(:,count) = avg_diopticMonocAlt(:,count) - median(avg_diopticMonocAlt(100:200,count));
    smooth_diopticMonocAlt(:,count) = smoothdata(bl_diopticMonocAlt(:,count),"gaussian",20);
    plot(tm_full,smooth_diopticMonocAlt(:,count)); hold on
end
vline(0); vline(800);
legend('PO RightEye -> PO LeftEye','NPO LeftEye -> NPO RightEye','NPO RightEye -> NPO LeftEye','PO LeftEye -> PO RightEye')
title('Monocular alternation - dioptic')
xlim([-200 1600])
box off

%% Monocular alternation dichoptic
nexttile(9)
count = 0;
for cond = [17,18,19,20]
    count = count + 1;
    % convert cell to array
    clear array_dichopticMonocAlt
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_dichopticMonocAlt(tm1,:,trl) = IDX(cond).MUAe{trl,1}(tm1,ch); % now we index the cell for the second 800ms
        array_dichopticMonocAlt(tm2_concat,:,trl) = IDX(cond).MUAe{trl,2}(tm2,ch); % now we index the cell for the second 800ms
    end
    % and average across whole electrode
    avg_dichopticMonocAlt(:,count) = median(median(array_dichopticMonocAlt,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    bl_dichopticMonocAlt(:,count) = avg_dichopticMonocAlt(:,count) - median(avg_dichopticMonocAlt(100:200,count));
    smooth_dichopticMonocAlt(:,count) = smoothdata(bl_dichopticMonocAlt(:,count),"gaussian",20);
    plot(tm_full,smooth_dichopticMonocAlt(:,count)); hold on
end
vline(0); vline(800);
legend('NPO RightEye -> PO LeftEye','PO LeftEye -> NPO RightEye','PO RightEye -> NPO LeftEye','NPO LeftEye -> PO RightEye')
title('Monocular alternation - dichoptic')
xlim([-200 1600])
box off

%% Figure title and axis
title(t,sessionLabel,'Interpreter','none')
xlabel(t,'Time (ms)')
ylabel(t,'Neural reponse (uV)')

