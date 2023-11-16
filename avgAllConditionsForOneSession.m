% Plot average of all conditions


%% Setup
sessionLabel = '211008_B_bmcBRFS001';

% Directories
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(outputDir)

% load data
load(strcat('sortedData_',sessionLabel,'.mat'))

%% start with monocular
close all

% Monocular
monoc_1 = [5, 11, 13, 19];
monoc_2 = [8, 10, 16, 18];
monoc_3 = [7, 9, 15, 17];
monoc_4 = [6, 12, 14, 20];

% convert from cell to double for single condition
count = 0;
for cond = monoc_1
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond(:,:,trl) = IDX(cond).MUAe{trl,1}; 
        avg_ofCond(:,count) = mean(mean(array_ofCond,3),2);
    end
end
avg_ofMonoc1 = mean(avg_ofCond,2);
clear cond count trl array_ofCond avg_ofCond

count = 0;
for cond = monoc_2
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond(:,:,trl) = IDX(cond).MUAe{trl,1}; 
        avg_ofCond(:,count) = mean(mean(array_ofCond,3),2);
    end
end
avg_ofMonoc2 = mean(avg_ofCond,2);
clear cond count trl array_ofCond avg_ofCond


count = 0;
for cond = monoc_3
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond(:,:,trl) = IDX(cond).MUAe{trl,1}; 
        avg_ofCond(:,count) = mean(mean(array_ofCond,3),2);
    end
end
avg_ofMonoc3 = mean(avg_ofCond,2);
clear cond count trl array_ofCond avg_ofCond


count = 0;
for cond = monoc_4
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond(:,:,trl) = IDX(cond).MUAe{trl,1}; 
        avg_ofCond(:,count) = mean(mean(array_ofCond,3),2);
    end
end
avg_ofMonoc4 = mean(avg_ofCond,2);
clear cond count trl array_ofCond avg_ofCond

% Smooth data
smooth_ofMonoc1 = smoothdata(avg_ofMonoc1);
smooth_ofMonoc2 = smoothdata(avg_ofMonoc2);
smooth_ofMonoc3 = smoothdata(avg_ofMonoc3);
smooth_ofMonoc4 = smoothdata(avg_ofMonoc4);

% Plot
figure
tiledlayout('flow')
nexttile
tm = -200:800;
plot(tm,smooth_ofMonoc1); hold on
plot(tm,smooth_ofMonoc2); hold on
plot(tm,smooth_ofMonoc3); hold on
plot(tm,smooth_ofMonoc4); hold on
legend('monoc1','monoc2','monoc3','monoc4')


%% binocular dioptic
for trl = 1:length(IDX(1).correctTrialIndex)
    array_bi_dioptic_PO(:,:,trl) = IDX(1).MUAe{trl,1}; 
    avg_bi_dioptic_PO = mean(mean(array_bi_dioptic_PO,3),2);
end
for trl = 1:length(IDX(2).correctTrialIndex)
    array_bi_dioptic_NPO(:,:,trl) = IDX(2).MUAe{trl,1}; 
    avg_bi_dioptic_NPO = mean(mean(array_bi_dioptic_NPO,3),2);
end

% Smooth data
smooth_bi_dioptic_PO = smoothdata(avg_bi_dioptic_PO);
smooth_bi_dioptic_NPO = smoothdata(avg_bi_dioptic_NPO);

% plot
nexttile
plot(tm,smooth_bi_dioptic_PO); hold on
plot(tm,smooth_bi_dioptic_NPO); hold on
legend('binoc dioptic PO','binoc dioptic NPO')

%% binocular dichoptic
for trl = 1:length(IDX(3).correctTrialIndex)
    array_bi_dichoptic_1(:,:,trl) = IDX(3).MUAe{trl,1}; 
    avg_bi_dichoptic_1 = mean(mean(array_bi_dichoptic_1,3),2);
end
for trl = 1:length(IDX(4).correctTrialIndex)
    array_bi_dichoptic_2(:,:,trl) = IDX(4).MUAe{trl,1}; 
    avg_bi_dichoptic_2 = mean(mean(array_bi_dichoptic_2,3),2);
end

% Smooth data
smooth_bi_dichoptic_1 = smoothdata(avg_bi_dichoptic_1);
smooth_bi_dichoptic_2 = smoothdata(avg_bi_dichoptic_2);

% plot
nexttile
plot(tm,smooth_bi_dichoptic_1); hold on
plot(tm,smooth_bi_dichoptic_2); hold on
legend('binoc dichoptic 1','binoc dichoptic 2')

%% BRFS-like dioptic
% Ok, you have a good formula. Now lets try to do this progromatically
nexttile
count = 0;
for cond = [5,6,7,8]
    count = count + 1;
    % convert cell to array
    for trl = 1:length(IDX(cond).correctTrialIndex)
        clear array_diopticAdapted
        array_diopticAdapted(:,:,trl) = IDX(cond).MUAe{trl,2}; % now we index the cell for the second 800ms
        % and average across whole electrode
        avg_diopticAdapted(:,count) = mean(mean(array_diopticAdapted,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    end
    smooth_diopticAdapted(:,count) = smoothdata(avg_diopticAdapted(:,count));
    plot(tm,smooth_diopticAdapted(:,count)); hold on
end
legend('diopticAdapted 1','diopticAdapted 2','diopticAdapted 3','diopticAdapted 4')

%% BRFS (dichoptic)
nexttile
count = 0;
for cond = [9,10,11,12]
    count = count + 1;
    % convert cell to array
    for trl = 1:length(IDX(cond).correctTrialIndex)
        clear array_dichopticAdapted
        array_dichopticAdapted(:,:,trl) = IDX(cond).MUAe{trl,2}; % now we index the cell for the second 800ms
        % and average across whole electrode
        avg_dichopticAdapted(:,count) = mean(mean(array_dichopticAdapted,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    end
    smooth_dichopticAdapted(:,count) = smoothdata(avg_dichopticAdapted(:,count));
    plot(tm,smooth_dichopticAdapted(:,count)); hold on
end
legend('dichopticAdapted 1','dichopticAdapted 2','dichopticAdapted 3','dichopticAdapted 4')

%% Monocular alternation dioptic
nexttile
count = 0;
for cond = [13,14,15,16]
    count = count + 1;
    % convert cell to array
    for trl = 1:length(IDX(cond).correctTrialIndex)
        clear array_diopticMonocAlt
        array_diopticMonocAlt(:,:,trl) = IDX(cond).MUAe{trl,2}; % now we index the cell for the second 800ms
        % and average across whole electrode
        avg_diopticMonocAlt(:,count) = mean(mean(array_diopticMonocAlt,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    end
    smooth_diopticMonocAlt(:,count) = smoothdata(avg_diopticMonocAlt(:,count));
    plot(tm,smooth_diopticMonocAlt(:,count)); hold on
end
legend('diopticMonocAlt 1','diopticMonocAlt 2','diopticMonocAlt 3','diopticMonocAlt 4')

%% Monocular alternation dichoptic
nexttile
count = 0;
for cond = [13,14,15,16]
    count = count + 1;
    % convert cell to array
    for trl = 1:length(IDX(cond).correctTrialIndex)
        clear array_dichopticMonocAlt
        array_dichopticMonocAlt(:,:,trl) = IDX(cond).MUAe{trl,2}; % now we index the cell for the second 800ms
        % and average across whole electrode
        avg_dichopticMonocAlt(:,count) = mean(mean(array_dichopticMonocAlt,3),2); % output is 1001x1, indexted into 1001 x 4 array for condition
    end
    smooth_dichopticMonocAlt(:,count) = smoothdata(avg_dichopticMonocAlt(:,count));
    plot(tm,smooth_dichopticMonocAlt(:,count)); hold on
end
legend('dichopticMonocAlt 1','dichopticMonocAlt 2','dichopticMonocAlt 3','dichopticMonocAlt 4')

