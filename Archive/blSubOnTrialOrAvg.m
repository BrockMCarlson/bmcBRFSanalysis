% Plot average of all monocular conditions
% show how taking the basleine at the trial level before averaging results
% in plots with a baseline below zero.

%% Setup
clear
% sessionLabel = '211008_B_bmcBRFS001';
% sessionLabel = '221206_J_bmcBRFS003';
sessionLabel = '221102_J_bmcBRFS001'; 

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

% convert from cell to double for single condition
count = 0;
for cond = monoc_1
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond1(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); % 1000 x 32
        avg_ofCond1(:,count) = median(median(array_ofCond1,3),2); 
    end
end
avg_ofMonoc1 = median(avg_ofCond1,2);
clear cond count trl 

count = 0;
for cond = monoc_2
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond2(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); 
        avg_ofCond2(:,count) = median(median(array_ofCond2,3),2); 
    end
end
avg_ofMonoc2 = median(avg_ofCond2,2);
clear cond count trl 


count = 0;
for cond = monoc_3
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond3(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); 
        avg_ofCond3(:,count) = median(median(array_ofCond3,3),2);
    end
end
avg_ofMonoc3 = median(avg_ofCond3,2);
clear cond count trl 


count = 0;
for cond = monoc_4
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        array_ofCond4(:,:,trl) = IDX(cond).MUAe{trl,1}(:,ch); 
        avg_ofCond4(:,count) = median(median(array_ofCond4,3),2);
    end
end
avg_ofMonoc4 = median(avg_ofCond4,2);
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
smooth_ofMonoc1 = smoothdata(blSub1,"gaussian");
smooth_ofMonoc2 = smoothdata(blSub2,"gaussian");
smooth_ofMonoc3 = smoothdata(blSub3,"gaussian");
smooth_ofMonoc4 = smoothdata(blSub4,"gaussian");

% Plot
figure
set(gcf,"Position",[17.6667 59 2.5313e+03 1.2813e+03])
t = tiledlayout(1,2);
set(0, 'DefaultLineLineWidth', 1.5);
nexttile
tm = -200:800;
plot(tm,smooth_ofMonoc1); hold on
plot(tm,smooth_ofMonoc2); hold on
plot(tm,smooth_ofMonoc3); hold on
plot(tm,smooth_ofMonoc4); hold on
ylim([-0.3 1.3])
vline(0)
legend('PO RightEye','PO LeftEye','NPO RightEye','NPO LeftEye')
title('Baseline subtraction after taking median of trials')
box off

%% bl Sub at trial level (needed for ANOVA)
% Monocular
monoc_1 = [5, 11, 13, 19]; % PO RightEye
monoc_2 = [8, 10, 16, 18]; % PO LeftEye
monoc_3 = [7, 9, 15, 17];  % NPO RightEye
monoc_4 = [6, 12, 14, 20]; % NPO LeftEye

% convert from cell to double for single condition
count = 0;
for cond = monoc_1
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        rawResponse = IDX(cond).MUAe{trl,1}(:,ch); % 1000 x 32
        baseLine = mean(rawResponse(100:200,:),1); % 1 x 32
        array_ofCond1(:,:,trl) = rawResponse - baseLine; % 1000 x 32 - 1 x 32, result should be 1000 x 32
        avg_ofCond1(:,count) = median(median(array_ofCond1,3),2); 
    end
end
avg_ofMonoc1 = median(avg_ofCond1,2);
clear cond count trl 

count = 0;
for cond = monoc_2
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        rawResponse = IDX(cond).MUAe{trl,1}(:,ch); 
        baseLine = mean(rawResponse(100:200,:),1); % 1 x 32
        array_ofCond2(:,:,trl) = rawResponse - baseLine; 
        avg_ofCond2(:,count) = median(median(array_ofCond2,3),2); 
    end
end
avg_ofMonoc2 = median(avg_ofCond2,2);
clear cond count trl 


count = 0;
for cond = monoc_3
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        rawResponse = IDX(cond).MUAe{trl,1}(:,ch); 
        baseLine = mean(rawResponse(100:200,:),1); % 1 x 32
        array_ofCond3(:,:,trl) = rawResponse - baseLine;
        avg_ofCond3(:,count) = median(median(array_ofCond3,3),2);
    end
end
avg_ofMonoc3 = median(avg_ofCond3,2);
clear cond count trl 


count = 0;
for cond = monoc_4
    count = count + 1;
    for trl = 1:length(IDX(cond).correctTrialIndex)
        rawResponse = IDX(cond).MUAe{trl,1}(:,ch); 
        baseLine = mean(rawResponse(100:200,:),1); % 1 x 32
        array_ofCond4(:,:,trl) = rawResponse - baseLine;
        avg_ofCond4(:,count) = median(median(array_ofCond4,3),2);
    end
end
avg_ofMonoc4 = median(avg_ofCond4,2);
clear cond count trl 

% % % Baseline subtractions
% % bl_1 = median(avg_ofMonoc1(100:200,1));
% % bl_2 = median(avg_ofMonoc2(100:200,1));
% % bl_3 = median(avg_ofMonoc3(100:200,1));
% % bl_4 = median(avg_ofMonoc4(100:200,1));
% % blSub1 = avg_ofMonoc1 - bl_1;
% % blSub2 = avg_ofMonoc2 - bl_2;
% % blSub3 = avg_ofMonoc3 - bl_3;
% % blSub4 = avg_ofMonoc4 - bl_4;

% Smooth data
smooth_ofMonoc1 = smoothdata(avg_ofMonoc1,"gaussian");
smooth_ofMonoc2 = smoothdata(avg_ofMonoc2,"gaussian");
smooth_ofMonoc3 = smoothdata(avg_ofMonoc3,"gaussian");
smooth_ofMonoc4 = smoothdata(avg_ofMonoc4,"gaussian");

% Plot
nexttile
tm = -200:800;
plot(tm,smooth_ofMonoc1); hold on
plot(tm,smooth_ofMonoc2); hold on
plot(tm,smooth_ofMonoc3); hold on
plot(tm,smooth_ofMonoc4); hold on
ylim([-0.3 1.3])
vline(0)
legend('PO RightEye','PO LeftEye','NPO RightEye','NPO LeftEye')
title('Baseline subtraction on trial level')
box off


