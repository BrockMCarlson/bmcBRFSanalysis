%% Laminar - stimulus type comparison


datetime

%% Setup
clear
% Directories
path1 = strcat(getenv('HOMEDRIVE'),getenv("HOMEPATH"));
path2 = 'Documents\GitHub\bmcBRFSanalysis';
codeDir = strcat(path1,filesep,path2);
cd(codeDir)
path3 = 'Documents\MATLAB\formattedDataOutputs\figures_231201';
plotDir = strcat(path1,filesep,path3);

%Import laminar assignments
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "Sheet1", [2, Inf]);

%% For loop
dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(dataDir)
allDataFiles = dir('**/*sortedData*.mat');

% signalTypeList = {'LFP_bb','LFP_delta','LFP_theta','LFP_alpha','LFP_beta1',...
%     'LFP_beta2','LFP_beta3','LFP_gamma1','LFP_gamma2',...
%     'CSD_bb','CSD_delta','CSD_theta','CSD_alpha','CSD_beta1','CSD_beta2',...
%     'CSD_beta3','CSD_gamma1','CSD_gamma2','MUAe'};
signalTypeList = {'LFP_bb','LFP_alpha','LFP_beta1',...
    'LFP_beta2','LFP_beta3','LFP_gamma1','LFP_gamma2',...
    'CSD_bb','CSD_alpha','CSD_beta1','CSD_beta2',...
    'CSD_beta3','CSD_gamma1','CSD_gamma2','MUAe'};

for file = 1:length(allDataFiles)
    
    % load data
    cd(dataDir)
    fileToLoad = allDataFiles(file).name;
    load(fileToLoad)
    sessionLabel = allDataFiles(file).name(12:end-4);

    for st = 1:length(signalTypeList)
        signalType = signalTypeList{st};
    
    granBtm = officLamAssign.Probe11stFold4c(file); % channel corresponding to the bottom of layer 4c
    
    
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
    if any(v1Ch > 32) || any(v1Ch < 1)
        columnNames = columnNames((v1Ch >= 1) & (v1Ch<=32));
        v1Ch        = v1Ch((v1Ch >= 1) & (v1Ch<=32));
    end
    
    
    %% test for monocular preference
    ch = v1Ch;
    clear array_ofMonoc1 array_ofMonoc2 array_ofMonoc3 array_ofMonoc4
    
    % Monocular
    monoc_1 = [5, 11, 13, 19]; % PO RightEye
    monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    
    % convert from cell to double and combine monocular conditions
    warning('taking abs() of IDX value')
    count = 0;
    for cond = monoc_1
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(IDX(cond).CSD_bb{trl,1}(:,ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(IDX(cond).CSD_bb{trl,1}(:,ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(IDX(cond).CSD_bb{trl,1}(:,ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:length(IDX(cond).correctTrialIndex)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(IDX(cond).CSD_bb{trl,1}(:,ch)); 
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
    tuned = nan(length(ch),1);
    prefMonoc = nan(length(ch),1);
    
    % FOR each individual unit 
    for i = 1:length(ch)
        % preallocate y based on trial count
        y = nan(maxTrls,4);
        % create an array of each trial's median response for each monoc 
        % condition (trl x 4 monoc)
        y(1:size(array_ofMonoc1,3),1) = median(squeeze(array_ofMonoc1(200:800,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc2,3),2) = median(squeeze(array_ofMonoc2(200:800,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc3,3),3) = median(squeeze(array_ofMonoc3(200:800,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc4,3),4) = median(squeeze(array_ofMonoc4(200:800,i,:)),1)'; % median of each trial after stim onset
        % Now we perform baseline subtraction
        % First we get the blAvg for this contact across all trials
        y_bl(1,1) = median(median(squeeze(array_ofMonoc1(100:200,i,:)),1)');
        y_bl(1,2) = median(median(squeeze(array_ofMonoc2(100:200,i,:)),1)');
        y_bl(1,3) = median(median(squeeze(array_ofMonoc3(100:200,i,:)),1)');
        y_bl(1,4) = median(median(squeeze(array_ofMonoc4(100:200,i,:)),1)');
        y_blSub = abs(y - median(y_bl));
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
    
    % Skip this file if no tuned units are found
    % % DATAOUT(file).numberOFUnits = sum(tuned);
    % % if sum(tuned) == 0
    % %     warning(strcat('No tuned channels on',sessionLabel))
    % %     continue
    % % end
    % % 
    % % % Creat channel array to use in subsequent steps - only plot tuned units.
    % % chTuned = ch(logical(tuned));
    % % 
    % % 


    %% Create array of preference-based data
    % concatenate two timecourses into array
    tm_full = -200:1600; % 1801 total timepoints
    tm1 = 1:801;
    tm2 = 1:1001;
    tm2_concat = 801:1801;

    % pre allocate
    maxTrlLength = max([length(IDX(9).correctTrialIndex),...
        length(IDX(10).correctTrialIndex),...
        length(IDX(11).correctTrialIndex),...
        length(IDX(12).correctTrialIndex)]);
    array_dichopticAdapted_pref = nan(1801,length(ch),maxTrlLength);
    array_dichopticAdapted_null = nan(1801,length(ch),maxTrlLength);
    for i = 1:length(ch)
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
            array_dichopticAdapted_pref(tm1,i,trl) = abs(IDX(prefCondOnFlash).CSD_bb{trl,1}(tm1,i)); % now we index the cell for the second 800ms
            array_dichopticAdapted_pref(tm2_concat,i,trl) = abs(IDX(prefCondOnFlash).CSD_bb{trl,2}(tm2,i)); % now we index the cell for the second 800ms
        end
        for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
            array_dichopticAdapted_null(tm1,i,trl) = abs(IDX(nullCondOnFlash).CSD_bb{trl,1}(tm1,i)); % now we index the cell for the second 800ms
            array_dichopticAdapted_null(tm2_concat,i,trl) = abs(IDX(nullCondOnFlash).CSD_bb{trl,2}(tm2,i)); % now we index the cell for the second 800ms
        end
    end
    
    % % % and average across whole electrode to see grand average
    % % avg_dichopticAdapted_pref = median(median(array_dichopticAdapted_pref,3,"omitmissing"),2,"omitmissing"); % average across trials and average across electrode
    % % bl_dichopticAdapted_pref(:,1) = avg_dichopticAdapted_pref(:,1) - median(avg_dichopticAdapted_pref(100:200,1));
    % % smooth_dichopticAdapted_pref(:,1) = smoothdata(bl_dichopticAdapted_pref(:,1),"gaussian",20);
    % % figure
    % % plot(tm_full,smooth_dichopticAdapted_pref); hold on
    % % avg_dichopticAdapted_null= median(median(array_dichopticAdapted_null,3,"omitmissing"),2,"omitmissing");% average across trials and average across electrode
    % % bl_dichopticAdapted_null(:,1) = avg_dichopticAdapted_null(:,1) - median(avg_dichopticAdapted_null(100:200,1));
    % % smooth_dichopticAdapted_null(:,1) = smoothdata(bl_dichopticAdapted_null(:,1),"gaussian",20);
    % % plot(tm_full,smooth_dichopticAdapted_null); hold on
    % % 
    % % vline(0); vline(800)
    % % legend('Preferred condition flashed','Same simulus, but null flashed')
    % % title('BRFS')
    % % xlim([-200 1600])
    % % box off
    
    %% Create stackedLinePlot
  
    ps_avg = median(array_dichopticAdapted_pref,3,"omitmissing"); % input is (tm,ch,trl)
    ns_avg = median(array_dichopticAdapted_null,3,"omitmissing");


    % Calculate as Percent Change
    %              X(t) - avgBl
    % %Ch = 100 * -------------
    %                 avgBl
    psBl = median(ps_avg(100:200,:));
    ps_PercentC = 100*((ps_avg-psBl)./psBl);
    DATAOUT(file).alternate_dichoptic_pref = bl_dichopticMonocAlt_pref;
    nsBl = median(ns_avg(100:200,:));
    ns_PercentC = 100*((ns_avg-nsBl)./nsBl);    

    % % % bl Subtract
    % % psBl = median(ps_avg(100:200,:));
    % % ps_blSub = abs(ps_avg - psBl);
    % % nsBl = median(ns_avg(100:200,:));
    % % ns_blSub = abs(ns_avg - nsBl);


    

    % smooth data
    ps_smooth = smoothdata(ps_PercentC,1,"gaussian",20);
    ns_smooth = smoothdata(ns_PercentC,1,"gaussian",20);

    % convert to table
    ps_table = array2table(ps_smooth);
    ns_table = array2table(ns_smooth);

    
    %Now convert ps_avg (a double array) into a table (input to timetable must
    %be a table)
    % The goal is to have 32 variables, each as a column, representing a
    % different depth. 
    Time = milliseconds(-200):milliseconds(1):milliseconds(1600);
    ps_TT = table2timetable(ps_table,'RowTimes',Time);
    preferredStimFlash = renamevars(ps_TT,ps_TT.Properties.VariableNames,columnNames);
    ns_TT = table2timetable(ns_table,'RowTimes',Time);
    nullStimFlash = renamevars(ns_TT,ns_TT.Properties.VariableNames,columnNames);
    
    
    %% stackedplot()
    close all
    stk = figure;
    set(stk,"Position",[1000 60.3333 560 1.2933e+03])
    s = stackedplot(preferredStimFlash,nullStimFlash);
    titleText = {sessionLabel,signalType};
    s.Title = titleText;
    s.LineWidth = 1;
    
    %save fig
    % cd(plotDir)
    % figName = strcat('stackedPlot_',sessionLabel,'_',signalType,'.png');
    % saveas(stk,figName)

    end

end