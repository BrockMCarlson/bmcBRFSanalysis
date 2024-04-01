%% Laminar - all electrodes
datetime

%% Setup
clear
% Directories
path1 = strcat(getenv('HOMEDRIVE'),getenv("HOMEPATH"));
path2 = 'Documents\GitHub\bmcBRFSanalysis';
codeDir = strcat(path1,filesep,path2);
cd(codeDir)
path3 = 'Documents\MATLAB\formattedDataOutputs\figures_231201';

%% For loop
if contains(path1,'Brock Carlson')
    dataDir = 'S:\bmcBRFS_sortedData_Nov23';
    plotDir = 'S:\formattedDataOutputs';
elseif contains(path1,'Neuropoixel')
    dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
    plotDir = strcat(path1,filesep,path3);
else
    error('wrong Home Dir')
end
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
    
        % bl Sub at average level (better for plotting)
        ch = 1:size(IDX(1).(signalType){1,1},2);
        if ch > 32
            ch = 1:32;
            warning(['Congrats! You made it to the 2 electrode days...' ...
                ' Time to adjust your code to include the 2nd probe!'])
        end
    
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
                array_ofMonoc1(:,:,count) = IDX(cond).(signalType){trl,1}(:,ch); % 1000 x 32
            end
        end
        clear cond count trl 
        count = 0;
        for cond = monoc_2
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc2(:,:,count) = IDX(cond).(signalType){trl,1}(:,ch); 
            end
        end
        clear cond count trl     
        count = 0;
        for cond = monoc_3
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc3(:,:,count) = IDX(cond).(signalType){trl,1}(:,ch); 
            end
        end
        clear cond count trl     
        count = 0;
        for cond = monoc_4
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc4(:,:,count) = IDX(cond).(signalType){trl,1}(:,ch); 
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
        
        % % % Skip this file if no tuned units are found
        % % DATAOUT(file).numberOFUnits = sum(tuned);
        % % if sum(tuned) == 0
        % %     warning(strcat('No tuned channels on',sessionLabel,signalType))
        % %     continue
        % % end
    
    
    
        % Create array of preference-based data
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
                array_dichopticAdapted_pref(tm1,i,trl) = IDX(prefCondOnFlash).(signalType){trl,1}(tm1,i); % now we index the cell for the second 800ms
                array_dichopticAdapted_pref(tm2_concat,i,trl) = IDX(prefCondOnFlash).(signalType){trl,2}(tm2,i); % now we index the cell for the second 800ms
            end
            for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
                array_dichopticAdapted_null(tm1,i,trl) = IDX(nullCondOnFlash).(signalType){trl,1}(tm1,i); % now we index the cell for the second 800ms
                array_dichopticAdapted_null(tm2_concat,i,trl) = IDX(nullCondOnFlash).(signalType){trl,2}(tm2,i); % now we index the cell for the second 800ms
            end
        end
        
    
        
        %% Create stackedLinePlot
        % Create matrix of all trials (time x ch x trial)
        % % IDX(CONDITION).(signalType){TRIAL,TRIAL PERIOD}(TIME,CHANNEL)
        % % psTrlLength = size(IDX(idxps).correctTrialIndex,1);
        % % for psTrl = 1:psTrlLength
        % %     ps_msXchXtrl(:,:,psTrl) = IDX(idxps).(signalType){psTrl,2}(:,v1Ch); % MUA output is time x ch
        % % end
        % % nsTrlLength = size(IDX(idxns).correctTrialIndex,1);
        % % for nsTrl = 1:nsTrlLength
        % %     ns_msXchXtrl(:,:,nsTrl) = IDX(idxns).(signalType){nsTrl,2}(:,v1Ch); % MUA output is time x ch
        % % end
        %trl avg
        ps_avg = median(array_dichopticAdapted_pref,3,"omitmissing"); % input is (tm,ch,trl)
        ns_avg = median(array_dichopticAdapted_null,3,"omitmissing");
    
        % Calculate as Percent Change
        %              X(t) - avgBl
        % %Ch = 100 * -------------
        %                 avgBl
    
    
    
        % smooth data
        ps_smooth = smoothdata(ps_avg',1,"gaussian",20);
        ns_smooth = smoothdata(ns_avg',1,"gaussian",20);
    
        % convert to table
        ps_table_1 = array2table(ps_smooth(:,1:16));
        ns_table_1 = array2table(ns_smooth(:,1:16));
        ps_table_2 = array2table(ps_smooth(:,17:32));
        ns_table_2 = array2table(ns_smooth(:,17:32));
        
        %Now convert ps_avg (a double array) into a table (input to timetable must
        %be a table)
        % The goal is to have 32 variables, each as a column, representing a
        % different depth. 
        Time = milliseconds(-200):milliseconds(1):milliseconds(800);
        % electrodes 1-16
        elName_1 = {'el1','el2','el3','el4','el5','el6','el7','el8','el9',...
            'el10','el11','el12','el13','el14','el15','el16'};
        ps_TT_1 = table2timetable(ps_table_1,'RowTimes',Time);
        preferredStimFlash_1 = renamevars(ps_TT_1,ps_TT_1.Properties.VariableNames,elName_1);
        ns_TT_1 = table2timetable(ns_table_1,'RowTimes',Time);
        nullStimFlash_1 = renamevars(ns_TT_1,ns_TT_1.Properties.VariableNames,elName_1);
        % electrodes 17-32
        elName_2 = {'el17','el18','el19','el20','el21','el22','el23','el24',...
            'el25','el26','el27','el28','el29','el30','el31','el32'};
        ps_TT_2 = table2timetable(ps_table_2,'RowTimes',Time);
        preferredStimFlash_2 = renamevars(ps_TT_2,ps_TT_2.Properties.VariableNames,elName_2);
        ns_TT_2 = table2timetable(ns_table_2,'RowTimes',Time);
        nullStimFlash_2 = renamevars(ns_TT_2,ns_TT_2.Properties.VariableNames,elName_2);    
        
        %% stackedplot()
        close all
        f = figure;
        set(f,"Position",[1000 54.3333 1213 1.1833e+03])
        tl = tiledlayout(3,2);
        
        nexttile(tl,[3 1])
        s_1 = stackedplot(preferredStimFlash_1,nullStimFlash_1);
        s_1.LineWidth = 1;
    
        nexttile(tl,[3 1])
        s_2 = stackedplot(preferredStimFlash_2,nullStimFlash_2);
        s_2.LineWidth = 1;
    
        titleText = {sessionLabel,signalType};
        title(tl,titleText,'Interpreter','none')
    
        
        %save fig
        cd(plotDir)
        figName = strcat('stackedPlot_',sessionLabel,'_',signalType,'.png');
        saveas(f,figName)
    
    
    end
    disp(strcat('completed_',sessionLabel))
end
datetime