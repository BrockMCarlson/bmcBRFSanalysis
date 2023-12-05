%% Laminar_grandAvg
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur


datetime

%% Setup
clear
% Directories
path1 = strcat(getenv('HOMEDRIVE'),getenv("HOMEPATH"));
path2 = 'Documents\GitHub\bmcBRFSanalysis';
codeDir = strcat(path1,filesep,path2);
cd(codeDir)
path3 = 'Documents\MATLAB\formattedDataOutputs\figures_231204';
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

    granBtm = officLamAssign.Probe11stFold4c(file); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',sessionLabel))
        continue
    end
    
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    columnNames = {'sg_1','sg_2','sg_3','sg_4','sg_5','g_1','g_2','g_3','g_4','g_5','ig_1','ig_2','ig_3','ig_4','ig_5'};
    if any(v1Ch > 32) || any(v1Ch < 1)
        warning('skipping session without full column for now')
        disp(sessionLabel)
        continue
        % % columnNames = columnNames((v1Ch >= 1) & (v1Ch<=32));
        % % v1Ch        = v1Ch((v1Ch >= 1) & (v1Ch<=32));
    end
    
    ch = v1Ch;

    for st = 1:length(signalTypeList)
        signalType = signalTypeList{st};
    
        %% bl Sub at average level (better for plotting)
        clear array_ofMonoc1 array_ofMonoc2 array_ofMonoc3 array_ofMonoc4
        
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
                array_ofMonoc1(:,:,count) = abs(IDX(cond).(signalType){trl,1}(:,ch)); % 1000 x 32
            end
        end
        clear cond count trl 
        
        count = 0;
        for cond = monoc_2
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc2(:,:,count) = abs(IDX(cond).(signalType){trl,1}(:,ch)); 
            end
        end
        clear cond count trl 
        
        
        count = 0;
        for cond = monoc_3
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc3(:,:,count) = abs(IDX(cond).(signalType){trl,1}(:,ch)); 
            end
        end
        clear cond count trl 
        
        
        count = 0;
        for cond = monoc_4
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc4(:,:,count) = abs(IDX(cond).(signalType){trl,1}(:,ch)); 
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
        for i = 1:length(ch)
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
                array_dichopticAdapted_pref(tm1,i,trl) = abs(IDX(prefCondOnFlash).(signalType){trl,1}(tm1,i)); % now we index the cell for the second 800ms
                array_dichopticAdapted_pref(tm2_concat,i,trl) = abs(IDX(prefCondOnFlash).(signalType){trl,2}(tm2,i)); % now we index the cell for the second 800ms
            end
            for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
                array_dichopticAdapted_null(tm1,i,trl) = abs(IDX(nullCondOnFlash).(signalType){trl,1}(tm1,i)); % now we index the cell for the second 800ms
                array_dichopticAdapted_null(tm2_concat,i,trl) = abs(IDX(nullCondOnFlash).(signalType){trl,2}(tm2,i)); % now we index the cell for the second 800ms
            end
        end
        
        
        %% Save data into structure array
      
        ps_avg = median(array_dichopticAdapted_pref,3,"omitmissing"); % input is (tm,ch,trl)
        ns_avg = median(array_dichopticAdapted_null,3,"omitmissing");
    
    
        % Calculate as Percent Change
        %              X(t) - avgBl
        % %Ch = 100 * -------------
        %                 avgBl
        psBl = median(ps_avg(100:200,:));
        ps_PercentC = 100*((ps_avg-psBl)./psBl);
        nsBl = median(ns_avg(100:200,:));
        ns_PercentC = 100*((ns_avg-nsBl)./nsBl);  
        
        % Save this output
        DATAOUT.(signalType).ps_PercentC(:,:,file) = ps_PercentC;
        DATAOUT.(signalType).ns_PercentC(:,:,file) = ns_PercentC;
       
      
        
    end
         disp(strcat('finished processing _',sessionLabel))

end

datetime

%% plot grand averages of whole probe
for st = 1:length(signalTypeList)
    signalType = signalTypeList{st};

    % convert 0 to NaN
    ps_NaNmatrix = DATAOUT.(signalType).ps_PercentC;
    ps_NaNmatrix(ps_NaNmatrix==0) = NaN;
    ns_NaNmatrix = DATAOUT.(signalType).ns_PercentC;
    ns_NaNmatrix(ns_NaNmatrix==0) = NaN;

    ps_grandAvg = mean(ps_NaNmatrix,3,"omitmissing");
    ns_grandAvg = mean(ns_NaNmatrix,3,"omitmissing");
    
    % smooth data
    ps_smooth_grandAvg = smoothdata(ps_grandAvg,1,"gaussian",20);
    ns_smooth_grandAvg = smoothdata(ns_grandAvg,1,"gaussian",20);
    
    % convert to table
    ps_table_grandAvg = array2table(ps_smooth_grandAvg);
    ns_table_grandAvg = array2table(ns_smooth_grandAvg);
    
    
    %Now convert ps_avg (a double array) into a table (input to timetable must
    %be a table)
    % The goal is to have 32 variables, each as a column, representing a
    % different depth. 
    Time = milliseconds(-200):milliseconds(1):milliseconds(1600);
    ps_TT_grandAvg = table2timetable(ps_table_grandAvg,'RowTimes',Time);
    preferredStimFlash_grandAvg = renamevars(ps_TT_grandAvg,ps_TT_grandAvg.Properties.VariableNames,columnNames);
    ns_TT_grandAvg = table2timetable(ns_table_grandAvg,'RowTimes',Time);
    nullStimFlash_grandAvg = renamevars(ns_TT_grandAvg,ns_TT_grandAvg.Properties.VariableNames,columnNames);
    
    
    % stackedplot()
    close all
    stk_grandAvg = figure;
    set(stk_grandAvg,"Position",[1000 60.3333 560 1.2933e+03])
    s_grandAvg = stackedplot(preferredStimFlash_grandAvg,nullStimFlash_grandAvg);
    titleText_grandAvg = {'_grandAvg',signalType};
    s_grandAvg.Title = titleText_grandAvg;
    s_grandAvg.LineWidth = 1;
    % s_grandAvg.DisplayLabels = ["% Change"];
    
    %save fig
    cd(plotDir)
    figName_grandAvg = strcat('stackedPlot_','_grandAvg_',signalType,'.png');
    saveas(stk_grandAvg,figName_grandAvg)
end

%% plot grand averages of laminar compartments

for st = 1:length(signalTypeList)
    signalType = signalTypeList{st};

    % convert 0 to NaN
    ps_NaNmatrix = DATAOUT.(signalType).ps_PercentC;
    ps_NaNmatrix(ps_NaNmatrix==0) = NaN;
    ns_NaNmatrix = DATAOUT.(signalType).ns_PercentC;
    ns_NaNmatrix(ns_NaNmatrix==0) = NaN;

    ps_grandAvg = mean(ps_NaNmatrix,3,"omitmissing");
    ns_grandAvg = mean(ns_NaNmatrix,3,"omitmissing");

    %Now convert grand average array into laminar compartments
    ps_S = mean(ps_grandAvg(:,1:5),2); 
    ps_G = mean(ps_grandAvg(:,6:10),2); 
    ps_I = mean(ps_grandAvg(:,11:15),2); 
    ns_S = mean(ns_grandAvg(:,1:5),2); 
    ns_G = mean(ns_grandAvg(:,6:10),2); 
    ns_I = mean(ns_grandAvg(:,11:15),2); 

    % smooth data
    ps_S_smooth = smoothdata(ps_S,1,"gaussian",20);
    ps_G_smooth = smoothdata(ps_G,1,"gaussian",20);
    ps_I_smooth = smoothdata(ps_I,1,"gaussian",20);
    ns_S_smooth = smoothdata(ns_S,1,"gaussian",20);
    ns_G_smooth = smoothdata(ns_G,1,"gaussian",20);
    ns_I_smooth = smoothdata(ns_I,1,"gaussian",20);
    
    % Ok, the data is together for plotting, now lets run statistics on
    % each laminar compartment to see if perceptual modulation occurs. The
    % goal here is to run a t-test to see if the average response between
    % 1200 and 1600ms significantly differs between ps and ns
    useIdx = squeeze(~isnan(ps_NaNmatrix(1,1,:))); 
    tInput_ps_S = reshape(squeeze(mean(ps_NaNmatrix(1400:1801,1:5,useIdx),1)),[],1);
    tInput_ps_G = reshape(squeeze(mean(ps_NaNmatrix(1400:1801,6:10,useIdx),1)),[],1);
    tInput_ps_I = reshape(squeeze(mean(ps_NaNmatrix(1400:1801,11:15,useIdx),1)),[],1);
    tInput_ns_S = reshape(squeeze(mean(ns_NaNmatrix(1400:1801,1:5,useIdx),1)),[],1);
    tInput_ns_G = reshape(squeeze(mean(ns_NaNmatrix(1400:1801,6:10,useIdx),1)),[],1);
    tInput_ns_I = reshape(squeeze(mean(ns_NaNmatrix(1400:1801,11:15,useIdx),1)),[],1);

    [h_S,p_S] = ttest2(tInput_ps_S,tInput_ns_S);
    [h_G,p_G] = ttest2(tInput_ps_G,tInput_ns_G);
    [h_I,p_I] = ttest2(tInput_ps_I,tInput_ns_I);

    % tiledLayout plot
    close all
    tm_full = -200:1600; % 1801 total timepoints
    lamCom = figure;
    set(gcf,"Position",[1000 503 560 734.6667])
    t = tiledlayout(3,1);
    nexttile
        plot(tm_full,ps_S_smooth,tm_full,ns_S_smooth,'LineWidth',1.5)
        vline(0)
        vline(800)
        xregion(1200,1600)
        ylabel({'Supragranular','% Change'})
        title({strcat('Significant perceptual modulation?_',string(h_S)),strcat('pVal =',string(p_S))},'Interpreter','none')
    nexttile
        plot(tm_full,ps_G_smooth,tm_full,ns_G_smooth,'LineWidth',1.5)
        vline(0)
        vline(800)
        xregion(1200,1600)
        ylabel({'Granular','% Change'})
        title({strcat('Significant perceptual modulation?_',string(h_G)),strcat('pVal =',string(p_G))},'Interpreter','none')
    nexttile
        plot(tm_full,ps_I_smooth,tm_full,ns_I_smooth,'LineWidth',1.5)
        vline(0)
        vline(800)
        xregion(1200,1600)
        ylabel({'Supragranular','% Change'})
        title({strcat('Significant perceptual modulation?_',string(h_I)),strcat('pVal =',string(p_I))},'Interpreter','none')

    titleText = {'Grand average of 135 multi-units per laminar compartment',signalType};
    title(t,titleText,'Interpreter','none')

    %save fig
    cd(plotDir)
    figName_lamCom = strcat('laminarCompartment_','_grandAvg_',signalType,'.png');
    saveas(lamCom,figName_lamCom)
end

