%% fig4
% The goal of this script is to average together data from each laminar
% compartment to see if differential perceptual modulations occur


datetime

%% Setup
clear
% Directories

codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
cd(codeDir)
outDir = 'S:\formattedDataOutputs';
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

%% For loop
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)



for file = 1:size(officLamAssign,1)
    
    % load data
    cd(dataDir)
    probeName = char(officLamAssign.SessionProbe(1,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');
    load(fileToLoad)



    
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
                array_ofMonoc1(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); % 1000 x 32
            end
        end
        clear cond count trl 
        
        count = 0;
        for cond = monoc_2
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc2(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); 
            end
        end
        clear cond count trl 
        
        
        count = 0;
        for cond = monoc_3
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc3(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); 
            end
        end
        clear cond count trl 
        
        
        count = 0;
        for cond = monoc_4
            for trl = 1:length(IDX(cond).correctTrialIndex)
                count = count + 1;
                array_ofMonoc4(:,:,count) = abs(IDX(cond).MUAe{trl,1}(:,ch)); 
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
                array_dichopticAdapted_pref(tm1,i,trl) = abs(IDX(prefCondOnFlash).MUAe{trl,1}(tm1,i)); % now we index the cell for the second 800ms
                array_dichopticAdapted_pref(tm2_concat,i,trl) = abs(IDX(prefCondOnFlash).MUAe{trl,2}(tm2,i)); % now we index the cell for the second 800ms
            end
            for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
                array_dichopticAdapted_null(tm1,i,trl) = abs(IDX(nullCondOnFlash).MUAe{trl,1}(tm1,i)); % now we index the cell for the second 800ms
                array_dichopticAdapted_null(tm2_concat,i,trl) = abs(IDX(nullCondOnFlash).MUAe{trl,2}(tm2,i)); % now we index the cell for the second 800ms
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
        DATAOUT.MUAe.ps_PercentC(:,:,file) = ps_PercentC;
        DATAOUT.MUAe.ns_PercentC(:,:,file) = ns_PercentC;
       

    
        % Calculate V1 ch boundaries
        granBtm = officLamAssign.Probe11stFold4c(file); % channel corresponding to the bottom of layer 4c
        v1Top_old = granBtm-9;
        v1Btm_old = granBtm+5;
        v1Ch_old = v1Top_old:v1Btm_old;

        numPenetrations = 30;
        laminarAligned = nan(1801,64,numPenetrations);
        alignDiff = 32-granBtm;
        v1Top_new = v1Top_old+alignDiff;
        v1Btm_new = v1Btm_old+alignDiff;
        v1Ch_new = v1Top_new:v1Btm_new;

        MUAe.ps_PercentC(v1Ch_new,:,file) = ps_PercentC(:,v1Ch_old);
        MUAe.ns_PercentC(v1Ch_new,:,file) = ns_PercentC(:,v1Ch_old);
        
    
         disp(strcat('finished processing _',sessionLabel))

end
load handel.mat
sound(y,1.15*Fs);
disp('Finished creating DATAOUT')
cd(outDir)
save('DATAOUT.mat',"DATAOUT")
datetime

%% plot grand averages of whole probe


    % convert 0 to NaN
    ps_NaNmatrix = DATAOUT.MUAe.ps_PercentC;
    ps_NaNmatrix(ps_NaNmatrix==0) = NaN;
    ns_NaNmatrix = DATAOUT.MUAe.ns_PercentC;
    ns_NaNmatrix(ns_NaNmatrix==0) = NaN;

    ps_grandAvg = median(ps_NaNmatrix,3,"omitmissing");
    ns_grandAvg = median(ns_NaNmatrix,3,"omitmissing");

    % Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
    ps_SEM = std(ps_NaNmatrix,0,3)/sqrt(size(ps_NaNmatrix,3)); 
    ns_SEM = std(ns_NaNmatrix,0,3)/sqrt(size(ns_NaNmatrix,3)); 

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
    titleText_grandAvg = {'_grandAvg'};
    s_grandAvg.Title = titleText_grandAvg;
    s_grandAvg.LineWidth = 1;
    % s_grandAvg.DisplayLabels = ["% Change"];
    
    %save fig
    cd(outDir)
    figName_grandAvg = strcat('stackedPlot_','_grandAvg_','.png');
    saveas(stk_grandAvg,figName_grandAvg)


%% plot grand averages of laminar compartments



    % convert 0 to NaN
    ps_NaNmatrix = DATAOUT.MUAe.ps_PercentC;
    ps_NaNmatrix(ps_NaNmatrix==0) = NaN;
    ns_NaNmatrix = DATAOUT.MUAe.ns_PercentC;
    ns_NaNmatrix(ns_NaNmatrix==0) = NaN;

    ps_grandAvg = median(ps_NaNmatrix,3,"omitmissing");
    ns_grandAvg = median(ns_NaNmatrix,3,"omitmissing");

    % Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
    ps_SEM = std(ps_NaNmatrix,0,3)/sqrt(size(ps_NaNmatrix,3)); 
    ns_SEM = std(ns_NaNmatrix,0,3)/sqrt(size(ns_NaNmatrix,3)); 

    %Now convert grand average array into laminar compartments
    ps_S = median(ps_grandAvg(:,1:5),2); 
    ps_G = median(ps_grandAvg(:,6:10),2); 
    ps_I = median(ps_grandAvg(:,11:15),2); 
    ns_S = median(ns_grandAvg(:,1:5),2); 
    ns_G = median(ns_grandAvg(:,6:10),2); 
    ns_I = median(ns_grandAvg(:,11:15),2); 


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
    tInput_ps_S = reshape(squeeze(median(ps_NaNmatrix(1400:1801,1:5,useIdx),1)),[],1);
    tInput_ps_G = reshape(squeeze(median(ps_NaNmatrix(1400:1801,6:10,useIdx),1)),[],1);
    tInput_ps_I = reshape(squeeze(median(ps_NaNmatrix(1400:1801,11:15,useIdx),1)),[],1);
    tInput_ns_S = reshape(squeeze(median(ns_NaNmatrix(1400:1801,1:5,useIdx),1)),[],1);
    tInput_ns_G = reshape(squeeze(median(ns_NaNmatrix(1400:1801,6:10,useIdx),1)),[],1);
    tInput_ns_I = reshape(squeeze(median(ns_NaNmatrix(1400:1801,11:15,useIdx),1)),[],1);

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

    titleText = {'Grand average of 135 multi-units per laminar compartment'};
    title(t,titleText,'Interpreter','none')

    %save fig
    cd(outDir)
    figName_lamCom = strcat('laminarCompartment_','_grandAvg_','.png');
    saveas(lamCom,figName_lamCom)



%% Notes
% What is currently working: 
% Step 0 - load in penetration list to analyze:
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
% Currently checking data quality for 211219_B - looks like noise

%% Quick pull out of data for individual session
% Note this was turned into the OrganizeData m file.
% Step 1 = load IDX
% Step 2 = chose NS and PS for BRFS
% Step 3, trial average
ps = 9;
ns = 10;
j = ps;
for trlNum = 1:length(IDX(j).correctTrialIndex)
    MUAe_ps(1:32,1:1001,trlNum) =  IDX(j).MUAe{trlNum,2}(1:1001,1:32)';
end
ps_avg = mean(MUAe_ps,3);

k = ns;
for trlNum = 1:length(IDX(k).correctTrialIndex)
    MUAe_ns(1:32,1:1001,trlNum) =  IDX(k).MUAe{trlNum,2}(1:1001,1:32)';
end
ns_avg = mean(MUAe_ns,3);

%% Laminar align data
cd('C:\Users\Brock Carlson\Box\Manuscripts\Maier')
load('DATAOUT.mat')
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
aligned_ps = nan(1801,64,size(officLamAssign,1));
aligned_ns = nan(1801,64,size(officLamAssign,1));
for i = 1:size(officLamAssign,1)
    % Calculate V1 ch boundaries
    granBtm = officLamAssign.Probe11stFold4c(i); % channel corresponding to the bottom of layer 4c
    v1Top_old = granBtm-9;
    v1Btm_old = granBtm+5;
    v1Ch_old = v1Top_old:v1Btm_old;
    
    numPenetrations = 30;
    alignDiff = 32-granBtm;
    v1Top_new = v1Top_old+alignDiff;
    v1Btm_new = v1Btm_old+alignDiff;
    v1Ch_new = v1Top_new:v1Btm_new;
    
    aligned_ps(:,v1Ch_new,i) = DATAOUT_ps(:,v1Ch_old,i);
    aligned_ns(:,v1Ch_new,i) = DATAOUT_ns(:,v1Ch_old,i);
end



%% Step 4, table organize
% Average data
ps_avg = median(DATAOUT_ps,3,"omitmissing"); % Avrage acros penetrations (other steps were done with median
ns_avg = median(DATAOUT_ns,3,"omitmissing");

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data));  
% % ps_SEM = std(DATAOUT_ps,0,3,"omitmissing")/sqrt(size(DATAOUT_ps,3)); 
% % ns_SEM = std(DATAOUT_ns,0,3,"omitmissing")/sqrt(size(DATAOUT_ns,3)); 
ps_SEM = std(DATAOUT_ps,0,3,"omitmissing")/sqrt(135); 
ns_SEM = std(DATAOUT_ns,0,3,"omitmissing")/sqrt(135); 


% smooth data
ps_smooth_mean = smoothdata(ps_avg,1,"gaussian",20);
ns_smooth_mean = smoothdata(ns_avg,1,"gaussian",20);
ps_smooth_sem = smoothdata(ps_SEM,1,"gaussian",20);
ns_smooth_sem = smoothdata(ns_SEM,1,"gaussian",20);

% convert to table
ps_table_1_mean = array2table(ps_smooth_mean(:,1:16));
ns_table_1_mean = array2table(ns_smooth_mean(:,1:16));
ps_table_2_mean = array2table(ps_smooth_mean(:,17:32));
ns_table_2_mean = array2table(ns_smooth_mean(:,17:32));

ps_table_1_sem = array2table(ps_smooth_sem(:,1:16));
ns_table_1_sem = array2table(ns_smooth_sem(:,1:16));
ps_table_2_sem = array2table(ps_smooth_sem(:,17:32));
ns_table_2_sem = array2table(ns_smooth_sem(:,17:32));

%Now convert ps_avg (a double array) into a table (input to timetable must
%be a table)
% The goal is to have 32 variables, each as a column, representing a
% different depth. 
Time = milliseconds(-200):milliseconds(1):milliseconds(1600);
% electrodes 1-16
elName_1 = {'el1','el2','el3','el4','el5','el6','el7','el8','el9',...
    'el10','el11','el12','el13','el14','el15','el16'};
ps_TT_1_mean = table2timetable(ps_table_1_mean,'RowTimes',Time);
preferredStimFlash_1_mean = renamevars(ps_TT_1_mean,ps_TT_1_mean.Properties.VariableNames,elName_1);
ns_TT_1_mean = table2timetable(ns_table_1_mean,'RowTimes',Time);
nullStimFlash_1_mean = renamevars(ns_TT_1_mean,ns_TT_1_mean.Properties.VariableNames,elName_1);
% electrodes 17-32
elName_2 = {'el17','el18','el19','el20','el21','el22','el23','el24',...
    'el25','el26','el27','el28','el29','el30','el31','el32'};
ps_TT_2_mean = table2timetable(ps_table_2_mean,'RowTimes',Time);
preferredStimFlash_2_mean = renamevars(ps_TT_2_mean,ps_TT_2_mean.Properties.VariableNames,elName_2);
ns_TT_2_mean = table2timetable(ns_table_2_mean,'RowTimes',Time);
nullStimFlash_2_mean = renamevars(ns_TT_2_mean,ns_TT_2_mean.Properties.VariableNames,elName_2);    

% SEM
ps_TT_1_sem = table2timetable(ps_table_1_sem,'RowTimes',Time);
preferredStimFlash_1_sem = renamevars(ps_TT_1_sem,ps_TT_1_sem.Properties.VariableNames,elName_1);
ns_TT_1_sem = table2timetable(ns_table_1_sem,'RowTimes',Time);
nullStimFlash_1_sem = renamevars(ns_TT_1_sem,ns_TT_1_sem.Properties.VariableNames,elName_1);
ps_TT_2_sem = table2timetable(ps_table_2_sem,'RowTimes',Time);
preferredStimFlash_2_sem = renamevars(ps_TT_2_sem,ps_TT_2_sem.Properties.VariableNames,elName_2);
ns_TT_2_sem = table2timetable(ns_table_2_sem,'RowTimes',Time);
nullStimFlash_2_sem = renamevars(ns_TT_2_sem,ns_TT_2_sem.Properties.VariableNames,elName_2);    

% Mean +/- SEM
% preferredStimFlash_1_sem
% nullStimFlash_1_sem
ps_mps_1 = preferredStimFlash_1_mean + preferredStimFlash_1_sem; %pref stim -- mean plus sem 1 
ps_mms_1 = preferredStimFlash_1_mean - preferredStimFlash_1_sem; %pref stim -- mean minus sem 1 
ns_mps_1 = nullStimFlash_1_mean + nullStimFlash_1_sem; %pref stim -- mean plus sem 1 
ns_mms_1 = nullStimFlash_1_mean - nullStimFlash_1_sem; %pref stim -- mean minus sem 1 

ps_mps_2 = preferredStimFlash_2_mean + preferredStimFlash_2_sem; %pref stim -- mean plus sem 1 
ps_mms_2 = preferredStimFlash_2_mean - preferredStimFlash_2_sem; %pref stim -- mean minus sem 1 
ns_mps_2 = nullStimFlash_2_mean + nullStimFlash_2_sem; %pref stim -- mean plus sem 1 
ns_mms_2 = nullStimFlash_2_mean - nullStimFlash_2_sem; %pref stim -- mean minus sem 1 


% Step 5, stackedplot()
close all
f = figure;
tl = tiledlayout(3,2);

nexttile(tl,[3 1])
s_1 = stackedplot(preferredStimFlash_1_mean,nullStimFlash_1_mean,...
    ps_mps_1,ps_mms_1,ns_mps_1,ns_mms_1);
s_1.LineWidth = 1;

nexttile(tl,[3 1])
s_2 = stackedplot(preferredStimFlash_2_mean,nullStimFlash_2_mean,...
    ps_mps_2,ps_mms_2,ns_mps_2,ns_mms_2);
s_2.LineWidth = 1;
