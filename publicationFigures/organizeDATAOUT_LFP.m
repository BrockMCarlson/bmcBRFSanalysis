%% Organize data

datetime

%% Setup
clear
% Directories

% codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures');
cd(codeDir)
% outDir = 'S:\formattedDataOutputs';
outDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
% dataDir = 'S:\bmcBRFS_sortedData_Nov23';
dataDir = 'D:\sortedData_240229';
cd(dataDir)
% officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);

%% For loop

DATAOUT_ps = nan(1801,32,size(officLamAssign,1));
DATAOUT_ns = nan(1801,32,size(officLamAssign,1));
DATAOUT_numTuned = nan(size(officLamAssign,1),1);
for file = 1:size(officLamAssign,1)
    
    %File specific info
    probeName = char(officLamAssign.SessionProbe(file,1));
    chStr = officLamAssign.ChtoUse(file,1);
    if strcmp(chStr,"1:32")
        ch = 1:32;
    elseif strcmp(chStr,"33:64")
        ch = 33:64;
    end

    % load data
    cd(dataDir)
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
        DATAOUT_numTuned(file,1) = sum(tuned);
        if sum(tuned) == 0
            warning(strcat('No tuned channels on_',probeName))
            continue
        end

        % Creat channel array to use in subsequent steps - only plot tuned units.
        chTuned = ch(logical(tuned));


    
    
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
                array_dichopticAdapted_pref(tm1,i,trl) = IDX(prefCondOnFlash).LFP_bb{trl,1}(tm1,i); % now we index the cell for the second 800ms
                array_dichopticAdapted_pref(tm2_concat,i,trl) = IDX(prefCondOnFlash).LFP_bb{trl,2}(tm2,i); % now we index the cell for the second 800ms
            end
            for trl = 1:length(IDX(nullCondOnFlash).correctTrialIndex)
                array_dichopticAdapted_null(tm1,i,trl) = IDX(nullCondOnFlash).LFP_bb{trl,1}(tm1,i); % now we index the cell for the second 800ms
                array_dichopticAdapted_null(tm2_concat,i,trl) = IDX(nullCondOnFlash).LFP_bb{trl,2}(tm2,i); % now we index the cell for the second 800ms
            end
        end
        % output is tm, ch, trl
        
        
        %% Save data into structure array
      
        % mean across trials
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
        DATAOUT_ps(:,:,file) = ps_PercentC;
        DATAOUT_ns(:,:,file) = ns_PercentC;
       
       % % 
       % % 
       % %  % Calculate V1 ch boundaries
       % % granBtm = officLamAssign.Probe11stFold4c(file); % channel corresponding to the bottom of layer 4c
       % %  v1Top = granBtm - 14;
       % %  v1Btm = granBtm + 10;
       % %  v1Ch = v1Top:v1Btm;
       % % 
       % %  LFP_bb.ps_PercentC(:,:,file) = ps_PercentC(:,v1Ch);
       % %  LFP_bb.ns_PercentC(:,:,file) = ns_PercentC(:,v1Ch);
       % % 
    
         disp(strcat('finished processing _',probeName))

end
load handel.mat
sound(y,1.15*Fs);
disp('Finished creating DATAOUT')
cd(outDir)
save('DATAOUT_LFP_trials.mat','DATAOUT_ps','DATAOUT_ns','DATAOUT_numTuned', '-v7.3')
datetime