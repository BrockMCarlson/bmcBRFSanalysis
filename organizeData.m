%% Organize Data
% The goal of this script is to load in all SortedData files and store just
% the MUA outputs (or LFP_bb on a seperate run) into a single matlab
% variable. 

% The output variable may be around 50 GB, but should dramatically increase
% ease of data workability. 

% Data should be organized on a per-penetration basis. So sessions with two
% probes should get seperated with this script

% 211027_B_ brfs001 and brfs002 should be concatenated
% 221202, 221206 and 221209 should perhaps not be concatenated as I think
% blake may have moved the stimulus on the screen between data recordings
% (I'm not quite sure...)

% The output format should be (tm x ch x trls) = (1001 x 32 x N) for each
% trial-type. Each row of the struct variable should be a different
% penetration. We should have a field for the data file name and a filed
% for the probe# (each penetration should get its own row). Each trial-type
% should then also have its own field. The output struct will be Num
% penetration x 1 struct with 22 fields.

%% What is this!? The goal here needs to be to quickly get the mean and variance for each condition. 

%%
clear
tic
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)

stringsForFields = {'diop PO','diop NO','dichop POL','Dichop POR',...
    'monoc adapt diop','','','','','','','','','','','','','','',''};
count = 0;
sortedDataFiles = dir('*sortedData*');
for i = 1:length(sortedDataFiles)
    if i == 5 
        fileToLoad1 = sortedDataFiles(5).name;
        fileToLoad2 = sortedDataFiles(6).name;
        brfs001 = load(fileToLoad1);
        brfs002 = load(fileToLoad2);
        error('file #5 must be concatenated with file #6')
    end
    if i == 6 % handled in previous step
        continue
    end
    if i == 14
        disp('file #14 does not have full data')
        continue
    end    
    fileToLoad = sortedDataFiles(i).name;
    load(fileToLoad)
    % Probe number
    if size(IDX(1).LFP_bb{1,1},2) == 32
        probeNumber = 1;
    elseif size(IDX(1).LFP_bb{1,1},2) == 64
        probeNumber = 2;
    end
    
    % Ok, lets get down to buisness
    for probes = 1:probeNumber
        count = count + 1;
        allMUAdata(count).fileName = sortedDataFiles(i).name;
        allMUAdata(count).probe = probes;
        conditionString = {IDX.conditionString}.';
        for j = 1:20
            clear holderMUAe
            for trlNum = 1:length(IDX(j).correctTrialIndex)
                if probes == 1
                    holderMUAe(1:1001,:,trlNum) =  IDX(j).MUAe{trlNum,1}(:,1:32);
                    holderMUAe(1001:1801,:,trlNum) = IDX(j).MUAe{trlNum,2}(201:end,1:32);
                elseif probes == 2
                    holderMUAe(1:1001,:,trlNum) =  IDX(j).MUAe{trlNum,1}(:,33:64);
                    holderMUAe(1001:1801,:,trlNum) = IDX(j).MUAe{trlNum,2}(201:end,33:64);
                end
            end
            allMUAdata(count).(string(j)) = holderMUAe;

        end
    end


end

%% Save Data
outputDir = 'S:\formattedDataOutputs';
save('allMUAdata.mat',allMUAdata,"-v7.3")
toc

%%


% MUAe, LFP_bb, and CSD_bb. Three types. 
% PS and NS varibles. 6 total variables
% Format: (Ch,tm,trls,penetration)
% Challenege, choosing ps and NS, and aligning ch across penetrations.
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

