%% Make all session var (LFP_trials/MUA_trials) from each session IDX
clear
close all
codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures\grand average');
cd(codeDir)

%% For loop
dataDir = 'D:\sortedData_240229';
% dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\sortedData_240229';
% % dataDir = 'S:\bmcBRFS_sortedData_Nov23';
cd(dataDir)
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);



for file = 1:size(officLamAssign,1)
    % load data
    cd(dataDir)
    probeName = char(officLamAssign.Session_probe_(file,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');
    load(fileToLoad)
    LFP_trials{file} = {IDX.LFP_bb}.';
    MUA_trials{file} = {IDX.MUAe}.';

    disp(strcat('Done with file number: ',string(file)))
end
%% save LFP_trials output
cd(dataDir)
save('LFP_trials.mat','LFP_trials','-v7.3')
save('MUA_trials.mat','MUA_trials','-v7.3')
save('prefOutput.mat','prefOutput','-v7.3')



