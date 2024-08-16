%% Make all session var (LFP_trials/MUA_trials) from each session IDX
clear
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir     = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'S:\TrialTriggeredLFPandMUA';
    plotDir     = '';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir     = '';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)

%% For loop
    cd(dataDir)
for file = 1:size(officLamAssign,1)
    % load data
    probeName = char(officLamAssign.Session_probe_(file,1));
    fileToLoad = strcat('trialTriggeredData_',probeName(1:19),'.mat');
    load(fileToLoad)
    LFP_trials{file,1} = {IDX.LFP_bb}.';
    MUA_trials{file,1} = {IDX.MUAe}.';

    disp(strcat('Done with file number: ',string(file)))
    memory
end
%% save LFP_trials output
save('LFP_trials.mat','LFP_trials','-v7.3')
save('MUA_trials.mat','MUA_trials','-v7.3')



