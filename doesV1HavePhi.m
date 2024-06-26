%doesV1HavePhi
% The goal of this script is to see if V1 has a reportable Phi value in
% MUAe.
% This script will take in MUAe data, generate a TPM, and then export that
% TPM for PhPhi application

clear
% Directories
path1 = strcat(getenv('HOMEDRIVE'),getenv("HOMEPATH"),filesep);
path2 = 'Documents\GitHub\bmcBRFSanalysis';
codeDir = strcat(path1,path2);
cd(codeDir)
path3 = 'Documents\MATLAB\formattedDataOutputs';
plotDir = strcat(path1,filesep,path3);
%sortedDataDir -- need to run new sorted data on home pc
dataDir = 'D:\';

officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);


%% Choose inidividual file to explore
cd(dataDir)
load('DATAOUT_trials.mat')
penetrationNum = 8;
granBtm = table2array(officLamAssign(8,"Probe11stFold4c"));


% Channels are chosen from laminarBoundaryCalculations.xlsx from SfN work
upperCh = granBtm-11:granBtm-6;
middleCh = granBtm-5:granBtm;
lowerCh = granBtm+1:granBtm+6;

chA = granBtm-8;
chB = granBtm;
chC = granBtm + 5;

dataType = {'CSD_bb','CSD_beta1','CSD_gamma1','LFP_bb','LFP_beta1','LFP_gamma1','MUAe'};
for dt = 1:7
    figure
    if dt == 1 || dt == 2 || dt == 3
        smallChSelection(1,:) = IDX.(dataType{dt})(:,chA)';
        smallChSelection(2,:) = input.(dataType{dt})(:,chB)';
        smallChSelection(3,:) = input.(dataType{dt})(:,chC)';
    elseif dt >= 4 
        smallChSelection(1,:) = input.(dataType{dt})(:,chA)';
        smallChSelection(2,:) = input.(dataType{dt})(:,chB)';
        smallChSelection(3,:) = input.(dataType{dt})(:,chC)';
    end
    medianCh = median(smallChSelection,2);
    
    tm = 1:500;

    tiledlayout("vertical",'TileSpacing','Compact')
    for i = 1:3
        %subplot(3,1,i)
        nexttile
        plot(tm,smallChSelection(i,3001:3500))
        hold on
        yline(medianCh(i))
        set(gca,"Box","off")
        ylabel('laminar channel, uV')
        if i < 3
            set(gca,'xtick',[])
        end
    end
    xlabel('Time, (ms)')
    sgtitle(dataType{dt},'Interpreter','None')
    
    binarizedData = smallChSelection > medianCh;
    nexttile
    spy(binarizedData(:,3001:3500))
    title('Binarized Data')
    ylim([0 4])
    yticklabels({'upper','middle','deep'})
    box off
    xlabel('Number of samples: 1 sample / ms')
    
    %How often is the data in the "on-state"?
    numberofCellsInData = size(binarizedData,1)*size(binarizedData,2);
    timesDataIsOn = sum(binarizedData,'all');
    percentOfSamplesInOnState = timesDataIsOn / numberofCellsInData;
    
    possibleStates = combvec([0 1],[0 1],[0 1]);
    
    indexOfStates = nan(size(possibleStates,2),length(smallChSelection));
    for i = 1:8
        state = possibleStates(:,i);
        for j = 1:size(binarizedData,2)
            indexOfStates(i,j) = isequal(binarizedData(:,j),state);
        end
    end

    
    nexttile
    spy(indexOfStates(:,3001:3500))
    set(gca,"Box","off")
    title('Index of state occupied')
    yticks([1 4 8])
    yticklabels({'1','4','8'})
    ylabel('State of substrate')
    
    
    % % nexttile
    % % spy(indexOfStates(:,3101:3125))
    % % set(gca,"Box","off")
    % % title({'Zoomed-in index of state occupied',dataType{dt}})
    % % 
    % % 
    % % 
    %% Now we can generate a state-by-state TPM
        % given the state at time t, what is the probability that binarizedData
        % will be in each of the possible states at time t+1.
    TPM = nan(8,8);
    for i = 1:8
        % for each state, we find the timepoints that that state occurs
        timeThatState_i_Occurs = find(indexOfStates(i,:));
        % We get the number of times this state occurs in the data
        numberOfStatePresentations = length(timeThatState_i_Occurs);
        % Then we get the state at time t+1
        timePointPlus1 = timeThatState_i_Occurs + 1;
        try
            stateAtTimePointPlus1 = indexOfStates(:,timePointPlus1);
        catch
            stateAtTimePointPlus1 = indexOfStates(:,timePointPlus1(1:end-1));
        end
        % Now we calculate the porportion that each state happens at t+1
        numberOfStatePresentationsAtTimePlus1 = sum(stateAtTimePointPlus1,2);
        proportionEachStateAppearsAtTimePlus1 = numberOfStatePresentationsAtTimePlus1./numberOfStatePresentations;
        TPM(i,:) = proportionEachStateAppearsAtTimePlus1;
    end
    
    nexttile([2 1])
    heatmap(TPM)
    title('TPM')
    colormap('parula');
    xlabel('State of substrate')
    ylabel('State of substrate')

    %% Save TPM output
%     saveName = strcat('TPM_bmcBRFS001_',dataType{dt});
%     save(saveName,"TPM")

end

