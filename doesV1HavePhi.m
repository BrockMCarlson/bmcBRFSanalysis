%doesV1HavePhi
% The goal of this script is to see if V1 has a reportable Phi value in
% MUAe.
% This script will take in MUAe data, generate a TPM, and then export that
% TPM for PhPhi application

cd('C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs')

load('filteredData_211008_B_bmcBRFS001.mat')

% Channels are chosen from laminarBoundaryCalculations.xlsx from SfN work
upperCh = 8:19;
middleCh = 20:25;
lowerCh = 26:32;

chA = 10;
chB = 20;
chC = 30;

smallChSelection(1,:) = MUAe(:,chA)';
smallChSelection(2,:) = MUAe(:,chB)';
smallChSelection(3,:) = MUAe(:,chC)';

medianCh = median(smallChSelection,2);

binarizedData = smallChSelection > medianCh;

possibleStates = combvec([0 1],[0 1],[0 1]);

for i = 1:8
    state = possibleStates(:,i);
    for j = 1:size(binarizedData,2)
        indexOfStates(i,j) = isequal(binarizedData(:,j),state);
    end
end

%% Now we can generate the TPM
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

%% Save TPM output
saveName = 'TPM_211008_B_bmcBRFS001';
save(saveName,"TPM")