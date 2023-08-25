%doesV1HavePhi
% The goal of this script is to see if V1 has a reportable Phi value in
% MUAe.
% This script will take in MUAe data, generate a TPM, and then export that
% TPM for PhPhi application

cd('S:\formattedDataOutputs')

load('filteredData_211008_B_bmcBRFS001.mat')
clearvars -except MUAe
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

% Now we can generate the TPM
TPM = nan(8,8);
for i = 1:8
    % given the state at time t, what is the probability that binarizedData
    % will be in each of the possible states at time t+1.
    