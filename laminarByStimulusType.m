%% Laminar - stimulus type comparison


clear

%% Setup
sessionLabel = '211008_B_bmcBRFS001';

% Directories
codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\Brock Carlson\Documents\MATLAB\FormattedDataOutputs';

cd(outputDir)
load(strcat('sortedData_',sessionLabel,'.mat'))

% Variables
cond = 1; % Using most excitatory stimulus, condition 1, 'Simult. Dioptic. PO'
probeLength = size(IDX(cond).LFP_bb{1,1},2);
sdftm = STIM(1).sdftm;
timeLength = length(sdftm);
trlLength = size(IDX(cond).LFP_bb,1);
xAxisTime = sdftm;
xAxisLim = [-100 400];
yAxisChannels = 1:32;
yAxisLim =      [1 32];

% Figure settings
set(0,'DefaultFigureWindowStyle','docked')
f = figure;
sgtitle(sessionLabel,'interpreter','none')
% % f.Position = [1 1 2502 1224];

%% loop stimulus types
stimType = {'LFP_bb','LFP_alphaBeta','LFP_gamma','CSD_gamma','MUAe'};

for i = 1:5
    ax = subplot(1,5,i);
    counter = 0;
    for cond = 1:20 % There are 20 conditions, 20 rows in IDX.
        trlLength = size(IDX(cond).(stimType{i}),1);
        for trl = 1:trlLength
            counter = counter + 1;
            allData(:,:,counter) = IDX(cond).(stimType{i}){trl,1};
        end
    end
    
    respMedian = median(allData,3)';
    baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
    baselinemedian = median(respMedian(:,baselineTimeIndex),2);
    blSubAvg = respMedian - baselinemedian;
    imagesc(xAxisTime,yAxisChannels,blSubAvg)
    colormap(ax,'turbo');
    e = colorbar;
    if i <= 3
        e.Label.String = "uV";
        e.Label.Rotation = 270;
    elseif i == 4
        e.Label.String = '';
    end
    vline(0)
    xlim(xAxisLim)
    xlabel('Time (ms)');
    
    box off
    set(ax,'ytick',yAxisChannels); 
    
    title(stimType{i},'Interpreter','none')
end