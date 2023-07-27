%% laminarLabeling.m
% Brock M. Carlson
% July 26th, 2023,
% This is a draft analysis for bmcBRFSanalysist to find the laminar 
% profile of a single bmcBRFS.ns2 file 

% LFP, CSD, PSD, Gamma/Beta, MUAe. All by depth

%% Setup
% Directories
clear
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';

cd(outputDir)
load('sortedData_211008_B_bmcBRFS001.mat')

% Variables
cond = 1; % Using most excitatory stimulus, condition 1, 'Simult. Dioptic. PO'
probeLength = size(IDX(cond).LFP_bb{1,1},2);
sdftm = STIM(1).sdftm;
timeLength = length(sdftm);
trlLength = size(IDX(cond).LFP_bb,1);
baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
xAxisTime = sdftm;
xAxisLim = [-100 400];
yAxisChannels = 1:32;

% Figure settings
set(0,'DefaultFigureWindowStyle','docked')
f = figure;
% % f.Position = [1 1 2502 1224];

%% LFP
subplot(1,5,1)
lfpAllTrials = nan(timeLength,probeLength,trlLength);
for trl = 1:trlLength
    lfpAllTrials(:,:,trl) = IDX(cond).LFP_bb{trl,1};
end
lfpRespmedian = median(lfpAllTrials,3)';
lfp_baselinemedian = median(lfpRespmedian(:,baselineTimeIndex),2);
lfp_blSubAvg = lfpRespmedian - lfp_baselinemedian;

counter = 0;
for contact = probeLength:-1:1
    counter = counter + 1;
    chResponse = lfp_blSubAvg(contact,:);
    yValueAdder = chResponse + 150 * counter;
    plot(sdftm,yValueAdder,'k')
    hold on
end

vline(0)

ylim([25 4900])
xlim(xAxisLim)
set(gca,'YTick', [])
box off

%% CSD
subplot(1,5,2)
% % csdAllTrials = nan(timeLength,probeLength,trlLength);
% % for trl = 1:trlLength
% %     csdAllTrials(:,:,trl) = IDX(cond).CSD_gamma{trl,1};
% % end

conditionLength = 20; % There are 20 conditions, 20 rows in IDX.
counter = 0;
for cond = 1:conditionLength
    trlLength = size(IDX(cond).CSD_gamma,1);
    for trl = 1:trlLength
        counter = counter + 1;
        csd_allLFPdata(:,:,counter) = IDX(cond).CSD_gamma{trl,1};
    end
end

csdRespmedian = median(csd_allLFPdata,3)';
csd_baselinemedian = median(csdRespmedian(:,baselineTimeIndex),2);
csd_blSubAvg = csdRespmedian - csd_baselinemedian;

imagesc(xAxisTime,yAxisChannels,csd_blSubAvg)
climit = max(abs(get(gca,'CLim'))*.8);
ydir = 'reverse';
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
colormap(flipud(jet));

vline(0)
xlim(xAxisLim)
box off
c = colorbar;


%% PSD

% % You are here - calculate PSD from IDX.LFP_bb triggered data. 

subplot(1,5,3)
conditionLength = 20; % There are 20 conditions, 20 rows in IDX.
counter = 0;
for cond = 1:conditionLength
    trlLength = size(IDX(cond).LFP_bb,1);
    for trl = 1:trlLength
        counter = counter + 1;
        psd_allLFPdata(:,:,counter) = IDX(cond).LFP_bb{trl,1};
    end
end

% FFT
Fs       = 1000; % Hz
chanN    = size(psd_allLFPdata,2); 
% loop through channels 
for ch = 1:chanN
    %loop through trials
    for trialNum = 1:size(psd_allLFPdata,3)
        clear x n Spec
        lfp_holder        = psd_allLFPdata(:,ch,trialNum);
    
        %%%%%%%%%%%%%%%%%
        % % % % % % % % % %  x =bandStopFiltLFP(lfp_holder); fix bandstop filter
        %%%%%%%%%%%%%%%%%
    
        n        = size(lfp_holder,1); % Number of data points
        % prep for psd 
            nfft     = 512; 
            window   = hanning(nfft); 
            nwind    = length(window); 
        if n < nwind    % zero-pad x if it has length less than the window length
            x(nwind)=0;  
            n=nwind;
        end
        noverlap = 1;
        k        = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
        index    = 1:nwind;
        % compute PSD
            x        = lfp_holder;   
            Spec     = zeros(nfft,1); 
        % loop through windows 
            for j=1:k
                xw    = window.*(x(index));
                index = index + (nwind - noverlap);
                Xx    = abs(fft(xw,nfft)).^2;
                Spec  = Spec + Xx;  
            end
        % Select first half
            if ~any(any(imag(x)~=0))   % check if x is complex 
                if rem(nfft,2)    % nfft odd
                    select = (1:(nfft+1)/2)';
                else
                    select = (1:nfft/2+1)';
                end
                Spec = Spec(select);
            else
                select = (1:nfft)';
            end
        freq_vector = (select - 1)*Fs/nfft;
        if ch == 1
            power = nan(chanN,size(Spec,1),size(psd_allLFPdata,3));
        end
        power(ch,:,trialNum) = Spec; 
    end
end
% permute power, initial design was (spec x chanN)
powerPermute = permute(power,[2 1 3]);

% average power
powerAvg = mean(powerPermute,3,"omitnan");


% cheat at your band stop filter - remove 60Hz artifact manually
idx60hz = find((freq_vector > 57.4 & freq_vector < 62.6 ));
powerAvg(idx60hz,:) = 0;
%remove top channel artifact?
    % warning('removing top channel artifact - necessary? BMC DEV!!')
    % powerAvg(:,1) = 0;

% normalize power @ each frequency relative to power across contacts 
 power_norm = nan(size(powerAvg)); 
 for ch = 1:size(powerAvg,2)
     for f = 1:size(powerAvg,1)
         power_norm(f,ch) = (powerAvg(f,ch) - mean(powerAvg(f,:)))./(mean(powerAvg(f,:))) * 100; % percent deviation from mean
     end
 end

set(gcf,'color','w'); 
imagesc(freq_vector,yAxisChannels,power_norm'); 
colormap('hot'); xlim([0 100]); 
xlabel('freq (Hz)'); ylabel('contact number'); 
set(gca,'tickdir','out','ytick',yAxisChannels); 

%% Gamma Beta Cross
subplot(1,5,4)

% Get the Gamma x Beta cross
% Beta is 12 - 20Hz (for our purposes)
% Gamma is 30-59,61:100
beta_index = (freq_vector > 12) & (freq_vector < 30);
gamma_index = (freq_vector > 30) ;
gamma_index(idx60hz) = false;
        
for j = yAxisChannels
    avgBeta(j,1) = mean(power_norm(beta_index,j));
    avgGamma(j,1) = mean(power_norm(gamma_index,j),"omitnan");
end

plot(avgBeta)
hold on
plot(fliplr(avgGamma))
view([90 -90])
set(gca,'xdir','reverse')
legend('Beta','Gamma','Location','best')


%% MUAe
subplot(1,5,5)

muaAllTrials = nan(timeLength,probeLength,trlLength);
for trl = 1:trlLength
    muaAllTrials(:,:,trl) = IDX(cond).MUAe{trl,1};
end
muaRespmedian = median(muaAllTrials,3)';
mua_baselinemedian = median(muaRespmedian(:,baselineTimeIndex),2);
mua_blSubAvg = muaRespmedian - mua_baselinemedian;
imagesc(xAxisTime,yAxisChannels,mua_blSubAvg)
vline(0)
xlim(xAxisLim)
box off