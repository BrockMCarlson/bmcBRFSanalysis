%% laminarLabeling.m
% Brock M. Carlson
% July 26th, 2023,
% This is a draft analysis for bmcBRFSanalysist to find the laminar 
% profile of a single bmcBRFS.ns2 file 

% LFP, CSD, PSD, Gamma/Beta, MUAe. All by depth

clear

%% Setup
sessionLabel = '211008_B_bmcBRFS001';

% Directories
codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
% outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
outputDir = 'S:\formattedDataOutputs';

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
f = figure;
sgtitle(sessionLabel,'interpreter','none')
% % f.Position = [1 1 2502 1224];

%% Average data from all trials

counter = 0;
for cond = 1:20 % There are 20 conditions, 20 rows in IDX.
    trlLength = size(IDX(cond).LFP_bb,1);
    for trl = 1:trlLength
        counter = counter + 1;
        lfpAllData(:,:,counter) = IDX(cond).LFP_bb{trl,1};
    end
end

lfpRespmedian = median(lfpAllData,3)';
baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
lfp_baselinemedian = median(lfpRespmedian(:,baselineTimeIndex),2);
lfp_blSubAvg = lfpRespmedian - lfp_baselinemedian;


%% LFP
ax1 = subplot(1,5,1);

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
xlabel('Time from stimulus onset (ms)');
set(gca,'YTick', [])
box off
title('LFP')


%% CSD
ax2 = subplot(1,5,2);

counter = 0;
for cond = 1:20
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
xlabel('Time (ms)');
ylabel('contact number'); 
set(ax2,'CLim',[-climit climit],'ytick',yAxisChannels); 
colormap(ax2,flipud(jet));

vline(0)
xlim(xAxisLim)
box off

c = colorbar;
c.Label.String = "(uV / mm)^2";
c.Label.Rotation = 270;


title('CSD')


%% PSD
ax3 = subplot(1,5,3);

% FFT
Fs       = 1000; % Hz
chanN    = size(lfpAllData,2); 
% loop through channels 
for ch = 1:chanN
    %loop through trials
    for trialNum = 1:size(lfpAllData,3)
        clear x n Spec
        lfp_holder        = lfpAllData(:,ch,trialNum);
    
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
                % fft --> absolute value --> squared.
                Xx    = abs(fft(xw,nfft)).^2; % Xx is power.
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
            power = nan(chanN,size(Spec,1),size(lfpAllData,3));
        end
        power(ch,:,trialNum) = Spec; 
    end
end
% permute power, initial design was (spec x chanN)
powerPermute = permute(power,[2 1 3]); %output is (freq x ch x trial)

% average power
powerAvg = mean(powerPermute,3,"omitnan");



%remove top channel artifact?
    warning('removing top channel artifact')
    powerAvg(:,1) = 0;


% normalize power. At each frequency, we divide the power of each channel
% by that of the channel with the highest power.
 power_norm = nan(size(powerAvg)); 
for freq = 1:size(powerAvg,1)
    powerOfAllChAtThisFreq = powerAvg(freq,:);
    [maxPower,chWithMaxPower] = max(powerOfAllChAtThisFreq);
    power_norm(freq,:) = powerOfAllChAtThisFreq ./ maxPower;
end

set(gcf,'color','w'); 
imagesc(freq_vector,yAxisChannels,power_norm'); 
colormap(ax3,'parula');
d = colorbar;
d.Label.String = "relative power";
d.Label.Rotation = 270;
d.Label.Position(1) = 4;


xlim([0 100]); 
xlabel('freq (Hz)'); 
set(ax3,'ytick',yAxisChannels)
set(ax3,'XTick',[5,12,40,100]); 
xlim([0 100]); 

title('PSD')

%% Theta-Alpha Gamma Cross
ax4 = subplot(1,5,4);
title('Theta-Alpha Gamma cross')

% Get the Gamma x Beta cross
% Beta is 12 - 20Hz (for our purposes)
% Gamma is 30-59,61:100
beta_index = (freq_vector > 5) & (freq_vector < 12);
gamma_index = (freq_vector > 40) & (freq_vector < 100);
% % % gamma_index(idx60hz) = false;
        
for j = yAxisChannels
    avgBeta(j,1) = mean(power_norm(beta_index,j));
    avgGamma(j,1) = mean(power_norm(gamma_index,j));
end

plot(avgBeta)
hold on
plot(fliplr(avgGamma))
set(ax4, 'XDir','reverse')
view([90 -90])
xlim(yAxisLim)
ylabel('Relative Power'); 
set(ax4,'xtick',yAxisChannels,'Box','off'); 
legend('5-12 Hz','40-100 Hz','Location','best')

title('Theta-Alpha Gamma cross')

%% MUAe
ax5 = subplot(1,5,5);
counter = 0;
for cond = 1:20 % There are 20 conditions, 20 rows in IDX.
    trlLength = size(IDX(cond).MUAe,1);
    for trl = 1:trlLength
        counter = counter + 1;
        muaAllData(:,:,counter) = IDX(cond).MUAe{trl,1};
    end
end

muaRespmedian = median(muaAllData,3)';
baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
mua_baselinemedian = median(muaRespmedian(:,baselineTimeIndex),2);
mua_blSubAvg = muaRespmedian - mua_baselinemedian;
imagesc(xAxisTime,yAxisChannels,mua_blSubAvg)
colormap(ax5,'turbo');
e = colorbar;
e.Label.String = "uV";
e.Label.Rotation = 270;

vline(0)
xlim(xAxisLim)
xlabel('Time (ms)');

box off
set(ax5,'ytick',yAxisChannels); 

title('MUAe')

