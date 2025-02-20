%% laminarLabeling.m
% Brock M. Carlson
% July 26th, 2023,
% This is a draft analysis for bmcBRFSanalysist to find the laminar 
% profile of a single bmcBRFS.ns2 file 

% LFP, CSD, PSD, Gamma/Beta, MUAe. All by depth

clear

%% Setup
disp('start time')
datetime
close
clearvars -except MUA_trials
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)
cd(dataDir)

%% Create dir loop
cd(dataDir)
allDataFiles = dir('**/*trialTriggeredData*.mat');

% % for file = 20:size(officLamAssign,1)
for file = 2
    
    % load data
    cd(dataDir)
    probeName = char(officLamAssign.Session_probe_(file,1));
    fileToLoad = strcat('trialTriggeredData_',probeName(1:19),'.mat');
    load(fileToLoad)
    chStr = string(officLamAssign.ChToUse(file,1));
    if strcmp(chStr,"1:32")
        yAxisChannels = 1:32;
        yAxisLim =      [1 32];
       contactToPlot = 32:-1:1;
    elseif strcmp(chStr,"33:64")
        yAxisChannels = 33:64;
        yAxisLim =      [33 64];
        contactToPlot = 64:-1:33;
    end
    
    
    
    % Variables
    sdftm = -200:1:1800;
    timeLength = length(sdftm);
    xAxisTime = sdftm;
    xAxisLim = [-100 400];
    
    
    
    
    %% Figure settings
    f = figure;
    set(gcf,"Position",[158 504 1237 306]);
    sgtitle(probeName,'interpreter','none')
    % % f.Position = [1 1 2502 1224];
    
        
    %% LFP
    
    counter = 0;
    clear lfpAllData
    % for cond = 1:4 % There are 4 binocular conditions, which we will use for lamainr analysis
    for cond = 1:2 
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
    
    
    
    ax1 = subplot(1,5,1);
    
    counter = 0;
    for contact = contactToPlot
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
    clear csd
    % calculate CSD from LFP
    for trlNum = 1:size(lfpAllData,3)
        csd(:,:,trlNum) = calcCSD_classic(lfpAllData(:,:,trlNum));
    end

    
    csdRespmedian = median(csd,3)';
    csd_baselinemedian = median(csdRespmedian(:,baselineTimeIndex),2);
    csd_blSubAvg = csdRespmedian - csd_baselinemedian;

    CSDf = filterCSD(csd_blSubAvg);

    
    imagesc(xAxisTime,yAxisChannels,CSDf)
    climit = max(abs(get(gca,'CLim'))*.5);
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
    % loop through channels 
    for chIDX = 1:32
        ch = yAxisChannels(chIDX);
        %loop through trials
        for trialNum = 1:size(lfpAllData,3)
            clear x n Spec
            lfp_holder  = lfpAllData(:,ch,trialNum);
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
            if chIDX == 1
                    power = nan(length(yAxisChannels),size(Spec,1),size(lfpAllData,3));
            end
            freq_vector = (select - 1)*Fs/nfft;
            power(chIDX,:,trialNum) = Spec; 
        end
    end
    % permute power, initial design was (spec x chanN)
    powerPermute = permute(power,[2 1 3]); %output is (freq x ch x trial)
    
    % average power
    powerAvg = mean(powerPermute,3,"omitnan");       
    
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
    
    
    xlim([0 150]); 
    xlabel('freq (Hz)'); 
    set(ax3,'ytick',yAxisChannels)
    % set(ax3,'XTick',[5,12,40,100]); 
    xlim([0 100]); 
    
    title('PSD')

    %% Gamma x Beta Gamma Cross
    % % ax4 = subplot(1,5,4);
    % % title('Alpha x Gamma cross')
    % % 
    % % % Get the Gamma x Beta cross
    % % % Beta is 12 - 20Hz (for our purposes)
    % % % Gamma is 30-59,61:100
    % % beta_index = (freq_vector >= 4) & (freq_vector <= 12);
    % % gamma_index = (freq_vector > 30) & (freq_vector < 100);
    % % % % % gamma_index(idx60hz) = false;
    % % 
    % % for j = 1:32
    % %     avgBeta(j,1) = mean(power_norm(beta_index,j));
    % %     avgGamma(j,1) = mean(power_norm(gamma_index,j));
    % % end
    % % 
    % % plot(avgBeta)
    % % hold on
    % % % plot(fliplr(avgGamma))
    % % plot(avgGamma)
    % % set(ax4, 'XDir','reverse')
    % % view([90 -90])
    % % xlim([1 32])
    % % xticks(1:32)
    % % xticklabels(yAxisChannels)
    % % ylabel('Relative Power'); 
    % % % set(ax4,'xtick',yAxisChannels,'Box','off'); 
    % % legend('4-12 Hz','30-100 Hz','Location','northeast')
    
    
  %% Coherence
    ax4 = subplot(1,5,4);
    title('Dioptic coherence')

    % Parameters for mscohere
    fs = 1000;  % Sampling frequency in Hz
    windowSize = 256;  % Window size for computing the coherence
    overlap = windowSize/2;  % Overlap between windows
    % tm_coher = 1:512; % Time window of data. Last 512ms of trial. 
    tm_coher = 490:1001; % Time window of data. Last 512ms of trial. 
    tm = 1:1001; % Time window of data. Last 512ms of trial. 

    % Get the number of trials for the chosen condition
    totalTrials = size(IDX(1).LFP_bb, 1) + size(IDX(2).LFP_bb, 1);
    
    
    % Initialize coherence matrix
    coherenceMatrix1 = nan(32,32, totalTrials);
    % Loop through all trials and compute coherence for each channel pair
    count = 0;
    for conditionNumber = 1:2
        for trialNumber = 1:size(IDX(conditionNumber).LFP_bb, 1)
            count = count + 1;
            for channel1 = yAxisChannels
                for channel2 = yAxisChannels(1):channel1
                    % Extract the LFP_bb data for the chosen channels and trial
                    lfpGammaData1 = IDX(conditionNumber).LFP_bb{trialNumber, 1}(tm, channel1);
                    lfpGammaData2 = IDX(conditionNumber).LFP_bb{trialNumber, 1}(tm, channel2);
                    %baseline subtract
                    bl1 = median(lfpGammaData1(1:200));
                    lfp_blsub_1 = lfpGammaData1 - bl1;
                    bl2 = median(lfpGammaData2(1:200));
                    lfp_blsub_2 = lfpGammaData2 - bl2;
        
                    % Compute coherence
                    [coherence, freq] = mscohere(lfp_blsub_1(tm_coher), lfp_blsub_2(tm_coher), windowSize, overlap, [], fs);
        
                    % Store coherence in the matrix
                    if strcmp(chStr,"1:32")
                        coherenceMatrix1(channel1, channel2, count) = median(coherence(2:15));  % You can use median or any other aggregation method
                    elseif strcmp(chStr,"33:64")
                        coherenceMatrix1(channel1-32, channel2-32, count) = median(coherence(2:15));  % You can use median or any other aggregation method
                    end
                end
            end
        end
    end
    
    
    % Average across trials and save output
    averageCoherenceMatrix = median(coherenceMatrix1,3);
    
    % Visualize the coherence matrix
    imagesc(averageCoherenceMatrix);
    colormap(ax4,'jet');
    set(ax4,'ytick',1:32)
    set(ax4,'yticklabels',yAxisChannels)
    e = colorbar;
    e.Label.String = "Inter-contact coherence";
    e.Label.Rotation = 270;
    xlabel('Channel');
    ylabel('Channel');
    title('Coherence, Dioptic');
    
    
    %% MUAe
    ax5 = subplot(1,5,5);
    counter = 0;
    clear muaAllData
    % for cond = 1:20 % There are 20 conditions, 20 rows in IDX.
    for cond = 1:2 % There are 20 conditions, 20 rows in IDX.
        trlLength = size(IDX(cond).MUAe,1);
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllData(:,:,counter) = IDX(cond).MUAe{trl,1}(:,yAxisChannels);
        end
    end
    
    muaRespmedian = median(muaAllData,3)';
    baselineTimeIndex = find(sdftm == -50):find(sdftm == 0); 
    mua_baselinemedian = median(muaRespmedian(:,baselineTimeIndex),2);
    mua_blSubAvg = muaRespmedian - mua_baselinemedian;
    imagesc(xAxisTime,yAxisChannels,mua_blSubAvg)
    colormap(ax5,'turbo');
    f = colorbar;
    f.Label.String = "dB - band-limited uV";
    f.Label.Rotation = 270;
    
    vline(0)
    xlim(xAxisLim)
    xlabel('Time (ms)');
    
    box off
    set(ax5,'ytick',yAxisChannels); 
    
    title('MUAe')

     %% Save figure
    cd(plotDir)
    figName_p = strcat('2025laminarFig_',probeName(1:end-1),'.png');
    saveas(gcf,figName_p)
    figName_s = strcat('2025laminarFig_',probeName(1:end-1),'.svg');
    saveas(gcf,figName_s)



end


