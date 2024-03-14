%% laminarLabeling.m
% Brock M. Carlson
% July 26th, 2023,
% This is a draft analysis for bmcBRFSanalysist to find the laminar 
% profile of a single bmcBRFS.ns2 file 

% LFP, CSD, PSD, Gamma/Beta, MUAe. All by depth

clear

%% Setup

% Directories
codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
dataDir = 'S:\bmcBRFS_sortedData_Nov23';
plotDir = 'S:\formattedDataOutputs';

%% Create dir loop
cd(dataDir)
sortedDataFiles = dir('*sortedData*');
for i = 1:length(sortedDataFiles)
    if i == 5 || i == 14
        disp('file #5 and #14 do not have full data')
        continue
    end
    fileToLoad = sortedDataFiles(i).name;
    load(fileToLoad)
    % Variables
    cond = 1; % Using most excitatory stimulus, condition 1, 'Simult. Dioptic. PO'
    sdftm = -200:1:800;
    timeLength = length(sdftm);
    trlLength = size(IDX(cond).LFP_bb,1);
    xAxisTime = sdftm;
    xAxisLim = [-100 400];
    if size(IDX(cond).LFP_bb{1,1},2) == 32
        probeNumber = 1;
    elseif size(IDX(cond).LFP_bb{1,1},2) == 64
        probeNumber = 2;
    end
    
    for probes = 1:probeNumber
        %% probe settings
        if probes == 1
            yAxisChannels = 1:32;
            yAxisLim =      [1 32];
        elseif probes == 2
            error('solve 2 probes plots')
        end

        %% Figure settings
        f = figure;
        set(gcf,"Position",[-1903 -89 1874 938]);
        sgtitle(sortedDataFiles(i).name(12:end-4),'interpreter','none')
        % % f.Position = [1 1 2502 1224];
        
        %% Average data from all trials
        
        counter = 0;
        clear lfpAllData
        for cond = 1:4 % There are 4 binocular conditions, which we will use for lamainr analysis
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
        if probes == 1
           contactToPlot = 32:-1:1;
        elseif probes == 2
           contactToPlot = 64:-1:33;
        end

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
        clear csd_allLFPdata
        try
        for cond = 1:20
            trlLength = size(IDX(cond).CSD_bb,1);
            for trl = 1:trlLength
                counter = counter + 1;
                csd_allLFPdata(:,:,counter) = IDX(cond).CSD_bb{trl,1};
            end
        end
        catch ME
            disp('no CSD_bb on file')
            disp(fileToLoad)
            continue
        end
        
        csdRespmedian = median(csd_allLFPdata,3)';
        csd_baselinemedian = median(csdRespmedian(:,baselineTimeIndex),2);
        csd_blSubAvg = csdRespmedian - csd_baselinemedian;

        CSDf = filterCSD(csd_blSubAvg);

        
        imagesc(xAxisTime,yAxisChannels,CSDf)
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
        
        %% Coherence
        ax4 = subplot(1,5,4);
        title('Coherence')
        conditionNumber = 1;
        
        % Get the number of trials for the chosen condition
        numTrials = size(IDX(conditionNumber).LFP_bb, 1);
        
        % Parameters for mscohere
        fs = 1000;  % Sampling frequency in Hz
        windowSize = 256;  % Window size for computing the coherence
        overlap = windowSize/2;  % Overlap between windows
        tm_coher = 490:1001; % Time window of data. Last 512ms of trial. 
        tm = 1:1001; % Time window of data. Last 512ms of trial. 
        
        % Initialize coherence matrix
        coherenceMatrix = nan(32,32, numTrials);
        % Loop through all trials and compute coherence for each channel pair
        for trialNumber = 1:numTrials
            for channel1 = 1:32
                for channel2 = 1:channel1
                    % Extract the LFP_bb data for the chosen channels and trial
                    lfpGammaData1 = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, channel1);
                    lfpGammaData2 = IDX(conditionNumber).LFP_bb{trialNumber, 2}(tm, channel2);
                    %baseline subtract
                    bl1 = median(lfpGammaData1(1:200));
                    lfp_blsub_1 = lfpGammaData1 - bl1;
                    bl2 = median(lfpGammaData2(1:200));
                    lfp_blsub_2 = lfpGammaData2 - bl2;
    
                    % Compute coherence
                    [coherence, freq] = mscohere(lfp_blsub_1(tm_coher), lfp_blsub_2(tm_coher), windowSize, overlap, [], fs);
        
                    % Store coherence in the matrix
                    coherenceMatrix(channel1, channel2, trialNumber) = median(coherence(2:15));  % You can use median or any other aggregation method
                end
            end
        end

        
        % Average across trials and save output
        averageCoherenceMatrix = median(coherenceMatrix,3);
        
        % Visualize the coherence matrix
        imagesc(averageCoherenceMatrix);
        colormap(ax4,'jet');
        e = colorbar;
        e.Label.String = "Inter-contact coherence";
        e.Label.Rotation = 270;
        xlabel('Channel');
        ylabel('Channel');
        title('Coherence, 4-55Hz');

        
        %% MUAe
        ax5 = subplot(1,5,5);
        counter = 0;
        clear muaAllData
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
        figName = strcat('laminarFig',sortedDataFiles(i).name(12:end-4),'_Probe#',string(probes),'.png');
        saveas(gcf,figName)

        close
    end
end


