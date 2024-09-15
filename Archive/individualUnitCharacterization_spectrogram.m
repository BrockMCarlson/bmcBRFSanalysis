%% individualUnitCharacterization.m
% Brock M. Carlson
% July 20th, 2023,
% This is a draft analysis for bmcBRFSanalysis. The questions we are
% interested in at an individual unit level are:
% 1. Visually responseive?
% 2. Feature Selective (Tuned to eye? Ori? Eye And Ori?)
% 3. Supressed by dichoptic stimuli?
% 4. Perceptually modulated? 
% These questions are all based on statistic run across presentation
% (trial) from a binned response window, for transient and sustained.

% Additionally, we want to find the laminar position of each unit.

% For continuous data, we want the following plots
    % 1. 2x2 eye v ori (potentially on the same ordinate) 
    % 2. dichoptic vs dichoptic
    % 3. Same stimulus, different history, Preferred vs null
    % 4. Same stimulus, different history, mixed vs mixed.

%% Setup
clear
close all
codeDir = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis'; 
cd(codeDir)
outputDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';

cd(outputDir)
fileName = 'sortedData_211008_B_bmcBRFS001.mat';
load(fileName)

% Variables
probeLength = size(IDX(1).LFP_bb{1,1},2);
sdftm = STIM(1).sdftm;
timeLength = length(sdftm);
baselineTimeIndex = find(sdftm == -100):find(sdftm == 0); 
xAxisTime = sdftm;
yAxisChannels = 1:32;

for chIdx = yAxisChannels
% for chIdx = 20
    %% Continuous line plots for individual units
    
    % Figure settings
    figure;
    
    % Feature Selectivity - Eye vs ori
    % First 800 ms of condition 5-8,9-12, 13-16, 17-20
    clear cond
    cond.monoc_PO_LE = [8 10 16 18];
    cond.monoc_PO_RE = [5 11 13 19];
    cond.monoc_NOP_LE = [6 12 14 20];
    cond.monoc_NPO_RE = [7 9 15 17];
    
    fields = fieldnames(cond);
    for i = 1:length(fields)
        conditions = cond.(fields{i});
        counter = 0;
        for j = 1:length(conditions)
            trlLength = size(IDX(conditions(j)).LFP_bb,1);
            for trl = 1:trlLength
                counter = counter + 1;
                muaAllTrls.(fields{i})(:,counter) = IDX(conditions(j)).LFP_bb{trl,1}(:,chIdx);
            end
        end
        muaConditionmean = mean(muaAllTrls.(fields{i}),2);
        mua_baselinemean = mean(muaConditionmean(baselineTimeIndex,1));
        mua_blSubAvg.(fields{i}) = muaConditionmean - mua_baselinemean;
        maxVals(i) = max(mua_blSubAvg.(fields{i}));
        minVals(i) = min(mua_blSubAvg.(fields{i}));
    end
    [sortedResponse,sortedResponseIdx] = sort(maxVals,'descend'); 
    [maxVal,maxValLoc] = max(maxVals); % Location: 1="PO_LE", 2="PO_RE",etc.
    [minVal,inValLoc] = min(maxVals); % Location: 1="PO_LE", 2="PO_RE",etc.

    fields = fieldnames(cond);
    preferedStimulus = fields{maxValLoc}(7:end);
    
    % Spectrogram settings
    % 100 samples at 1000Hz is 100ms. 
    % With the number of overlap we can set the temporal resolution for
    % which we can detect signal changes. i.e. noverlap = 99 for win = 100
    % means that the fft will be recalculated for a moving step of every
    % 1ms in the data.
    window = 64; % window, integer | vector | []
        % The optimum window length will depend on your application.
        % If your application is such that you need time domain information 
        % to be more accurate, reduce the size of your windows. 
        % If the application demands frequency domain information to be 
        % more specific, then increase the size of the windows. 
        % If you want to detect "events" in your EEG signal with a 
        % resolution of 10ms, then this should be your window length. 
        % This will give you a frequency resolution of about 100 Hz.
    noverlap = window-1; % Number of overlapped samples, positive integer | []
    nDFT = 513; %  number of DFT points
        % spectrogram uses 1024/2 + 1 = 513 points
    fs = 1000; % sample rate, 1 Hz (default) | positive scalar 
    freqrange = "onesided"; % Frequency range for PSD estimate, "onesided" | "twosided" | "centered
        % one-sided transforms, are the default for real signals
    spectrumtype = "power"; % Power spectrum scaling, "psd" (default) | "power"
        % The spectrogram function has a fourth argument that corresponds 
        % to the segment-by-segment power spectrum or power spectral 
        % density. The ps argument is already squared and includes the 
        % normalization factor. 
        % For one-sided spectrograms of real signals, you still have to 
        % include the extra factor of 2. Set the scaling argument of the 
        % function to "power".
    freqloc = "yaxis"; % Frequency display axis, "xaxis" (default) | "yaxis"
    
    % Plot 4 spectrograms in order of monocular preference in the LFP_bb
    for i = 1: length(sortedResponseIdx)
        subplot(4,2,i)
        clear x
        x = mua_blSubAvg.(fields{sortedResponseIdx(i)}); %input signal, vector
        spectrogram(x,window,noverlap,nDFT,fs,spectrumtype,freqloc)
        xticks([0 200 400 600 800 1000])
        xticklabels({'-200', '0', '200', '400', '600', '800'})
        vline(200)
        ylim([0 100])
        title(fields{sortedResponseIdx(i)},'interpreter','none')
    end
    
    titleText = {strcat('Channel Index = ',string(chIdx));strcat('Preferred Response = ',preferedStimulus)};
    sgtitle(titleText,'interpreter','none')
    
    

    
    %% Physical Alternation
    % Chose condition based on ocular preferences
    % cond.monoc_PO_LE = [8 10 16 18];
    % cond.monoc_PO_RE = [5 11 13 19];
    % cond.monoc_NOP_LE = [6 12 14 20];
    % cond.monoc_NPO_RE = [7 9 15 17];
    if maxValLoc == 1 % PO LE
        nullToPreferred = 17;
        preferredToNull = 18;
    elseif maxValLoc == 2 % PO RE
        nullToPreferred = 20;
        preferredToNull = 19;
    elseif maxValLoc == 3 % NPO_LE
        nullToPreferred = 19;
        preferredToNull = 20;
    elseif maxValLoc == 4 % NPO_RE
        nullToPreferred = 18;
        preferredToNull = 17;
    end
    
    % Null to preferred
    % Average data for individual conditions
        trlLength = size(IDX(nullToPreferred).LFP_bb,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_NtP_adapter(:,counter) = IDX(nullToPreferred).LFP_bb{trl,1}(:,chIdx); %first 800
            muaAllTrls_NtP_flash(:,counter) = IDX(nullToPreferred).LFP_bb{trl,2}(:,chIdx); % second 800
        end
        muaConditionmean_NtP(:,1) = mean(muaAllTrls_NtP_adapter,2);
        muaConditionmean_NtP(:,2) = mean(muaAllTrls_NtP_flash,2);
        mua_baselinemean_NtP = mean(muaConditionmean_NtP(baselineTimeIndex,:));
        mua_blSubAvg_NtP = muaConditionmean_NtP - mua_baselinemean_NtP;
    % Calculate spectrogram from average data
        x_adapter = mua_blSubAvg_NtP(:,1); %input signal, vector, first time period
        [s_NtP_adapter,f,t] = spectrogram(x_adapter,window,noverlap,nDFT,fs,spectrumtype,freqloc);
        x_flash = mua_blSubAvg_NtP(:,2); %input signal, vector, second time period
        [s_NtP_flash,~,~] = spectrogram(x_flash,window,noverlap,nDFT,fs,spectrumtype,freqloc);


    % preferredToNull
    % Average data for individual conditions
        trlLength = size(IDX(preferredToNull).LFP_bb,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_PtN(:,counter,1) = IDX(preferredToNull).LFP_bb{trl,1}(:,chIdx); %first 800
            muaAllTrls_PtN(:,counter,2) = IDX(preferredToNull).LFP_bb{trl,2}(:,chIdx); %second 800
        end
        muaConditionmean_PtN = squeeze(mean(muaAllTrls_PtN,2)); % average across trials
        mua_baselinemean_PtN = mean(muaConditionmean_PtN(baselineTimeIndex,:));
        mua_blSubAvg_PtN = muaConditionmean_PtN - mua_baselinemean_PtN(1);
    % Calculate spectrogram from average data
        x_adapter = mua_blSubAvg_PtN(:,1); %input signal, vector, first time period
        [s_PtN_adapter,~,~] = spectrogram(x_adapter,window,noverlap,nDFT,fs,spectrumtype,freqloc);
        x_flash = mua_blSubAvg_PtN(:,2); %input signal, vector, second time period
        [s_PtN_flash,~,~] = spectrogram(x_flash,window,noverlap,nDFT,fs,spectrumtype,freqloc);


    % Subtract the prefered response from the null response so that any
    % suppression is seen as 
    s_sub_adapter   = abs(s_NtP_adapter - s_PtN_adapter)';% null - prefered
    s_sub_flash     = abs(s_NtP_flash   - s_PtN_flash)';% null - prefered
    
    % Plot Data
    subplot(4,2,5)
    imagesc(t,f,s_sub_adapter)
    set(gca,'YDir','normal')
    xticks([0 .200 .400 .600 .800 1.000])
    xticklabels({'-200', '0', '200', '400', '600', '800'})
    vline(200)
    ylim([0 100])
    vline(.2)
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    e = colorbar;
    e.Label.String = "abs(STFT)";
        title('Pref - null subtraction plot, Physical Alternation, monocular adaptation')
    
    subplot(4,2,6)
    imagesc(t,f,s_sub_flash)
    set(gca,'YDir','normal')
    xticks([0 .200 .400 .600 .800 1.000])
    xticklabels({'-200', '0', '200', '400', '600', '800'})
    vline(200)
    ylim([0 100])
    vline(.2)
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    e = colorbar;
    e.Label.String = "abs(STFT)";
        title('Pref - null subtraction plot, Physical Alternation, monocular flash')
    
    %% Same Stim, different history, preferred vs null (BRFS)
    % Chose condition based on ocular preferences
    % cond.monoc_PO_LE = [8 10 16 18];
    % cond.monoc_PO_RE = [5 11 13 19];
    % cond.monoc_NOP_LE = [6 12 14 20];
    % cond.monoc_NPO_RE = [7 9 15 17];
    if maxValLoc == 1 % PO LE
        nullToPreferred = 9;
        preferredToNull = 10;
    elseif maxValLoc == 2 % PO RE
        nullToPreferred = 12;
        preferredToNull = 11;
    elseif maxValLoc == 3 % NPO_LE
        nullToPreferred = 11;
        preferredToNull = 12;
    elseif maxValLoc == 4 % NPO_RE
        nullToPreferred = 10;
        preferredToNull = 9;
    end
    
    % Null to preferred
    % Average data for individual conditions
        trlLength = size(IDX(nullToPreferred).LFP_bb,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_NtP_adapter(:,counter) = IDX(nullToPreferred).LFP_bb{trl,1}(:,chIdx); %first 800
            muaAllTrls_NtP_flash(:,counter) = IDX(nullToPreferred).LFP_bb{trl,2}(:,chIdx); % second 800
        end
        muaConditionmean_NtP(:,1) = mean(muaAllTrls_NtP_adapter,2);
        muaConditionmean_NtP(:,2) = mean(muaAllTrls_NtP_flash,2);
        mua_baselinemean_NtP = mean(muaConditionmean_NtP(baselineTimeIndex,:));
        mua_blSubAvg_NtP = muaConditionmean_NtP - mua_baselinemean_NtP;
    % Calculate spectrogram from average data
        x_adapter = mua_blSubAvg_NtP(:,1); %input signal, vector, first time period
        [s_NtP_adapter,f,t] = spectrogram(x_adapter,window,noverlap,nDFT,fs,spectrumtype,freqloc);
        x_flash = mua_blSubAvg_NtP(:,2); %input signal, vector, second time period
        [s_NtP_flash,~,~] = spectrogram(x_flash,window,noverlap,nDFT,fs,spectrumtype,freqloc);


    % preferredToNull
    % Average data for individual conditions
        trlLength = size(IDX(preferredToNull).LFP_bb,1);
        counter = 0;
        for trl = 1:trlLength
            counter = counter + 1;
            muaAllTrls_PtN(:,counter,1) = IDX(preferredToNull).LFP_bb{trl,1}(:,chIdx); %first 800
            muaAllTrls_PtN(:,counter,2) = IDX(preferredToNull).LFP_bb{trl,2}(:,chIdx); %second 800
        end
        muaConditionmean_PtN = squeeze(mean(muaAllTrls_PtN,2)); % average across trials
        mua_baselinemean_PtN = mean(muaConditionmean_PtN(baselineTimeIndex,:));
        mua_blSubAvg_PtN = muaConditionmean_PtN - mua_baselinemean_PtN(1);
    % Calculate spectrogram from average data
        x_adapter = mua_blSubAvg_PtN(:,1); %input signal, vector, first time period
        [s_PtN_adapter,~,~] = spectrogram(x_adapter,window,noverlap,nDFT,fs,spectrumtype,freqloc);
        x_flash = mua_blSubAvg_PtN(:,2); %input signal, vector, second time period
        [s_PtN_flash,~,~] = spectrogram(x_flash,window,noverlap,nDFT,fs,spectrumtype,freqloc);


    % Subtract the prefered response from the null response so that any
    % suppression is seen as 
    s_sub_adapter   = abs(s_NtP_adapter - s_PtN_adapter)';% null - prefered
    s_sub_flash     = abs(s_PtN_flash   - s_NtP_flash)';% null - prefered
    
    % Plot Data
    subplot(4,2,7)
    imagesc(t,f,s_sub_adapter)
    set(gca,'YDir','normal')
    xticks([0 .200 .400 .600 .800 1.000])
    xticklabels({'-200', '0', '200', '400', '600', '800'})
    vline(200)
    ylim([0 100])
    vline(.2)
    title('Pref - null subtraction plot, BRFS, monocular adaptation')
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    e = colorbar;
    e.Label.String = "abs(STFT)";
    
    subplot(4,2,8)
    imagesc(t,f,s_sub_flash)
    set(gca,'YDir','normal')
    xticks([0 .200 .400 .600 .800 1.000])
    xticklabels({'-200', '0', '200', '400', '600', '800'})
    vline(200)
    ylim([0 100])
    vline(.2)
    e = colorbar;
    e.Label.String = "abs(STFT)";
    title('Pref - null subtraction plot, BRFS, same stimulus different history')
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
    

    %% Save output
    cd(outputDir)
    pdfOutputName = strcat(fileName(12:end-4),'_individualUnitClassification_spectrogram.pdf');
    exportgraphics(gcf,pdfOutputName,"Append",true)

end