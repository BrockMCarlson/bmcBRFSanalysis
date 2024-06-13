% fft on raw LFP for each condition
%IDX(1).LFP_bb{trl,event} = 1201x32 double
dataDir = 'D:\sortedData_240229';
outDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir\evpAndFft';



%% Create dir loop
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
for file = 20:size(officLamAssign,1)
    clearvars -except dataDir outDir officLamAssign file
    % load data
    cd(dataDir)
    probeName = char(officLamAssign.SessionProbe(file,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');
    load(fileToLoad)
    % % chStr = officLamAssign.ChtoUse(file,1);
    % % if strcmp(chStr,"1:32")
    % %     yAxisChannels = 1:32;
    % %     yAxisLim =      [1 32];
    % %    contactToPlot = 32:-1:1;
    % % elseif strcmp(chStr,"33:64")
    % %     yAxisChannels = 33:64;
    % %     yAxisLim =      [33 64];
    % %     contactToPlot = 64:-1:33;
    % % end

    %EVP
    for cond = 1:20
        for trl = 1:size(IDX(cond).LFP_bb,1)
            avgEVPAlongProbeForEachTrl_1(:,trl) = mean(IDX(cond).LFP_bb{trl,1},2);
            avgEVPAlongProbeForEachTrl_2(:,trl) = mean(IDX(cond).LFP_bb{trl,2},2);
        end
        if cond <5 
            avgEvpAcrossTrls_1 = mean(avgEVPAlongProbeForEachTrl_1,2);
            blAvg_1 = mean(avgEvpAcrossTrls_1(100:225,:),1);
            blSubEVP_1 = avgEvpAcrossTrls_1 - blAvg_1;
          
            avgEvpAcrossTrls_2 = mean(avgEVPAlongProbeForEachTrl_2,2);
            blSubEVP_2 = avgEvpAcrossTrls_2 - blAvg_1; 
        
            concatData = cat(1,blSubEVP_1(1:end-401),blSubEVP_2);
            concatTm = -200:1800;
            avgEVPtimeXcond(:,cond) = concatData;
        else
            avgEvpAcrossTrls_1 = mean(avgEVPAlongProbeForEachTrl_1,2);
            blAvg_1 = mean(avgEvpAcrossTrls_1(100:225,:),1);
            blSubEVP_1 = avgEvpAcrossTrls_1 - blAvg_1;
          
            avgEvpAcrossTrls_2 = mean(avgEVPAlongProbeForEachTrl_2,2);
            blAvg_2 = mean(avgEvpAcrossTrls_2(100:225,:),1);
            blSubEVP_2 = avgEvpAcrossTrls_2 - blAvg_2; 
        
            concatData = cat(1,blSubEVP_1(1:end-401),blSubEVP_2);
            concatTm = -200:1800;
            avgEVPtimeXcond(:,cond) = concatData;
        end
    end
    
    %plot EVP
    close all
    f = figure;
    set(gcf,"Position",[1 41 2560 1.3273e+03])
    subplot(3,2,1)
    plot(concatTm,avgEVPtimeXcond(:,1)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,2)); hold on;
    xlim([-200 1800])
    % ylim([-250 100])
    vline(0)
    vline(1600)
    trlNum = [length(IDX(1).correctTrialIndex), length(IDX(2).correctTrialIndex)];
    titleText = strcat('Dioptic: trl N = [',num2str(trlNum),']');
    title(titleText)
    legend({IDX(1).conditionString, IDX(2).conditionString},'Location','northoutside')
    
    subplot(3,2,2)
    plot(concatTm,avgEVPtimeXcond(:,3)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,4)); hold on;
    xlim([-200 1800])
    % ylim([-250 100])
    vline(0)
    vline(1600)
    trlNum = [length(IDX(3).correctTrialIndex), length(IDX(4).correctTrialIndex)];
    titleText = strcat('Dioptic: trl N = [',num2str(trlNum),']');
    title(titleText)
    legend({IDX(3).conditionString, IDX(4).conditionString},'Location','northoutside')
    
    subplot(3,2,3)
    plot(concatTm,avgEVPtimeXcond(:,5)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,6)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,7)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,8)); hold on;
    xlim([-200 1800])
    % ylim([-200 100])
    vline(0)
    vline(800)
    vline(1600)
    trlNum = [length(IDX(5).correctTrialIndex), length(IDX(6).correctTrialIndex),...
        length(IDX(7).correctTrialIndex), length(IDX(8).correctTrialIndex)];
    titleText = strcat('Dioptic: trl N = [',num2str(trlNum),']');
    title(titleText)
    legend({IDX(5).conditionString, IDX(6).conditionString,...
        IDX(7).conditionString, IDX(8).conditionString},'Location','northoutside')
    
    subplot(3,2,4)
    plot(concatTm,avgEVPtimeXcond(:,9)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,10)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,11)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,12)); hold on;
    xlim([-200 1800])
    % ylim([-200 100])
    vline(0)
    vline(800)
    vline(1600)
    trlNum = [length(IDX(9).correctTrialIndex), length(IDX(10).correctTrialIndex),...
        length(IDX(11).correctTrialIndex), length(IDX(12).correctTrialIndex)];
    titleText = strcat('Dioptic: trl N = [',num2str(trlNum),']');
    title(titleText)
    legend({IDX(9).conditionString, IDX(10).conditionString,...
        IDX(11).conditionString, IDX(12).conditionString},'Location','northoutside')
    
    subplot(3,2,5)
    plot(concatTm,avgEVPtimeXcond(:,13)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,14)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,15)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,16)); hold on;
    xlabel("time (ms) from stimulus onset")
    ylabel("Voltage (uV)")
    xlim([-200 1800])
    % ylim([-200 100])
    vline(0)
    vline(800)
    vline(1600)
    trlNum = [length(IDX(13).correctTrialIndex), length(IDX(14).correctTrialIndex),...
        length(IDX(15).correctTrialIndex), length(IDX(16).correctTrialIndex)];
    titleText = strcat('Dioptic: trl N = [',num2str(trlNum),']');
    title(titleText)
    legend({IDX(13).conditionString, IDX(14).conditionString,...
        IDX(15).conditionString, IDX(16).conditionString},'Location','northoutside')
    
    subplot(3,2,6)
    plot(concatTm,avgEVPtimeXcond(:,17)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,18)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,19)); hold on;
    plot(concatTm,avgEVPtimeXcond(:,20)); hold on;
    xlabel("time (ms) from stimulus onset")
    xlim([-200 1800])
    % ylim([-200 100])
    vline(0)
    vline(800)
    vline(1600)
    trlNum = [length(IDX(17).correctTrialIndex), length(IDX(18).correctTrialIndex),...
        length(IDX(19).correctTrialIndex), length(IDX(20).correctTrialIndex)];
    titleText = strcat('Dioptic: trl N = [',num2str(trlNum),']');
    title(titleText)
    legend({IDX(17).conditionString, IDX(18).conditionString,...
        IDX(19).conditionString, IDX(20).conditionString},'Location','northoutside')
    
    sgtitle(probeName,"Interpreter","none")
    cd(outDir)
    saveas(f,strcat(probeName,'evp.png'))
    close all
    
    
    %FFT
    for cond = 1:20
        for trl = 1:size(IDX(cond).LFP_bb,1)
            testData = IDX(cond).LFP_bb{trl,2}(500:1000,:);
            Y = fft(testData);
            avgPowerAlongProbeForEachTrl(:,trl) = mean(Y(1:100,:),2);
        end
        avgPowAcrossTrls = mean(avgPowerAlongProbeForEachTrl,2);
        avgFFTfreqXcond(:,cond) = avgPowAcrossTrls;
    end
    
    %plot fft
    f = figure;
    set(gcf,"Position",[1 41 2560 1.3273e+03])
    subplot(3,2,1)
    semilogx(abs(avgFFTfreqXcond(:,1))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,2))); hold on;
        xlim([2 100])
        title("Complex Magnitude of fft Spectrum")
        xlabel("f (Hz)")
        ylabel("|fft(X)|")
        title(strcat('condition_',num2str(cond)),"Interpreter","none")
    title('Dioptic')
    
    subplot(3,2,2)
    semilogx(abs(avgFFTfreqXcond(:,3))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,4))); hold on;
        xlim([2 100])
        title("Complex Magnitude of fft Spectrum")
        xlabel("f (Hz)")
        ylabel("|fft(X)|")
        title(strcat('condition_',num2str(cond)),"Interpreter","none")
    title('Dichoptic')
    
    subplot(3,2,3)
    semilogx(abs(avgFFTfreqXcond(:,5))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,6))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,7))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,8))); hold on;
        xlim([2 100])
        title("Complex Magnitude of fft Spectrum")
        xlabel("f (Hz)")
        ylabel("|fft(X)|")
        title(strcat('condition_',num2str(cond)),"Interpreter","none")
    title('BRFS-like - dioptic')
    
    subplot(3,2,4)
    semilogx(abs(avgFFTfreqXcond(:,9))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,10))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,11))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,12))); hold on;
        xlim([2 100])
        title("Complex Magnitude of fft Spectrum")
        xlabel("f (Hz)")
        ylabel("|fft(X)|")
        title(strcat('condition_',num2str(cond)),"Interpreter","none")
    title('BRFS - paired percetpual modulations')
    
    subplot(3,2,5)
    semilogx(abs(avgFFTfreqXcond(:,13))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,14))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,15))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,16))); hold on;
        xlim([2 100])
        title("Complex Magnitude of fft Spectrum")
        xlabel("f (Hz)")
        ylabel("|fft(X)|")
        title(strcat('condition_',num2str(cond)),"Interpreter","none")
    title('Physical alternation -dioptic')
    
    subplot(3,2,6)
    semilogx(abs(avgFFTfreqXcond(:,17))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,18))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,19))); hold on;
    semilogx(abs(avgFFTfreqXcond(:,20))); hold on;
        xlim([2 100])
        title("Complex Magnitude of fft Spectrum")
        xlabel("f (Hz)")
        ylabel("|fft(X)|")
        title(strcat('condition_',num2str(cond)),"Interpreter","none")
    title('Physical alternation - dichoptic - BRFS perceptual analog')
    sgtitle(probeName,"Interpreter","none")
    cd(outDir)
    saveas(f,strcat(probeName,'fft.png'))
    close all
end
