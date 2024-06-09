% fft on raw LFP for each condition
% cd('D:\sortedData_240229')
% load('sortedData_211008_B_bmcBRFS001.mat')
%IDX(1).LFP_bb{trl,event} = 1201x32 double
close all
clearvars -except IDX
penetrationID = '211008_B_bmcBRFS001_probe1';
%EVP
for cond = 1:20
    for trl = 1:size(IDX(cond).LFP_bb,1)
        avgEVPAlongProbeForEachTrl_1(:,trl) = mean(IDX(cond).LFP_bb{trl,1},2);
        avgEVPAlongProbeForEachTrl_2(:,trl) = mean(IDX(cond).LFP_bb{trl,2},2);
    end
    avgEvpAcrossTrls_1 = mean(avgEVPAlongProbeForEachTrl_1,2);
    avgEvpAcrossTrls_2 = mean(avgEVPAlongProbeForEachTrl_2,2);
    concatData = cat(1,avgEvpAcrossTrls_1(1:end-401),avgEvpAcrossTrls_2);
    concatTm = -200:1800;

    %baseline subtract
    blAvg = mean(concatData(100:225,:),1);
    blSubEVP = concatData - blAvg;
    avgEVPtimeXcond(:,cond) = blSubEVP;
    
end
%plot EVP
figure
subplot(3,2,1)
plot(concatTm,avgEVPtimeXcond(:,1)); hold on;
plot(concatTm,avgEVPtimeXcond(:,2)); hold on;
xlim([-200 1800])
ylim([-250 100])
vline(0)
vline(1600)
title('Dioptic')

subplot(3,2,2)
plot(concatTm,avgEVPtimeXcond(:,3)); hold on;
plot(concatTm,avgEVPtimeXcond(:,4)); hold on;
xlim([-200 1800])
ylim([-250 100])
vline(0)
vline(1600)
title('Dichoptic')

subplot(3,2,3)
plot(concatTm,avgEVPtimeXcond(:,5)); hold on;
plot(concatTm,avgEVPtimeXcond(:,6)); hold on;
plot(concatTm,avgEVPtimeXcond(:,7)); hold on;
plot(concatTm,avgEVPtimeXcond(:,8)); hold on;
xlim([-200 1800])
ylim([-200 100])
vline(0)
vline(800)
vline(1600)
title('BRFS-like - dioptic')

subplot(3,2,4)
plot(concatTm,avgEVPtimeXcond(:,9)); hold on;
plot(concatTm,avgEVPtimeXcond(:,10)); hold on;
plot(concatTm,avgEVPtimeXcond(:,11)); hold on;
plot(concatTm,avgEVPtimeXcond(:,12)); hold on;
xlim([-200 1800])
ylim([-200 100])
vline(0)
vline(800)
vline(1600)
title('BRFS - paired percetpual modulations')

subplot(3,2,5)
plot(concatTm,avgEVPtimeXcond(:,13)); hold on;
plot(concatTm,avgEVPtimeXcond(:,14)); hold on;
plot(concatTm,avgEVPtimeXcond(:,15)); hold on;
plot(concatTm,avgEVPtimeXcond(:,16)); hold on;
xlabel("time (ms) from stimulus onset")
ylabel("Voltage (uV)")
xlim([-200 1800])
ylim([-200 100])
vline(0)
vline(800)
vline(1600)
title('Physical alternation -dioptic')

subplot(3,2,6)
plot(concatTm,avgEVPtimeXcond(:,17)); hold on;
plot(concatTm,avgEVPtimeXcond(:,18)); hold on;
plot(concatTm,avgEVPtimeXcond(:,19)); hold on;
plot(concatTm,avgEVPtimeXcond(:,20)); hold on;
xlabel("time (ms) from stimulus onset")
xlim([-200 1800])
ylim([-200 100])
vline(0)
vline(800)
vline(1600)
title('Physical alternation - dichoptic - BRFS perceptual analog')
sgtitle(penetrationID,"Interpreter","none")



%FFT
figure
for cond = 1:20
    for trl = 1:size(IDX(cond).LFP_bb,1)
        testData = IDX(cond).LFP_bb{trl,2}(500:1000,:);
        Y = fft(testData);
        avgPowerAlongProbeForEachTrl(:,trl) = mean(Y(1:100,:),2);
    end

    figure
    subplot(2,1,1)
    avgPowerAcrossTrial = mean(avgPowerAlongProbeForEachTrl,2);
    semilogx(abs(avgPowerAcrossTrial))
    xlim([0 100])
    title("Complex Magnitude of fft Spectrum")
    xlabel("f (Hz)")
    ylabel("|fft(X)|")
    title(strcat('condition_',num2str(cond)),"Interpreter","none")
    
end