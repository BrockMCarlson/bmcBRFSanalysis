% example CCH hypothesis figs
% Created using chatGPT 3.5
% January 17th, 2024

% Example MATLAB script for auto-correlation and cross-correlation

% Generate sample data
t = 0:0.01:5; % Time vector
signal1 = sin(2*pi*1*t); % Example signal 1 (sinusoidal)
signal2 = 0.5*sin(2*pi*1*t + pi/4); % Example signal 2 (sinusoidal with phase shift)

% Auto-correlation
auto_corr_signal1 = xcorr(signal1, 'coeff');
auto_corr_signal2 = xcorr(signal2, 'coeff');

% Cross-correlation
cross_corr_signals = xcorr(signal1, signal2, 'coeff');

% Plotting
figure;

subplot(3,1,1);
plot(t, signal1);
title('Signal 1');
xlabel('Time');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, signal2);
title('Signal 2');
xlabel('Time');
ylabel('Amplitude');

subplot(3,1,3);
plot(-length(auto_corr_signal1)/2:length(auto_corr_signal1)/2-1, auto_corr_signal1);
hold on;
plot(-length(auto_corr_signal2)/2:length(auto_corr_signal2)/2-1, auto_corr_signal2);
title('Auto-correlation');
xlabel('Lag');
ylabel('Correlation Coefficient');
legend('Signal 1', 'Signal 2');

figure;
plot(-length(cross_corr_signals)/2:length(cross_corr_signals)/2-1, cross_corr_signals);
title('Cross-correlation');
xlabel('Lag');
ylabel('Correlation Coefficient');

