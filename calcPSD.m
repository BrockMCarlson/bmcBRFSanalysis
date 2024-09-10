function [freq_vector,powerPermute,power_norm] = calcPSD(DATA,Fs,nfft,noverlap)
%% Power spectral density
% Brock Carlson - Maier Lab - 9/10/24

%% PSD
if isempty(vargin)
    Fs  = 1000; % Hz
    nfft = 512;
    noverlap = nfft/2;
end


% FFT
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