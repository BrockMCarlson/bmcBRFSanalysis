function [PSD,freq] = calcPSD(DATA,tmOfAnalysis) %input is tm x ch
%% Power spectral density
% Brock Carlson - Maier Lab - 9/10/24

%% PSD
nfft = 512;
PSD = nan(257,size(DATA,2));
n        = size(DATA,1); % Number of data points
Fs  = 1000; % Hz
noverlap = nfft/2;
window   = hanning(nfft); 
nwind    = length(window); 
if n < nwind    % zero-pad x if it has length less than the window length
    x(nwind)=0;  
    n=nwind;
end
k        = fix((n-noverlap)/(nwind-noverlap))-1;	% Number of windows

% FFT
% loop through channels 
for ch = 1:size(DATA,2)
    clear x Spec
    % compute PSD
    x        = DATA(tmOfAnalysis,ch);   
    % fft --> absolute value --> squared.
    Spec    = abs(fft(x,nfft)).^2; % Spec is power.
    % Select first half
    % If the data is real, the PSD will only take the first half of the 
    % FFT results due to the symmetry of the Fourier transform 
    % (for real signals).
    % If the data is complex, the PSD will take the full FFT, and the 
    % number of rows in PSD will be nfft.
    if ~any(imag(x)~=0)   % check if x is real 
        if rem(nfft,2)    % nfft odd - include the nyquist (ie, length 257)
            select = (1:(nfft+1)/2)';
        else
            select = (1:nfft/2+1)';
        end
        Spec = Spec(select);
    else %if x is complex
        error(['Expected inputs are raw LFP voltages (which are not complex)' ...
            'Have you put filtered data into this function?'])
        select = (1:nfft)'; % Take the full spectrum, because it is asymetrical.
    end

    freq = (select - 1)*Fs/nfft;
    PSD(:,ch) = Spec'; 

end

Normalize! You must normalize the output. !~##~!#
