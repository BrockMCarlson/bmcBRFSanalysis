function CSD = calcCSD_classic(DAT)
% DAT is in (ms x chan) (i.e. 10000ms x 32ch)


totchan = size(DAT,2);

% electrode contact spacing in mm:
d = .1; % This is the electrode spacing

%padarray would also work here for simplicity.
EVPvak = padarray(DAT,[0 1],NaN,'replicate');
calcCSDonThis = permute(EVPvak,[2 1]); %dim chould be chan x ms
 
% Calculate CSD
%a -- chan above
%b -- f(x) voltage
%c -- chan below
CSD = nan(size(DAT));
for i = 2:totchan+1  %loop through channels of padded LFP (i.e 2-33)
    a = calcCSDonThis(i-1,:); %output is (ms x trial)
    b = calcCSDonThis(i,:); % This is your target channel
    c = calcCSDonThis(i+1,:);
    CSD(:,i-1) = (-1)*(a - 2*b + c)./(d^2);
end

%See page 3 schroeder et al (1998) for more detail.

%Vaknin correction for CSD analysis (Max Sherman [Brown U] implementation) 
% Allows CSD to be performed on all N contacts instead of N-2 contacts 
% See Vaknin G, DiScenna P, Teyler T. A method for calculating current 
% source density (CSD) analysis without resorting to recording sites 
% outside the sampling volume. Journal of neuroscience methods. 
% 1988;24(2):131-135. doi:10.1016/0165-0270(88)90056-8