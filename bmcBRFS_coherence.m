% bmcBRFS_coherence.m
% Brock Carlson - 1/9/24
% The goal of this code is to analyze the bmcBRFS data in the same manner
% as Fries, Roelfsema, Engel, Konig, & Singer (1997) PNAS, which suggests
% that differential changes in oscilitory patterning at the early stages of
% visual processing determines which signals are perceived, as opposed to
% neural discharge rates.
% CCH - cross-correlation histogram
% RMA - relative modulation amplitude
% STA - spike triggered averages
% SFC - spike-field coherence.
% We aim to first perform this analysis on individual units, individual
% sessions, laminar compartments, and then full data averages.


%% Setup
clear
% Directories
path1 = strcat(getenv('HOMEDRIVE'),getenv("HOMEPATH"));
path2 = 'Documents\GitHub\bmcBRFSanalysis';
codeDir = strcat(path1,filesep,path2);
cd(codeDir)
path3 = 'Documents\MATLAB\formattedDataOutputs\figures_231201';
plotDir = strcat(path1,filesep,path3);

%Import laminar assignments
officLamAssign = importLaminarAssignments("C:\Users\neuropixel\..." + ...
    "Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx",...
    "Sheet1", [2, Inf]);

%% Calculate CCH (cross-correlation histogram) from MUA
% Plotted as +/- 80ms vs RMA
% Detailed methods available in detail at Konig (1994) J Neurosci. Methods
% 1. Calculate correlogram
% 2. Fit a damped cosine wave (Gabor function) to the correlogram
% 3. Ensure that the function accounts for...
%   3a. at least 15% of the data variance 
%   3b. the z-scores of significant peaks must be >2

% Calculate RMA (Response modulation amplitude - a percentage) from MUA
% The strength of synchronization and the regularity of oscillations
% RMA are calculated for the central and the first satellite peak,
% respectively.
% RMA = amplitude of respective peak (measured from the offset caused by
% accidntial coincidences) divided by the offset and multiplied by 100.

%% Calculate STA (Spike-triggered averages) from LFP
% Plotted as +/- 128ms vs LFP uV
% 1. band-pass filter raw data between 1-100 Hz
% 2. average LFPs within a window of +/- 128ms centered on each triggered
% spike.
% Methods detailed in Gray & Singer (1989) PNAS


%% Calculate SFC (spike-field coherence) 
% Plotted as freq from 0 to 102 Hz in 3.9Hz bins vs 0-0.1 SFC
% You should see peak between 39 and 63 Hz
% For each of the LFP segments used for computation of STAs...
% 1. calculate the power spectrum
% 2. Average these spectra (obtaining the spike-triggered power spectrum)
% 3. SFC (computed as a ratio) = power spectrum of STA / spike-triggered
% power spectrum.
