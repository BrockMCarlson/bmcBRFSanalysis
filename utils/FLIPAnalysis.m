function [startinglowfreq,endinglowfreq,startinghighfreq,endinghighfreq,goodnessvalue,superficialchannel,deepchannel,highfreqmaxchannel,lowfreqmaxchannel,crossoverchannel] = FLIPAnalysis(nonnormpowmat,laminaraxis,freqaxis,setfreqbool,alphafreqrange, gammafreqrange)
    arguments 
        nonnormpowmat (:, :) double
        laminaraxis (1, :) double
        freqaxis (1, :) double
        setfreqbool (1, 1) double
        % set default values if not passed by user
        alphafreqrange (1, 2) double = [10 19];
        gammafreqrange (1, 2) double = [75 150];
    end
        
% FLIPANALYSIS takes a non-normalized power matrix as an input and computes
% the spectrolaminar cross and returns pertinent values.
% 
% Required Inputs:
%   nonnormpowmat: the non normalized power matrix that you want to analyze
%       using this script. The FLIP algorithm requires the X axis to
%       represent frequency and the Y axis to represent laminar position
%
%   laminaraxis: (1 x n) row vector where n represents the number of probe
%       contacts represented in the power matrix. Example vector for a 32
%       channel probe with 100 micron spacing: 0:0.1:3.1
%       
%   freqaxis: (1 x n) row vector where n represents the number of frequency
%       bins represented in the power matrix. Example vector for a power
%       matrix using frequencies 0 Hz to 150 Hz with frequency bins of 1
%       Hz: 1:150
%
%   setfreqbool: boolean value that, if true, sets frequency window to the
%       optimal values for alpha/beta and gamma found in the Ubiquitous
%       Spectolaminar Motif paper (10-19 Hz for alpha/beta and 75-150 Hz
%       for gamma). If false, the algorithm will look at every possible 
%       frequency range to optimize the alpha/beta and gamma goodness of 
%       fit values (sacrificing efficiency).
%
%   Optional Parameters:
%       alphafreqrange: Gives the user the option to set the alpha/beta
%       frequency range to a desired value. Format is a 1x2 array where
%       element 1,1 is the lower bound and 1,2 is the upper bound. i.e. 
%       [10 19] would utilize 10-19 Hz to determine the alpha/beta region.
%
%       gammafreqrange: Gives the user the option to set the gamma
%       frequency range to a desired value.Format is a 1x2 array where
%       element 1,1 is the lower bound and 1,2 is the upper bound. i.e. 
%       [75 150] would utilize 75-150 Hz to determine the gamma region.
%
%   Example Usage:
%       FLIPANALYSIS(nonnormpowmat,laminaraxis,freqaxis,setfreqbool)
%       runs the FLIP algorithm on the optimal frequency ranges found in 
%       the paper (on the condition that setfreqbool == TRUE).
%
%       FLIPANALYSIS(nonnormpowmat,laminaraxis,freqaxis,setfreqbool)
%       runs the FLIP algorithm on every possible range of frequencies to 
%       determine which returns the best k value value (on the condition 
%       that setfreqbool == FALSE).
%
%       FLIPANALYSIS(nonnormpowmat,laminaraxis,freqaxis,setfreqbool, alphafreqrange, gammafreqrange)
%       allows the user to determine the frequency ranges over which the
%       FLIP algorithm will analyze. FLIP expects these ranges to be input
%       in a 1x2 array format: [lower range, upper range] (on the condition
%       that setfreqbool == TRUE).

% validate user inputs
% % % freqSteps = diff(freqaxis);
% % % if (all(freqSteps ~= freqSteps(1)))
% % %     error("Error using FLIPAnalysis: frequency values must be consistently spaced");
% % % elseif (length(freqaxis) ~= size(nonnormpowmat, 2))
% % %     error("Error using FLIPAnalysis: lengths of power matrix and frequency axis are inconsistent")
% % % elseif freqaxis(1) > freqaxis(length(freqaxis))
% % %     error("Error using FLIPAnalysis: frequency axis must be in ascending order")
% % % else
% % %     freqbinsize = freqSteps(1);
% % % end
% % % 
% % % laminarSteps = diff(laminaraxis);
% % % if (all(laminarSteps ~= laminarSteps(1)))
% % %     error("Error using FLIPAnalysis: laminar values must be consistently spaced");
% % % elseif (length(laminaraxis) ~= size(nonnormpowmat, 1))
% % %     error("Error using FLIPAnalysis: lengths of power matrix and laminar axis are inconsistent")
% % % elseif laminaraxis(1) > laminaraxis(length(laminaraxis))
% % %     error("Error using FLIPAnalysis: laminar axis must be in ascending order")
% % % elseif laminaraxis(1) ~= 0
% % %     error("First element of laminar axis must be 0")
% % % else
% % %     if max(laminaraxis) < 0.6
% % %         error("Error using FLIPAnalysis: laminar axis must be larger than 0.6 mm")
% % %     end
% % %     ind = 1;
% % %     while laminaraxis(ind) < 0.6
% % %         ind = ind + 1;
% % %     end
% % %     minrange = ind;
% % % end
% % % 
% % % if ~(isvector(alphafreqrange) && numel(alphafreqrange) == 2 && size(alphafreqrange, 1) == 1)
% % %     error("Size of alpha/beta frequency range must be a 1x2 row vector")
% % % elseif ~(isvector(gammafreqrange) && numel(gammafreqrange) == 2 && size(gammafreqrange, 1) == 1)
% % %     error("Size of gamma frequency range must be a 1x2 row vector")
% % % end

fillmissing(nonnormpowmat, "constant", 0);
freqSteps = diff(freqaxis);
freqbinsize = freqSteps(1);

ind = 1;
while laminaraxis(ind) < 0.6
    ind = ind + 1;
end
minrange = ind;

% check to determine whether frequency values were set by user
if ~(setfreqbool)
    % Create a data structure to store the results of the VFLIP algorithm
    pooled_results = NaN(1,11,1);
    bki=1;
    nanboolean=~isnan(squeeze(nonnormpowmat));
    startrow= find(nanboolean ==1,1);
    if sum(nanboolean)==0
        startinglowfreq =  nan;
        endinglowfreq =  nan;
        startinghighfreq =  nan;
        endinghighfreq =  nan;
        goodnessvalue = nan;
        startrow=nan;
        superficialchannel=nan;
        deepchannel=nan;
        highfreqmaxchannel=nan;
        lowfreqmaxchannel=nan;
        crossoverchannel=nan;
    else

        powerspectra_image = nonnormpowmat;

        %setting different frequency bins to iterate over
        low_freq = ceil([1,10, 20, 30, 40, 50,60, 70] / freqbinsize);
        gamma_freq = ceil([30, 40, 50, 60, 70, 80,90,100,110,120,130,140,150 ] / freqbinsize);
        pooling_index = 0;

        pooled_results = NaN(1, 11, 1529);

        %iterating through all frequency combinations
        f = waitbar(0, "Finding Optimal Frequencies");
        for low_freq_startindex = 1:length(low_freq)-1
            for low_freq_endindex = low_freq_startindex+1:length(low_freq)

                for gamma_freq_startindex = 1:length(gamma_freq)-1
                    if gamma_freq(gamma_freq_startindex) < low_freq(low_freq_endindex)
                        continue
                    else

                        for gamma_freq_endindex = gamma_freq_startindex+1:length(gamma_freq)
                            alpha = low_freq(low_freq_startindex):low_freq(low_freq_endindex);
                            gamma = gamma_freq(gamma_freq_startindex):gamma_freq(gamma_freq_endindex);
                            % for each range of alpha/beta and gamma windows,
                            % calculate the optimal channel range that produces the highest goodness value
                            [goodness_k_value,upright_boolean,superficialchannel,deepchannel,gammamaxchannel,betamaxchannel,final_crossoverchannel] = evaluateFrequencyRange(powerspectra_image,alpha,gamma,minrange);
                            pooling_index = pooling_index+1;
                            %inputting the goodness value, frequency bin,
                            %information to the pooled_results matrix
                            pooled_results(bki,1,pooling_index) =low_freq(low_freq_startindex)*freqbinsize;
                            pooled_results(bki,2,pooling_index) =low_freq(low_freq_endindex)*freqbinsize;


                            pooled_results(bki,3,pooling_index) =gamma_freq(gamma_freq_startindex)*freqbinsize;
                            pooled_results(bki,4,pooling_index) =gamma_freq(gamma_freq_endindex)*freqbinsize;

                            pooled_results(bki,5,pooling_index) =goodness_k_value;
                            pooled_results(bki,6,pooling_index) =startrow;

                            %setting nans for non-identifiable probes
                            if isempty(superficialchannel)
                                superficialchannel=nan;
                            end
                            if isempty(deepchannel)
                                deepchannel=nan;
                            end

                            if isempty(final_crossoverchannel)
                                final_crossoverchannel = nan;
                            end


                            if isempty(gammamaxchannel)
                                gammamaxchannel = nan;
                            end

                            if isempty(betamaxchannel)
                                betamaxchannel = nan;
                            end
                            %inputting other information about each frequency
                            %pairing

                            pooled_results(bki,7,pooling_index) =superficialchannel;
                            pooled_results(bki,8,pooling_index) =deepchannel;
                            pooled_results(bki,9,pooling_index) =gammamaxchannel;
                            pooled_results(bki,10,pooling_index) =betamaxchannel;
                            pooled_results(bki,11,pooling_index) =final_crossoverchannel;

                        end
                    end
                end
            end
            waitbar(low_freq_startindex/(length(low_freq)-1), f, sprintf('Progress: %d %%', floor((low_freq_startindex/(length(low_freq)-1))*100)))
        end

        probe1results = squeeze(pooled_results(bki,:,:));


        if  isnan(max(abs(probe1results(5,:)))) || max(abs(probe1results(5,:)))==0
            startinglowfreq =  nan;
            endinglowfreq =  nan;
            startinghighfreq =  nan;
            endinghighfreq =  nan;
            goodnessvalue = nan;
            startrow=nan;
            superficialchannel=nan;
            deepchannel=nan;
            highfreqmaxchannel=nan;
            lowfreqmaxchannel=nan;
            crossoverchannel=nan;

        else
            %finidng the optimal frequency pairing that results in the highest
            %magnitude g value.


            indexnum = find(abs(probe1results(5,:))== max(abs(probe1results(5,:))));
            startinglowfreq =  probe1results(1,indexnum);
            endinglowfreq =  probe1results(2,indexnum);
            startinghighfreq =  probe1results(3,indexnum);
            endinghighfreq =  probe1results(4,indexnum);
            goodnessvalue = probe1results(5,indexnum);
            startrow=probe1results(6,indexnum);
            superficialchannel=probe1results(7,indexnum);
            deepchannel=probe1results(8,indexnum);
            highfreqmaxchannel=probe1results(9,indexnum);
            lowfreqmaxchannel=probe1results(10,indexnum);
            crossoverchannel=probe1results(11,indexnum);
        end
    end

else
    alpha = alphafreqrange(1):alphafreqrange(2);
    gamma = gammafreqrange(1):gammafreqrange(2);
    [goodness_k_value,upright_boolean,superficialchannel,deepchannel,gammamaxchannel,betamaxchannel,final_crossoverchannel] = evaluateFrequencyRange(nonnormpowmat,alpha,gamma,minrange);
    startinglowfreq = alphafreqrange(1);
    endinglowfreq = alphafreqrange(2);
    startinghighfreq = gammafreqrange(1);
    endinghighfreq = gammafreqrange(2);
    goodnessvalue = goodness_k_value;
    highfreqmaxchannel = gammamaxchannel;
    lowfreqmaxchannel = betamaxchannel;
    crossoverchannel = final_crossoverchannel;
end

%output
%startinglowfreq is the starting value of the optimal low frequency band
%endinglowfreq is the ending value of the optimal low frequency band
%startinghighfreq is the starting value of the optimal high frequency band
%endinghighfreq is the ending value of the optimal high frequency band
%goodnessvalue is the optimal G goodness value detected by v flip
%superficial channel is the superficial channel of the optimal band
%deep channel is the deep channel of the optimal band? 
%highfreqmaxchannel is the channel where the high freq band is at a max
%lowfreqmaxchannel is the channel where the low freq band is at a max
%crossoverchannel is the channel where the low and high bands cross over
