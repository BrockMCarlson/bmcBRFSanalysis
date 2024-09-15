% I am looking for any units that show significant difference between
% physically identical stimuli. This can either be a modulation between
% condition 18 and 20 or a modulation between condition 22 and 24

% Note, in order to get proper percentages we need to take a few things
% into consideration...
% 1. every unit in IDX is what Bahmani et al (2014) would consider "SM" or
% "sensory stimulation modulation" units. Out of  393 MUA, 365 (92%) were
% visually responsie (VR). Of the VR units, 76% were SM, (275/364 VR). Only
% 24% of the VR units were PM (perceptual stimulus modulation) (88/364).
% However, I believe they should only consider units tuned to eye and
% orientation, as we did in the iScience paper. If they had done this
% comparison, 88/275, they would have found 32% PM. TLDR; If you find ~30%
% units are PM from this sample, this is, unfortunately, not groundbreaking
% news.
% 2. I have several units that will now have conditions (18 & 20 | 22 & 24)
% so we will need to figure out how to account for total percentage units
% from those that were actually shown the conditions upon which to detect
% perceptual modulation.

% Despite these two caveats, I expect to find 75% of units or more are PM.

% We should do this on the sustained period (because the Bahanmi period was
% 500-1000 ms)

ROOTDIR = 'C:\Users\neuropixel\Documents\MATLAB\BMC_AdaptdcosFigures\';
cd(ROOTDIR)
load(strcat(ROOTDIR,filesep,'IDX_iScienceSubmission.mat'))
% IDX breakdown
% Total number of multi-unit channels
number_Total = length(IDX.allV1) + length(ERR);
number_SensoryMod = length(IDX.allV1);

% Now we find the number of units who have the correct conditions presented
% to even evaluate if they show perceptual modulation
count_ConditionsPresent = 0;
count_ConditionsAbsent  = 0;
for i = 1: length(IDX.allV1)
    if (~isnan(IDX.allV1(i).RESP_avg{18}(2)) && ~isnan(IDX.allV1(i).RESP_avg{20}(2))) ...
            || (~isnan(IDX.allV1(i).RESP_avg{22}(2)) && ~isnan(IDX.allV1(i).RESP_avg{24}(2)))
        count_ConditionsPresent = count_ConditionsPresent + 1;
        idxNumberOfAvailableUnitForComparison(count_ConditionsPresent,1) = i;
    else
        count_ConditionsAbsent = count_ConditionsAbsent + 1;
    end
end
% It looks like that for 39 units the full PM condition is present. For
% each of these units, both the prefered and the orthogonal conditions are
% present

% Run T-test on units with condition present to see if they show perceptual
% modulation
%Transient
for j = 1:count_ConditionsPresent
        unitToEvaluate = idxNumberOfAvailableUnitForComparison(j);
        prefFlash = IDX.allV1(unitToEvaluate).RESP_alltrls{18}(1,:);
        nullFlash = IDX.allV1(unitToEvaluate).RESP_alltrls{20}(1,:);
        [h(j,1),p(j,1)] = ttest2(prefFlash,nullFlash);
end
number_PerceptualMod_transient = sum(h);
percent_PMofAvailableUnits_transient = number_PerceptualMod_transient/count_ConditionsPresent


%Sustained
for j = 1:count_ConditionsPresent
        unitToEvaluate = idxNumberOfAvailableUnitForComparison(j);
        prefFlash = IDX.allV1(unitToEvaluate).RESP_alltrls{18}(2,:);
        nullFlash = IDX.allV1(unitToEvaluate).RESP_alltrls{20}(2,:);
        [h(j,1),p(j,1)] = ttest2(prefFlash,nullFlash);
end
number_PerceptualMod_sustained = sum(h);
percent_PMofAvailableUnits_sustained = number_PerceptualMod_sustained/count_ConditionsPresent

%fullTime
for j = 1:count_ConditionsPresent
        unitToEvaluate = idxNumberOfAvailableUnitForComparison(j);
        prefFlash = IDX.allV1(unitToEvaluate).RESP_alltrls{18}(3,:);
        nullFlash = IDX.allV1(unitToEvaluate).RESP_alltrls{20}(3,:);
        [h(j,1),p(j,1)] = ttest2(prefFlash,nullFlash);
end
number_PerceptualMod_fullTime = sum(h);
percent_PMofAvailableUnits_fullTime = number_PerceptualMod_fullTime/count_ConditionsPresent