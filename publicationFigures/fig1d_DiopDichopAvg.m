%% fig1
% MUA grand average

%% Setup
disp('start time')
datetime
close
clearvars -except MUA_trials
workingPC = 'home'; % options: 'home', 'office'
if strcmp(workingPC,'home')
    codeDir = 'C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir = 'S:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\Brock Carlson\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
elseif strcmp(workingPC,'office')
    codeDir     = 'C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures';
    dataDir    = 'D:\TrialTriggeredLFPandMUA';
    plotDir = 'C:\Users\neuropixel\Box\Manuscripts\Maier\plotDir';
    officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
end
cd(codeDir)
cd(dataDir)

if ~exist('MUA_trials','var')
    tic
    load('MUA_trials.mat') % format is MUA_trials{penetration,1}{cond,1}{trial,flash}
    toc
end



%%
averageMUAMatrix_BRFSps = nan(2001,15,31);
averageLFPMatrix_BRFSns = nan(2001,15,31);
for penetration = 1:size(MUA_trials,1)
    
    probeName = char(officLamAssign.Session_probe_(penetration,1));
    fileToLoad = strcat('sortedData_',probeName(1:19),'.mat');

    granBtm = officLamAssign.stFold4c(penetration); % channel corresponding to the bottom of layer 4c
    if isnan(granBtm)
        warning(strcat('no sink found on _',fileToLoad))
        continue
    end
    % Calculate V1 ch boundaries
    v1Top = granBtm - 9;
    v1Btm = granBtm + 5;
    v1Ch = v1Top:v1Btm;
    % limit ch to cortex only
    if strcmp(string(officLamAssign.ChToUse(penetration)),"1:32")
        if any(v1Ch > 32) || any(v1Ch < 1)
            warning('skipping session without full column for now')
            disp(fileToLoad)
            continue
        end
    elseif strcmp(string(officLamAssign.ChToUse(penetration)),"33:64")
        if any(v1Ch > 64) || any(v1Ch < 33)
            warning('skipping session without full column for now')
            disp(fileToLoad)
            continue
        end
    end

    %% Obtain Monocular preference
    % bl Sub at average level (better for plotting)
    clear array_ofMonoc1 array_ofMonoc2 array_ofMonoc3 array_ofMonoc4
    
    % Monocular
    monoc_1 = [5, 11, 13, 19]; % PO RightEye
    monoc_2 = [8, 10, 16, 18]; % PO LeftEye
    monoc_3 = [7, 9, 15, 17];  % NPO RightEye
    monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
    
    % convert from cell to double and combine monocular conditions
    count = 0;
    for cond = monoc_1
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc1(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); % 1000 x 32
        end
    end
    clear cond count trl 
    
    count = 0;
    for cond = monoc_2
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc2(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_3
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc3(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    count = 0;
    for cond = monoc_4
        for trl = 1:size(MUA_trials{penetration,1}{cond,1},1)
            count = count + 1;
            array_ofMonoc4(:,:,count) = abs(MUA_trials{penetration,1}{cond,1}{trl,1}(:,v1Ch)); 
        end
    end
    clear cond count trl 
    
    
    % Test for monocular preference with anova
    % The anova function ignores NaN values, <undefined> values, empty 
    % characters, and empty strings in y. If factors or tbl contains NaN or 
    % <undefined> values, or empty characters or strings, the function ignores 
    % the corresponding observations in y. The ANOVA is balanced if each factor 
    % value has the same number of observations after the function disregards 
    % empty or NaN values. Otherwise, the function performs an unbalanced 
    % ANOVA. (We will be doing an unbalanced ANOVA).
    
    % preallocate array with nan - running an unbalanced ANOVA)
    A = [size(array_ofMonoc1,3),...
        size(array_ofMonoc2,3),...
        size(array_ofMonoc3,3),...
        size(array_ofMonoc4,3)];
    maxTrls = max(A);
    minTrls = min(A);
    tuned = nan(length(v1Ch),1);
    prefMonoc = nan(length(v1Ch),1);
    
    % FOR each individual unit 
    for i = 1:length(v1Ch)
        % preallocate y based on trial count
        y = nan(maxTrls,4);
        % create an array of each trial's median response for each monoc 
        % condition (trl x 4 monoc)
        y(1:size(array_ofMonoc1,3),1) = median(squeeze(array_ofMonoc1(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc2,3),2) = median(squeeze(array_ofMonoc2(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc3,3),3) = median(squeeze(array_ofMonoc3(200:450,i,:)),1)'; % median of each trial after stim onset
        y(1:size(array_ofMonoc4,3),4) = median(squeeze(array_ofMonoc4(200:450,i,:)),1)'; % median of each trial after stim onset
        % Now we perform baseline subtraction
        % First we get the blAvg for this contact across all trials
        y_bl(1,1) = median(median(squeeze(array_ofMonoc1(100:200,i,:)),1)');
        y_bl(1,2) = median(median(squeeze(array_ofMonoc2(100:200,i,:)),1)');
        y_bl(1,3) = median(median(squeeze(array_ofMonoc3(100:200,i,:)),1)');
        y_bl(1,4) = median(median(squeeze(array_ofMonoc4(100:200,i,:)),1)');
        y_blSub = y - median(y_bl);
        p = anova1(y_blSub,[],'off');
        % % % if p < .05
            tuned(i,1) = true;
            % Now we find the maximum response
            [M,maxRespIdx] = max(median(y_blSub,1,"omitmissing"));
            prefMonoc(i,1) = maxRespIdx;
            % And assign the null condition
            if maxRespIdx == 1
                nullRespIdx = 4;
            elseif maxRespIdx == 2
                nullRespIdx = 3;
            elseif maxRespIdx == 3
                nullRespIdx = 2;
            elseif maxRespIdx == 4
                nullRespIdx = 1;
            end
            nullMonoc(i,1) = nullRespIdx;
            % % % monoc_1 = [5, 11, 13, 19]; % PO RightEye
            % % % % monoc_2 = [8, 10, 16, 18]; % PO LeftEye
            % % % % monoc_3 = [7, 9, 15, 17];  % NPO RightEye
            % % % % monoc_4 = [6, 12, 14, 20]; % NPO LeftEye
            if prefMonoc(i,1) == 1
                prefCondOnFlash(i,1) = 12;
                nullCondOnFlash(i,1) = 11;
            elseif prefMonoc(i,1) == 2
                prefCondOnFlash(i,1) = 9;
                nullCondOnFlash(i,1) = 10;
            elseif prefMonoc(i,1) == 3
                prefCondOnFlash(i,1) = 10;
                nullCondOnFlash(i,1) = 9;
            elseif prefMonoc(i,1) == 4
                prefCondOnFlash(i,1) = 11;
                nullCondOnFlash(i,1) = 12;
            end
        % % % else
        % % %     tuned(i,1) = false;
        % % %     prefMonoc(i,1) = NaN;
        % % %     nullMonoc(i,1) = NaN;
        % % %     prefCondOnFlash(i,1) = NaN;
        % % %     nullCondOnFlash(i,1) = NaN;
        % % % end

    end
   

   

    %% MUA in percent change
    % BRFS pref vs null
    % Get the number of trials for the chosen condition
    for i = 1:length(v1Ch)
        if ~tuned(i,1)
            continue
        end
        numTrials_ps = size(MUA_trials{penetration,1}{1,1},1);
        MUAflashOut_ps = nan(2001,numTrials_ps);
        for trl = 1:numTrials_ps
            MUAflashOut_ps(:,trl) = MUA_trials{penetration,1}{1,1}{trl,1}(:,v1Ch(i));
        end
        
        % median across trials
        ps_avg = median(MUAflashOut_ps,2,"omitmissing"); % input is (tm,trl)    
        % Calculate as Percent Change
        %              X(t) - avgBl
        % %Ch = 100 * -------------
        %                 avgBl
        psBl = median(ps_avg(1:200,:));
        if psBl == 0
            psBl = .1;
        end
        ps_PercentC = 100*((ps_avg-psBl)./psBl);
        averageMUAMatrix_BRFSps(:,i,penetration) = ps_PercentC;

    
        % Get the number of trials for the chosen condition
        numTrials_ns = size(MUA_trials{penetration,1}{3 ,1},1);
        MUAflashOut_ns = nan(2001,numTrials_ns);
        for trl = 1:numTrials_ns
            MUAflashOut_ns(:,trl) = MUA_trials{penetration,1}{3,1}{trl,1}(:,v1Ch(i));
        end
        % Average across trials and save output
        ns_avg = median(MUAflashOut_ns,2,"omitmissing"); % input is (tm,trl)    
        nsBl = median(ns_avg(1:200,:));
        if nsBl == 0
            nsBl = .1;
        end
        ns_PercentC = 100*((ns_avg-nsBl)./nsBl);  
        averageMUAMatrix_BRFSns(:,i,penetration) = ns_PercentC; % Average across trl. averageMUAMatrix is (tm x ch x x penetration)
    end
    disp(strcat('Done with file number: ',string(penetration)))
end
%%
%%

%% Organize into compartments, median and std across contacts+penetration
% Reshape
useIdx = squeeze(~isnan(averageMUAMatrix_BRFSps(1,1,:))); 
ps_reshaped = reshape(averageMUAMatrix_BRFSps(:,:,useIdx),[2001,405]);
ns_reshaped = reshape(averageMUAMatrix_BRFSns(:,:,useIdx),[2001,405]);

% Average across penetrations
psAvg = smoothdata(mean(ps_reshaped,2,"omitmissing"),"gaussian",20);
nsAvg = smoothdata(mean(ns_reshaped,2,"omitmissing"),"gaussian",20);

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
contactNum = size(ps_reshaped,2);
psSEM = std(ps_reshaped,0,2,"omitmissing")./sqrt(contactNum); 
nsSEM = std(ns_reshaped,0,2,"omitmissing")./sqrt(contactNum); 

% Mean +/- SEM
ps_S_avgPlusSEM = psAvg + psSEM; % Dioptic -- median plus sem
ps_S_avgMinusSEM = psAvg - psSEM; % Dioptic -- median minus sem

ns_S_avgPlusSEM = nsAvg + nsSEM; 
ns_S_avgMinusSEM = nsAvg - nsSEM;  


%%
%% Figure generation! 
% Tiled layout plot
close all
tm_full = (-200:1800)'; % 1801 total timepoints

% Open figure
lamCom = figure;
hold on; % Ensure all plots are displayed on the same figure

% ----------------------
% PLOT FILLED SEM SHADES FOR DIOPTIC AND DICHOPTIC
% ----------------------

% Plot the shaded region for dioptic (psAvg) in **lighter blue**
fill([tm_full; flipud(tm_full)], ...
    [ps_S_avgPlusSEM; flipud(ps_S_avgMinusSEM)], ...
    [0.2, 0.4, 1], 'FaceAlpha', 0.4, 'EdgeColor', [0 0 0.5]); % Light blue shade, darker border

% Plot the shaded region for dichoptic (nsAvg) in **lighter red**
fill([tm_full; flipud(tm_full)], ...
    [ns_S_avgPlusSEM; flipud(ns_S_avgMinusSEM)], ...
    [1, 0.4, 0.4], 'FaceAlpha', 0.4, 'EdgeColor', [0.5 0 0]); % Light red shade, darker border

% ----------------------
% PLOT MAIN LINES FOR DIOPTIC AND DICHOPTIC
% ----------------------
plot(tm_full, psAvg, 'color', [0, 0, 1], 'LineWidth', 1.5); % Bold Blue for dioptic
plot(tm_full, nsAvg, 'color', [1, 0, 0], 'LineWidth', 1.5); % Bold Red for dichoptic

% ----------------------
% PLOT COSMETICS
% ----------------------
ylim([0 45])
xlim([-200 1800])

% Create bold black vertical lines at 0 ms and 1633 ms (you can change the second value)
plot([0 0], ylim, 'k', 'LineWidth', 2); % Vertical line at 0 ms
plot([1633 1633], ylim, 'k', 'LineWidth', 2); % Vertical line at 1633 ms

% Create a bold black horizontal line at 0
plot(xlim, [0 0], 'k', 'LineWidth', 2); % Horizontal line at 0

title('Grand Average')
xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('% change from baseline', 'FontSize', 14, 'FontWeight', 'bold')
box("off")
legend({'Dioptic', 'Dichoptic'}) % Ensure the legend matches the displayed plots

% Set custom x-ticks and y-ticks for better readability
set(gca, 'XTick', [0 800 1600], 'YTick', [0 15 30 45]);

% Make axes bold
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

% Significance asterisks for 100ms bins

bin_width = 50; % ms
time_bins = 0:bin_width:max(tm_full); % Divide time into 100 ms bins
num_bins = length(time_bins) - 1; % Total number of bins for multiple comparisons correction

% Calculate the Bonferroni-adjusted significance threshold
original_threshold = 0.1; % Original p-value threshold for significance
bonferroni_threshold = original_threshold / num_bins; % Adjust for multiple comparisons

% Set y_pos to a constant position above the highest point in the plot
y_pos = max(max(psAvg), max(nsAvg)) - 0.2 * range([psAvg(:); nsAvg(:)]); % Adjust y_pos slightly above max

for i = 1:num_bins
    % Find indices for this 100 ms bin
    bin_indices = find(tm_full >= time_bins(i) & tm_full < time_bins(i+1));
    
    % Extract raw data for the current bin
    ps_bin_data = mean(ps_reshaped(bin_indices, :), 1, 'omitnan'); % Mean across time bin
    ns_bin_data = mean(ns_reshaped(bin_indices, :), 1, 'omitnan');
    
    % Run a t-test (can be replaced with ranksum if data is non-normal)
    [~, p] = ttest2(ps_bin_data, ns_bin_data,'Tail','right');
    pout(i) = p;
    
    % If the p-value is below the Bonferroni-adjusted threshold, plot an asterisk above the two lines
    if p < original_threshold
        % Get the x-position for the middle of the bin (this ensures a scalar)
        x_pos = mean(time_bins(i:i+1)); 
        
        % Plot the asterisk at a consistent vertical position
        text(x_pos, y_pos, '*', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'center');
    end
end




% % %%
% % % ----------------------
% % % PLOT THE DIFFERENCE BETWEEN DICHOPTIC AND DIOPTIC
% % % ----------------------
% % 
% % % Compute the difference (dichoptic - dioptic)
% % difference = nsAvg - psAvg;
% % 
% % % Calculate the standard error of the mean (SEM) for the difference
% % diffSEM = sqrt(psSEM.^2 + nsSEM.^2); % SEM propagation for difference
% % 
% % % 1. Difference plot with Y-axis limits from -3 to 2, and labeling peak values
% % figure;
% % hold on;
% % 
% % % Plot shaded region for ±SEM around the difference
% % fill([tm_full; flipud(tm_full)], ...
% %     [difference + diffSEM; flipud(difference - diffSEM)], ...
% %     [0.8 0.8 0.8], 'EdgeColor', 'none'); % Gray shaded region
% % 
% % % Plot the difference line
% % plot(tm_full, difference, 'k', 'LineWidth', 1.5); % Black line for the difference
% % 
% % % Set Y-axis limits explicitly from -3 to 2
% % ylim([-3 2]);
% % xlim([-200 1800])
% % 
% % % Add labels to indicate peak values (max and min)
% % y_max = max(difference); % Find the maximum value of the difference
% % y_min = min(difference); % Find the minimum value of the difference
% % 
% % % Add text to show max and min values at the right edge of the plot
% % text(tm_full(end) + 50, y_max, sprintf('Peak = %.2f', y_max), 'FontSize', 12, 'VerticalAlignment', 'bottom');
% % text(tm_full(end) + 50, y_min, sprintf('Trough = %.2f', y_min), 'FontSize', 12, 'VerticalAlignment', 'top');
% % 
% % % 2. Plot cosmetics
% % % Significance asterisks for 100ms bins (for the difference plot)
% % bin_width = 50; % ms
% % time_bins = 0:bin_width:max(tm_full); % Divide time into 100 ms bins
% % significance_threshold = 0.05; % p-value threshold for significance
% % 
% % % Set y_pos to a constant position above the highest point in the difference plot
% % y_pos = max(difference) - 0.2 * range(difference); % Adjust y_pos slightly above max of the difference
% % 
% % for i = 1:length(time_bins)-1
% %     % Find indices for this 100 ms bin
% %     bin_indices = find(tm_full >= time_bins(i) & tm_full < time_bins(i+1));
% % 
% %     % Extract the difference data within this bin
% %     diff_data = difference(bin_indices);
% % 
% %     % Run a one-sample t-test (testing if the difference is significantly different from 0)
% %     [~, p] = ttest(diff_data, 0); % t-test against 0
% % 
% %     % If the p-value is below the threshold, plot an asterisk above the difference line
% %     if p < significance_threshold
% %         % Get the x-position for the middle of the bin (this ensures a scalar)
% %         x_pos = mean(time_bins(i:i+1)); 
% % 
% %         % Plot the asterisk at a consistent vertical position
% %         text(x_pos, y_pos, '*', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment', 'center');
% %     end
% % end
% % 
% % 
% % % Create bold black vertical lines at 0 ms and 1600 ms
% % plot([0 0], ylim, 'k', 'LineWidth', 2); % Vertical line at 0 ms
% % plot([1600 1600], ylim, 'k', 'LineWidth', 2); % Vertical line at 1600 ms
% % 
% % % Create a bold black horizontal line at 0
% % plot(xlim, [0 0], 'k', 'LineWidth', 2); % Horizontal line at 0
% % 
% % % Set custom x-ticks and y-ticks for better readability
% % set(gca, 'XTick', [0 800 1600], 'YTick', [-3 0 2]);
% % 
% % % Set font size and make the axes bold
% % set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% % xlabel('Time (ms)', 'FontSize', 14, 'FontWeight', 'bold');
% % ylabel('Difference (Dichoptic - Dioptic)', 'FontSize', 14, 'FontWeight', 'bold');
% % 
% % % Add title
% % title('Difference Between Dichoptic and Dioptic Traces', 'FontSize', 16, 'FontWeight', 'bold');
% % 
% % hold off;


%% save fig
answer = questdlg('Would you like to save these figures?', ...
    'Save Figures', ...
    'Yes', 'No', 'No');

% Handle response
switch answer
    case 'Yes'
        disp('Saving figures to plotDir...')
        cd(plotDir)
        
        % Save diopticVsDichoptic figure
        figName_lamCom_png = strcat('fig1d','_DiopticVsDichoptic_','.png');
        saveas(lamCom, figName_lamCom_png)
        figName_lamCom_svg = strcat('fig1d','_DiopticVsDichoptic_','.svg');
        saveas(lamCom, figName_lamCom_svg)
        
        % % % Save difference figure
        % % figName_difference_png = strcat('fig1d','_Difference_','.png');
        % % saveas(gcf, figName_difference_png) % Assuming the difference plot is the current figure
        % % figName_difference_svg = strcat('fig1d','_Difference_','.svg');
        % % saveas(gcf, figName_difference_svg)
        
    case 'No'
        disp(plotDir)
        disp(codeDir)
        disp('Please see plotDir for last saved versions.')
end

