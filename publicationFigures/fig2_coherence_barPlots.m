%% Plot ANOVA results as grouped bar plots and save outputs to plotDir

% Assuming aov_cross_1, aov_same_1, aov_cross_2, and aov_same_2 contain the ANOVA results

% Prepare data for bar plot
categories = {'GxS', 'IxS', 'IxG', 'SxS', 'GxG', 'IxI'};

% Grouped data for Dioptic vs Dichoptic
bars_dioptic = [nanmean(holder_cross_1); nanmean(holder_same_1)];
errors_dioptic = [nanstd(holder_cross_1) ./ sqrt(sum(~isnan(holder_cross_1))); 
                  nanstd(holder_same_1) ./ sqrt(sum(~isnan(holder_same_1)))];

% Grouped data for BRFS
bars_brfs = [nanmean(holder_cross_2); nanmean(holder_same_2)];
errors_brfs = [nanstd(holder_cross_2) ./ sqrt(sum(~isnan(holder_cross_2))); 
               nanstd(holder_same_2) ./ sqrt(sum(~isnan(holder_same_2)))];

% Combine into one larger grouped dataset
all_bars = [bars_dioptic; bars_brfs];
all_errors = [errors_dioptic; errors_brfs];

% Group names
group_labels = {'Cross-Dioptic', 'Within-Dioptic', 'Cross-BRFS', 'Within-BRFS'};

% Ensure all inputs have the same length
num_conditions = size(all_bars, 1);
num_categories = size(all_bars, 2);

% Use gramm toolbox to create grouped bar plot
f = figure;
g = gramm('x', repmat(categories, 1, num_conditions), ...
          'y', all_bars(:), ...
          'ymin', all_bars(:) - all_errors(:), ...
          'ymax', all_bars(:) + all_errors(:), ...
          'color', repelem(group_labels, num_categories));

% Use gramm's built-in grouping
g.geom_bar('width', 0.4, 'dodge', 0.7);
g.geom_interval('geom', 'errorbar', 'dodge', 0.7);
g.set_names('x', 'Comparison Groups', 'y', 'Mean Coherence Difference', 'color', 'Condition');
g.set_title('ANOVA Results as Grouped Bar Plots');
g.axe_property('Grid', 'on');

% Draw the plot
g.draw();

%% Save output to plotDir
% Get the current directory
originalDir = pwd;

try
    % Change to the plotDir directory (as defined in the original script)
    cd(plotDir);
    
    % Save the figure in multiple formats
    saveas(f, 'anova_grouped_bar_plots.png');
    saveas(f, 'anova_grouped_bar_plots.svg');
    
    disp('Grouped bar plots successfully saved to plotDir.');
catch ME
    % If an error occurs, display a message and return to the original directory
    disp('An error occurred while saving the bar plots.');
    disp(ME.message);
end

% Return to the original directory
cd(originalDir);