%% Plot ANOVA results as grouped bar plots and save outputs to plotDir

% Prepare data for bar plot
x = {'Cross-Dioptic'; 'Within-Dioptic'}; % Each observation has its own row

% Compute bar heights (means) and errors (standard errors)
y = [nanmean(holder_cross_1); nanmean(holder_same_1)]; % 2x1 vector of bar heights
ymin = y - [nanstd(holder_cross_1) ./ sqrt(sum(~isnan(holder_cross_1))); 
            nanstd(holder_same_1) ./ sqrt(sum(~isnan(holder_same_1)))]; % 2x1
ymax = y + [nanstd(holder_cross_1) ./ sqrt(sum(~isnan(holder_cross_1))); 
            nanstd(holder_same_1) ./ sqrt(sum(~isnan(holder_same_1)))]; % 2x1

% --- Create gramm plot ---

f = figure;
g = gramm('x', x, 'y', y, 'ymin', ymin, 'ymax', ymax); % All inputs are column vectors of the same size
g.geom_bar('width', 0.4); % Remove dodge for simple bars
g.geom_interval('geom', 'errorbar'); % Add error bars
g.set_names('x', 'Comparison Groups', 'y', 'Mean Coherence Difference');
g.set_title('Cross-Dioptic vs Within-Dioptic ANOVA Results');
g.axe_property('Grid', 'on');

% Draw the plot
g.draw();


% % %% Save output to plotDir
% % %save fig
% % answer = questdlg('Would you like to save this figure?', ...
% % 	'Y', ...
% % 	'N');
% % % Handle response
% % switch answer
% %     case 'Yes'
% %         disp('alright, saving figure to plotdir')
% %         sgtitle('Coherence Penetration Average')
% %         cd(plotDir)
% %         saveas(f, 'fig2d_anova_grouped_bar_plots.png');
% %         saveas(f, 'fig2d_anova_grouped_bar_plots.svg');
% %     case 'No'
% %         cd(plotDir)
% %         disp('please see plotdir for last save')
% % end
% % 
