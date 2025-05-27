%% bmcBRFS_fig4g_coherence_tMap

%% plot BRFS
grandAverageCoherence_BRFS = median(averageCoherenceMatrix_BRFS,4,"omitmissing"); % average across penetration
coh_pref  = squeeze(grandAverageCoherence_BRFS(:,:,1,:)); % (15 x 15 x N)
coh_null = squeeze(grandAverageCoherence_BRFS(:,:,2,:)); % (15 x 15 x N)
[~, ~, ~, stats] = ttest(coh_pref, coh_null, 'dim', 3); % operates across penetrations
tMap = stats.tstat; % t-values matrix (15 x 15)
figure;
imagesc(tMap);
colormap('jet'); % blue-white-red diverging colormap (you can also use other perceptual maps)
colorbar;
title('Paired t-score map (Dioptic - Dichoptic)');
xlabel('Channel'); ylabel('Channel');
caxis([-max(abs(tMap(:))) max(abs(tMap(:)))]); % symmetric color scaling
[~, p] = ttest(coh_diop, coh_dichop, 'dim', 3);
sigMask = p < 0.05;
hold on;
contour(sigMask, 1, 'k', 'LineWidth', 1); % outline significant regions


%% Save output
%save fig
toc
answer = questdlg('Would you like to save this figure?', ...
	'Y', ...
	'N');
% Handle response
switch answer
    case 'Yes'
        disp('alright, saving figure to plotdir')
        sgtitle('Coherence tMap')
        cd(plotDir)
        saveName = strcat('fig4_coherence_tMap.png');
        saveas(f,saveName) 
        saveName = strcat('fig3_coherence.svg');
        saveas(f,saveName) 
    case 'No'
        cd(plotDir)
        disp('please see plotdir for last save')
end

