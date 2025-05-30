%% bmcBRFS_fig4g_coherence_tMap

%% plot BRFS
% averageCoherenceMatrix is (ch1 x ch2 x cond x penetration)
coh_diop  = squeeze(averageCoherenceMatrix_diopDichop(:,:,1,:)); % (15 x 15 x N) coh_pref
coh_dichop = squeeze(averageCoherenceMatrix_diopDichop(:,:,2,:)); % (15 x 15 x N) coh_null
[~, ~, ~, stats] = ttest(coh_diop, coh_dichop,'Dim',3); % operates across penetrations
tMap = stats.tstat; % t-values matrix (15 x 15)
figure;
imagesc(tMap);
colormap('jet'); % blue-white-red diverging colormap (you can also use other perceptual maps)
colorbar;
title('Paired t-score map (Dioptic - Dichoptic)');
xlabel('Channel'); ylabel('Channel');
caxis([-max(abs(tMap(:))) max(abs(tMap(:)))]); % symmetric color scaling
[~, p] = ttest(coh_diop, coh_dichop, 'dim', 3);
sigMask = p < (0.05/112);
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

