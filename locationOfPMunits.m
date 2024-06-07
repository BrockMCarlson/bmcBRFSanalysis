%% Setup
clear
close all
% Directories

codeDir = strcat('C:\Users\Brock Carlson\Documents\GitHub\bmcBRFSanalysis');
% codeDir = strcat('C:\Users\neuropixel\Documents\GitHub\bmcBRFSanalysis\publicationFigures');
cd(codeDir)
outDir = 'S:\formattedDataOutputs';
% outDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs\figures_240404';
dataDir = 'S:\';
% dataDir = 'C:\Users\neuropixel\Documents\MATLAB\formattedDataOutputs';
cd(dataDir)
officLamAssign = importLaminarAssignments("C:\Users\Brock Carlson\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);
% officLamAssign = importLaminarAssignments("C:\Users\neuropixel\Box\Manuscripts\Maier\officialLaminarAssignment_bmcBRFS.xlsx", "AnalysisList", [2, Inf]);


%% load DATAOUT
cd(dataDir)
load("DATAOUT_trials.mat")

%% statistics
% Find perceptually modulated units. 
for i = 1:32 % penetrations
    if i == 12 || i == 13 || i == 32
        continue
    end
    trlNum = size(DATAOUT_ps{i},3);
    clear h_trans p_trans h_susta p_susta
    for j = 1:32 % individual channel on each penetration
        clear trans_ps trans_ns susta_ps susta_ns
        for k = 1:trlNum
            trans_ps(k) = mean(DATAOUT_ps{i}(1051:1151,j,k));
            trans_ns(k) = mean(DATAOUT_ns{i}(1051:1151,j,k));
            susta_ps(k) = mean(DATAOUT_ps{i}(1400:1801,j,k));
            susta_ns(k) = mean(DATAOUT_ns{i}(1400:1801,j,k));
        end
        [h_trans(j,1),p_trans(j,1)] = ttest2(trans_ns,trans_ps,"Tail","left");
        [h_susta(j,1),p_susta(j,1)] = ttest2(susta_ns,susta_ps,"Tail","left");
        dOut_trans = meanEffectSize(trans_ps,trans_ns,Effect="cohen");
        d_trans(j,1) = dOut_trans.Effect;
        dOut_susta = meanEffectSize(susta_ps,susta_ns,Effect="cohen");
        d_susta(j,1) = dOut_susta.Effect;
    end
    tuned_trans(:,i) = h_trans;
    tuned_susta(:,i) = h_susta; % output is in 32ch x 31 penetrations
    effect_trans(:,i) = d_trans;
    effect_susta(:,i) = d_susta;
end
figure
spy(tuned_trans)
title('significnt difference in transient')

figure 
spy(tuned_susta)
title('significant different in sustained')

%% Laminar align data
aligned_100_trans_h = zeros(100,size(officLamAssign,1));
aligned_100_susta_h = zeros(100,size(officLamAssign,1));
aligned_100_trans_d = nan(100,size(officLamAssign,1));
aligned_100_susta_d = nan(100,size(officLamAssign,1));
for i = 1:size(officLamAssign,1)-1

    % Calculate V1 ch boundaries
    granBtm = officLamAssign.Probe11stFold4c(i); % channel corresponding to the bottom of layer 4c

    v1Top_new = 50-granBtm+1;
    v1Btm_new = 50+(32-granBtm);
    v1Ch_new = v1Top_new:v1Btm_new;
    
    aligned_100_trans_h(v1Ch_new,i) = tuned_trans(1:32,i);
    aligned_100_susta_h(v1Ch_new,i) = tuned_susta(1:32,i); 
    aligned_100_trans_d(v1Ch_new,i) = effect_trans(1:32,i);
    aligned_100_susta_d(v1Ch_new,i) = effect_susta(1:32,i);
end
figure
subplot(1,2,1)
spy(aligned_100_trans_h)
hline(35)
hline(45)
hline(50)
hline(60)
xlabel('penetration')
ylabel('electrode depth')
ylim([34 61])
title('significnt difference in transient')

subplot(1,2,2)
mes_trans = mean(aligned_100_trans_d,2,"omitmissing");
plot(mes_trans)
xlim([34 61])
ylim([-.1 1])
view([90 90])
vline(35)
vline(45)
vline(50)
vline(60)
xlabel('electrode depth')
ylabel('effect size')
title('effect size by depth in transient window')

figure 
subplot(1,2,1)
spy(aligned_100_susta_h)
hline(35)
hline(45)
hline(50)
hline(60)
xlabel('penetration')
ylabel('electrode depth')
ylim([34 61])
title('significant different in sustained')

subplot(1,2,2)
mes_susta = mean(aligned_100_susta_d,2,"omitmissing");
plot(mes_susta)
xlim([34 61])
ylim([-.1 1])
view([90 90])
vline(35)
vline(45)
vline(50)
vline(60)
xlabel('electrode depth')
ylabel('effect size')
title('effect size by depth in sustained window')

%% grand averages of tuned units
% Laminar align continuous data

data_ps = nan(32,1801,31);
data_ns = nan(32,1801,31);
aligned_100_data_ps = nan(100,1801,31);
aligned_100_data_ns = nan(100,1801,31);
for i = 1:size(officLamAssign,1)-1 % penetrations
    if i == 12 || i == 13 || i == 32
        continue
    end
    for j = 1:32 % individual channel on each penetration
        if tuned_susta(j,i) == 1  % tuned_susta is in 32ch x 31 penetrations
            data_ps(j,:,i) = mean(DATAOUT_ps{i}(:,j,:),3,"omitmissing");
            data_ns(j,:,i) = mean(DATAOUT_ns{i}(:,j,:),3,"omitmissing");
        end
    end

    % Calculate V1 ch boundaries
    granBtm = officLamAssign.Probe11stFold4c(i); % channel corresponding to the bottom of layer 4c

    v1Top_new = 50-granBtm+1;
    v1Btm_new = 50+(32-granBtm);
    v1Ch_new = v1Top_new:v1Btm_new;
    
    aligned_100_data_ps(v1Ch_new,:,i) = data_ps(:,:,i);
    aligned_100_data_ns(v1Ch_new,:,i) = data_ns(:,:,i); 

end


ps_grandAvg = mean(aligned_100_data_ps,3,"omitmissing")';
ns_grandAvg = mean(aligned_100_data_ns,3,"omitmissing")';

% Calculate variance (Using SEM) SEM = std(data)/sqrt(length(data)); 
penetrationNumber = 28;
ps_SEM = std(aligned_100_data_ps,0,3,"omitmissing")'/sqrt(penetrationNumber); 
ns_SEM = std(aligned_100_data_ns,0,3,"omitmissing")'/sqrt(penetrationNumber); 

% Smooth Data
ps_smooth_mean = smoothdata(ps_grandAvg,1,"gaussian",20);
ns_smooth_mean = smoothdata(ns_grandAvg,1,"gaussian",20);
ps_smooth_sem = smoothdata(ps_SEM,1,"gaussian",20);
ns_smooth_sem = smoothdata(ns_SEM,1,"gaussian",20);

% Mean +/- SEM
ps_mps = ps_smooth_mean + ps_smooth_sem; %pref stim -- mean plus sem 1 
ps_mms = ps_smooth_mean - ps_smooth_sem; %pref stim -- mean minus sem 1 
ns_mps = ns_smooth_mean + ns_smooth_sem; %pref stim -- mean plus sem 1 
ns_mms = ns_smooth_mean - ns_smooth_sem; %pref stim -- mean minus sem 1 

%% Now convert grand average array into laminar compartments

ps_S_mean = mean(ps_smooth_mean(:,41:45),2,"omitmissing"); 
ps_G_mean = mean(ps_smooth_mean(:,46:50),2,"omitmissing");
ps_I_mean = mean(ps_smooth_mean(:,51:55),2,"omitmissing");
ns_S_mean = mean(ns_smooth_mean(:,41:45),2,"omitmissing");
ns_G_mean = mean(ns_smooth_mean(:,46:50),2,"omitmissing"); 
ns_I_mean = mean(ns_smooth_mean(:,51:55),2,"omitmissing");

% Mean plus sem
ps_S_mps = mean(ps_mps(:,41:45),2,"omitmissing");
ps_G_mps = mean(ps_mps(:,46:50),2,"omitmissing");
ps_I_mps = mean(ps_mps(:,51:55),2,"omitmissing");
ns_S_mps = mean(ns_mps(:,41:45),2,"omitmissing");
ns_G_mps = mean(ns_mps(:,46:50),2,"omitmissing");
ns_I_mps = mean(ns_mps(:,51:55),2,"omitmissing");

% Mean minus sem
ps_S_mms = mean(ps_mms(:,41:45),2,"omitmissing");
ps_G_mms = mean(ps_mms(:,46:50),2,"omitmissing");
ps_I_mms = mean(ps_mms(:,51:55),2,"omitmissing");
ns_S_mms = mean(ns_mms(:,41:45),2,"omitmissing");
ns_G_mms = mean(ns_mms(:,46:50),2,"omitmissing");
ns_I_mms = mean(ns_mms(:,51:55),2,"omitmissing");



%% Figure generation! 
% tiledLayout plot
tm_full = -200:1600; % 1801 total timepoints
lamCom = figure;
set(gcf,"Position",[1000 123.6667 757.6667 1.1140e+03])
t = tiledlayout(3,1);
nexttile
    plot(tm_full,ps_S_mean,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_S_mps,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_S_mms,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_S_mean,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_S_mps,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_S_mms,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    % ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Supragranular','% Change'})
    % title({strcat('Significant transient modulation?_',string(h_S_trans),'_pVal =',string(p_S_trans)),...
    %    strcat('Significant sustained modulation?_',string(h_S_susta),'_pVal =',string(p_S_susta))},'Interpreter','none')
nexttile
    plot(tm_full,ps_G_mean,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_G_mps,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_G_mms,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_G_mean,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_G_mps,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_G_mms,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    % ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Granular','% Change'})
    % title({strcat('Significant transient modulation?_',string(h_G_trans),'_pVal =',string(p_G_trans)),...
    %    strcat('Significant sustained modulation?_',string(h_G_susta),'_pVal =',string(p_G_susta))},'Interpreter','none')
nexttile
    plot(tm_full,ps_I_mean,'color',[230/255 97/255 1/255],'LineWidth',1.5); hold on
    plot(tm_full,ps_I_mps,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ps_I_mms,'color',[230/255 97/255 1/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_I_mean,'color',[94/255 60/255 153/255],'LineWidth',1.5); 
    plot(tm_full,ns_I_mps,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':'); 
    plot(tm_full,ns_I_mms,'color',[94/255 60/255 153/255],'LineWidth',1,'Linestyle',':');
    % ylim([0 40])
    vline(0)
    vline(800)
    xregion(850,1000)
    xregion(1200,1600)
    ylabel({'Infragranular','% Change'})
    % title({strcat('Significant transient modulation?_',string(h_I_trans),'_pVal =',string(p_I_trans)),...
    %    strcat('Significant sustained modulation?_',string(h_I_susta),'_pVal =',string(p_I_susta))},'Interpreter','none')

titleText = {'Grand average of units that show tuning in the sustained window'};
title(t,titleText,'Interpreter','none')

%save fig
cd(outDir)
figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.png');
saveas(lamCom,figName_lamCom)
figName_lamCom = strcat('MUA_laminarCompartment_','_grandAvg_','.svg');
saveas(lamCom,figName_lamCom)




