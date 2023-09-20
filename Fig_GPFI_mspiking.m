%calculate total mitral cell spikes (and average per mc) vs inhibitory
%level and mc spike rate over the trial in 19 ms bins


% mcspikerate = zeros(19,8001); %will cut out the first 200 ms
% mcspikerate_std = zeros(19,8001);
% mcspikerate_1ms = zeros(19,800);
% mcspikerate_1ms_std = zeros(19,800);

clearvars

inptlvls = 0.25:0.25:4.75;%this round, only doing 1.75 to 4.25 in steps of .5
                            %so indices 7, 9, 11, 13, 15, 17
inptindxs = 7:2:19;
ti_lvls = [0 0.9:0.3:2.1]; %3 lvls
mGABA_lvls = 0.18:0.05:1.13; %20 levels, only used 10 (2:2:20)
mGABAindxs = 2:2:20;

load('glomeruli50.mat')

Nmc = 1793;

%first load null trials

sil_mcs = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
sil_mcs_std = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
sil_mcs_temp = zeros(5,numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));

mcedges = 1:4:149;
mcbartics = -1:4:147;
mspike_h_temp2 = zeros(5,numel(mcedges)-1);
mean_hists = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs),numel(mcedges));
std_hists = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs),numel(mcedges));

mSpikes_temp = zeros(5,numel(inptindxs),numel(ti_lvls),numel(mGABAindxs),8001); %5 trials, 5 ti lvls, 9 k lvls, 8001 time steps 

pop_mean_temp = zeros(5,1);
pop_std_temp = zeros(5,1);
pop_mean = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
pop_std = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));

mGABAcounter = 0;
for mGABAindx = mGABAindxs
    mGABAcounter = mGABAcounter+1;
    inptcounter = 0;
    for inptindx = inptindxs
        inptcounter = inptcounter+1;

        for tiindx = [1,2,5]
            if tiindx ==1
                tiindx2 = 0;
            else
                tiindx2 = tiindx;
            end
    %         if tiindx ==5 && inptindx~=7
    %             continue
    %         end
            for n = 1:5  %trial n
                fname = sprintf('LFP50_GPFI_mGABAlvl%g_0gloms_inputlvl%g_tilvl%g_trial%g.mat',...
                    mGABAindx,inptindx,tiindx,n);
                load(fname)
                mSpikes_temp(n,inptcounter,tiindx,mGABAcounter,:) = sum(mSpikeTrain(:,2001:end));
                
                sil_mcs_temp1 = sum(mSpikeTrain(:,2001:end),2);
                sil_mcs_temp(n,inptcounter,tiindx,mGABAcounter) = numel(sil_mcs_temp1(sil_mcs_temp1==0));
        
                %histogram stuff
                mspike_h_temp1 = sum(mSpikeTrain(:,2001:end),2)/0.8; %divide by .8 s to get Hz
                pop_mean_temp(n) = mean(mspike_h_temp1); %average across the population
                pop_std_temp(n) = std(mspike_h_temp1); %standard deviation across the population
                mspike_h_temp2(n,:) = histogram(mspike_h_temp1,mcedges).Values;
            end
    %     mcspikerate(m,:) = mean(mSpikes_temp(:,m,:)); %this is the average number of mc spikes at each time point
    %     mcspikerate_std(m,:) = std(mSpikes_temp(:,m,:));
            sil_mcs(inptcounter,tiindx,mGABAcounter) = mean(sil_mcs_temp(:,inptcounter,tiindx,mGABAcounter));
            sil_mcs_std(inptcounter,tiindx,mGABAcounter) = std(sil_mcs_temp(:,inptcounter,tiindx,mGABAcounter));
        
            mean_hists(inptcounter,tiindx,mGABAcounter,:) = [sil_mcs(inptcounter,tiindx,mGABAcounter) mean(mspike_h_temp2)];
            std_hists(inptcounter,tiindx,mGABAcounter,:) = [sil_mcs_std(inptcounter,tiindx,mGABAcounter) std(mspike_h_temp2)];
            
            pop_mean(inptcounter,tiindx,mGABAcounter) = mean(pop_mean_temp);
            pop_std(inptcounter,tiindx,mGABAcounter) = sqrt(sum(pop_std_temp.^2)/numel(pop_std_temp)); %average the variances then take sqrt to get std again
        end
    end
end

% for i = 1:800
%     %mcspikerate_1ms(:,i) = sum(mcspikerate(:,10*i-9:10*i),2)/0.001; %this sums the spikes in 1 ms bins
%     rel_data = sum(mSpikes_temp(:,:,10*i-9:10*i),3)/0.001;
%     mcspikerate_1ms(:,i) = mean(rel_data);
%     mcspikerate_1ms_std(:,i) = std(rel_data);
% end
%%
% mspikesTI = reshape(mean(sum(mSpikes_temp,5)),[numel(inptindxs),numel(ti_lvls),numel(mGABAindxs)]); %total number of MC spikes for each level of inh
% mspikesTI_std = reshape(std(sum(mSpikes_temp,5)),[numel(inptindxs),numel(ti_lvls),numel(mGABAindxs)]); %associated standard deviations, 19 inh lvls
mspikesTI = squeeze(mean(sum(mSpikes_temp,5))); %total number of MC spikes for each level of inh
mspikesTI_std = squeeze(std(sum(mSpikes_temp,5))); %associated standard deviations, 19 inh lvls
%% Plotting vs max mGABA

% figure('Position',[0, 0,1200,600])
% for tiindx = 1:6
%     subplot(2,3,tiindx)
%     yvalues = inptlvls(inptindxs);
%     xvalues = mGABA_lvls(mGABAindxs);
%     h = heatmap(xvalues,yvalues,squeeze(mspikesTI(:,tiindx,:))/(Nmc*0.8));
%     h.Title=sprintf('Spont Avg MC Spike Rate, g_{Tonic} = %g nS',ti_lvls(tiindx));
%     h.YLabel='Min Mean Input';
%     h.XLabel='GABA Conductance, nS';
%     h.YDisplayData = flip(yvalues);
%     h.ColorLimits = [min(mspikesTI/(Nmc*0.8),[],'all')*1.1,...
%         max(mspikesTI/(Nmc*0.8),[],"all")*1.1];
%     colormap hsv
% end
% figname = sprintf('Fig_Gradedonly_mspiking.png');
% fig = gcf;
% exportgraphics(fig,figname,'Resolution',600)

figure('Position',[0, 0,1120,420])
tiled = tiledlayout(1,3,"TileSpacing","compact");
ticounter = 0;
for tiindx = [1,2,5]
    ticounter = ticounter + 1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(pop_mean(:,tiindx,:)));
    h.Title=['\rm ' sprintf('g_{tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(pop_mean(:,[1,2,5],:),[],'all'),...
        max(pop_mean(:,[1,2,5],:),[],"all")];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap hot
end
h.ColorbarVisible = 'on'; %only on for the last map
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
title(tiled,'Average MC Firing Rate, Graded And Firing Inhibition','fontweight','bold')
figname = sprintf('Fig_GPFI_mspiking.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution',600)


figure('Position',[0, 0,1120,420])
tiled = tiledlayout(1,3,"TileSpacing","compact");
for tiindx = [1,2,5]
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,1-squeeze(sil_mcs(:,tiindx,:))/Nmc);
    h.Title=['\rm ' sprintf('g_{tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    % h.ColorLimits = [min(1 - sil_mcs(:,[1,2,5],:)/Nmc,[],'all'),...
    %     max(1 - sil_mcs(:,[1,2,5],:)/Nmc,[],"all")];
    h.ColorLimits = [0.5,1];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap hot
end
h.ColorbarVisible = 'on'; %only on for the last map
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
title(tiled,'MC Spike Participation, Graded And Firing Inhibition','fontweight','bold')
figname = sprintf('Fig_GPFI_msilent.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution',600)

%% Plot Histograms for each level of each thing
% figure('Position',[0,-50,560,420*2])
% tiled = tiledlayout(1,1,"TileSpacing","tight");
% for tiindx = 3
%     for inptit = 3
%         inptindx = inptindxs(inptit);
%         for mGABAit=1
%             nexttile
%             mGABAindx = mGABAindxs(mGABAit);
%             meanhist_plot = squeeze(mean_hists(inptit,tiindx,mGABAit,:))/Nmc;
%             stdhist_plot = squeeze(std_hists(inptit,tiindx,mGABAit,:))/Nmc;
%             bar(mcbartics,meanhist_plot);
%             hold on
%             er = errorbar(mcbartics,meanhist_plot,stdhist_plot);
%             er.Color = [0 0 0];
%             er.LineStyle = 'none';
%             g_gaba_txt = sprintf('g_{GABA,MC} = %g nS',mGABA_lvls(mGABAindx));
%             text(35,0.5,g_gaba_txt)
%             ylim([0,.62])
%             xlim([-3,80])
%             xticks(0:10:80)
%             if mGABAit ~=5
%                 set(gca,'XTickLabel',[])
%             end
%             hold off
%         end
%     end
% end
% ylabel(tiled,'Fraction of MCs')
% xlabel(tiled,'Firing Rate, Hz')
% title(tiled,sprintf('Input = %g, g_{Tonic} = %g nS',inptlvls(inptindx),ti_lvls(tiindx)),...
%     'fontweight','bold')
% fontsize(18,'points')
% figname = sprintf('Fig_mSpikeRate_Hists_mGABA1_GPFI.png');
% fig = gcf;
% exportgraphics(fig,figname,'Resolution',600)

% load("mSpiking_GPFI_full.mat")
figure('Position',[0,-50,1120,420])
tiled = tiledlayout(1,2,"TileSpacing","compact");
% tiled = tiledlayout(1,2);
tiindx = 5;
inptit = [7,2];
mGABAit = [1,10];
for itit = 1:2
    inptindx = inptindxs(inptit(itit));
    nexttile
    mGABAindx = mGABAindxs(mGABAit(itit));
    meanhist_plot = squeeze(mean_hists(inptit(itit),tiindx,mGABAit(itit),:))/Nmc;
    stdhist_plot = squeeze(std_hists(inptit(itit),tiindx,mGABAit(itit),:))/Nmc;
    bar(mcbartics,meanhist_plot);
    hold on
    er = errorbar(mcbartics,meanhist_plot,stdhist_plot);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    g_gaba_txt = {[sprintf('Min \\mu_r = %g',inptlvls(inptindx))],...
                  [sprintf('g_{GABA,MC} = %g nS',mGABA_lvls(mGABAindx))]};
    text(35,0.375,g_gaba_txt)
    ylim([0,.45])
    xlim([-3,80])
    xticks(0:10:80)
    if itit ~=1
        set(gca,'YTickLabel',[])
    end
    hold off
end


ylabel(tiled,'Fraction of MCs')
xlabel(tiled,'Firing Rate, Hz')
title(tiled,sprintf('MC Firing Rate Distribution, g_{Tonic} = %g nS',ti_lvls(tiindx)),...
    'fontweight','bold')
fontsize(18,'points')
figname = sprintf('Fig_mSpikeRate_Hists_GPFI_full.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution',600)