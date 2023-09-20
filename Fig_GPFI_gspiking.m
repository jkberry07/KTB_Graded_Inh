%Calculate percentage of GCs that spiked above a certain level for various
%values of GC GABA conductancemGABAindx

clearvars

inptlvls = 0.25:0.25:4.75;%this round, only doing 1.75 to 4.25 in steps of .5
                            %so indices 7, 9, 11, 13, 15, 17
ti_lvls = [0 0.9:0.3:2.1]; %6 lvls
mGABA_lvls = 0.18:0.05:1.13; %20 levels, only used 10 (2:2:20)

mGABAindxs = 2:2:20;
inptindxs = 7:2:19;

gSpikesSI = zeros(5,numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
gSpikes_avg = zeros(length(inptindxs),numel(ti_lvls),numel(mGABAindxs)); %avg percentage of GCs that spiked more than 1, 2,
% 3, 4, ... 16 times during the last 800 ms of the trial. 16 spikes/.8 s =
% 20Hz
gSpikes_std = zeros(length(inptindxs),numel(ti_lvls),numel(mGABAindxs)); %associated standard deviations, 16 spike levels,
% 11 conductance levels
NgcSample = 26895; %all GCs
% NgcSample = 5000; %use this for a smaller sample

%first load null trials
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
%             if tiindx ==5 && inptindx~=7
%                 continue
%             end
            gSpikes_temp = zeros(1,5); % 5 trials to be averaged over
            for n = 1:5  %trial n
                fname = sprintf('LFP50_GPFI_mGABAlvl%g_0gloms_inputlvl%g_tilvl%g_trial%g.mat',...
                    mGABAindx,inptindx,tiindx,n);
                load(fname)
                % Ngc = length(gSpikes);
                % gcSample = randperm(Ngc,NgcSample); %select 5000 GCs randomly
                gcSample = 1:NgcSample;
                for i = gcSample
                    gSpikesSI(n,inptcounter,tiindx,mGABAcounter) = gSpikesSI(n,inptcounter,tiindx,mGABAcounter)...
                        + length(gSpikes{i}(gSpikes{i}>2000)); %this will give the total number of GC spikes at each parameter point
                    for k=1 %for each level of spiking
                        if length(gSpikes{i}(gSpikes{i}>2000))>= k %only counting after 200 ms
                            gSpikes_temp(n) = gSpikes_temp(n) + 1;
                        end
                    end
                end
            end
            gSpikes_avg(inptcounter,tiindx,mGABAcounter) = mean(gSpikes_temp);
            gSpikes_std(inptcounter,tiindx,mGABAcounter) = std(gSpikes_temp);
        end
    end
end
save("gSpikes_GPFI_full.mat","gSpikesSI","gSpikes_avg","gSpikes_std","inptindxs","inptlvls",...
    "NgcSample","mGABA_lvls","mGABAindxs","ti_lvls")
%% 
% load("gSpikes_GPFI_full.mat")
% Plot gSpikes_avg at all levels
Ngc = 26895;
gSpikesSIavg = squeeze(mean(gSpikesSI))/(Ngc*0.8); %To get average GC firing rate
gSpikesSIstd = squeeze(std(gSpikesSI))/(Ngc*0.8);

% figure
% hold on
% % lbls = cell(1,numel(ti_lvls));
% figname = sprintf('gSpikes_VDSI1_mGABA1.png');
% for tiindx = 5
%     for inptindx = 7
%         errorbar(mGABA_lvls,squeeze(gSpikesSIavg(inptindx,tiindx,:)/Ngc),...
%         squeeze(gSpikesSIstd(inptindx,tiindx,:)/Ngc))
% %        lbls{tiindx} = sprintf('g_{tonic} = %g',ti_lvls(tiindx));
%     end
% end
% ylabel(append('Average spikes/GC (',num2str(Ngc),' GCs)'))
% xlabel('Maximum MC GABA Conductance (nS)')
% title(sprintf('Spontaneous Average Spikes/GC vs MC GABA Conductance, input = %g, g_{tonic} = %g nS',...
%     inptlvls(inptindx),ti_lvls(tiindx)))
% % legend(lbls,'NumColumns',2,'Location','best')
% saveas(gcf,figname)
% hold off

figure('Position',[0, 0,1120,420])
tiled = tiledlayout(1,3,"TileSpacing","compact");
for tiindx = [1,2,5]
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(gSpikes_avg(:,tiindx,:))/NgcSample);
    h.Title=['\rm ' sprintf('g_{tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(gSpikes_avg(:,[1,2,5],:)/NgcSample,[],'all'),...
        max(gSpikes_avg(:,[1,2,5],:)/NgcSample,[],"all")];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap parula
end
h.ColorbarVisible = 'on'; %only on for the last map
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
title(tiled,'GC Spike Participation, Graded and Firing Inhibition','fontweight','bold')
figname = sprintf('Fig_GPFI_gspiking.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')


figure('Position',[0, 0,1120,420])
tiled = tiledlayout(1,3,"TileSpacing","compact");
for tiindx = [1,2,5]
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(gSpikesSIavg(:,tiindx,:)));
    h.Title=['\rm ' sprintf('g_{tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(gSpikesSIavg(:,[1,2,5],:),[],'all'),...
        max(gSpikesSIavg(:,[1,2,5],:),[],"all")];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap parula
end
h.ColorbarVisible = 'on'; %only on for the last map
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
title(tiled,'Average GC Firing Rate, Graded and Firing Inhibition','fontweight','bold')
figname = sprintf('Fig_GPFI_gfiringrate.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')

% figure
% hold on
% %histdata = zeros(16,19);
% % lbls = cell(1,numel(ti_lvls));
% for i = 5
%     for inptindx = 7
%         errorbar(mGABA_lvls, squeeze(gSpikes_avg(inptindx,i,:))/Ngc,...
%             squeeze(gSpikes_std(inptindx,i,:))/Ngc)
%     %     lbls{i} = sprintf('g_{tonic} = %g',ti_lvls(i));
%     end
% end
% %histdata(16,:) = gSpikes_avg(16,:);
% 
% % legend(lbls,'NumColumns',2,'Location','best')
% ylabel('Fraction of GCs that spiked at least once')
% xlabel('Maximum MC GABA Conductance (nS)')
% title(sprintf('GC Spiking vs MC GABA Conductance, input = %g, g_{tonic} = %g nS',...
%     inptlvls(inptindx),ti_lvls(tiindx)))
% saveas(gcf,sprintf('gSpikes_VDSI1_mGABA1_to20Hz.png'))
% hold off

%% Plot it again, but vs g_tonic

% figure
% hold on
% %histdata = zeros(16,19);
% lbls = cell(1,numel(inptlvls));
% for i = 1:numel(inptlvls)
%     errorbar(ti_lvls, gSpikes_avg(:,i)/Ngc, gSpikes_std(:,i)/Ngc)
%     lbls{i} = sprintf('Input = %g',inptlvls(i));
% end
% %histdata(16,:) = gSpikes_avg(16,:);
% 
% legend(lbls,'NumColumns',2,'Location','best')
% ylabel('% of GCs that spiked at least once')
% xlabel('g_{Tonic}')
% title(append('GC Spiking vs g_{Tonic}'))
% saveas(gcf,sprintf('gSpikes_VDSI1_vsTI_to20Hz.png'))
% hold off


%% 

% 
% for i=1:19
%     figure
%     histdata(:,i) = histdata(:,i)./sum(histdata(:,i));
%     bar(histdata(:,i))
%     xlabel('Fired exactly x times (except last is 16 or more')
%     ylabel('Percentage of firing GCs')
%     fname2 = append('GCfiringbar_SIwithTonicInhlvl_',num2str(i),'.png');
%     saveas(gcf,fname2)
% end

%Make histogram for each level of inhibition, how many GCs spiked how many
%times? How many spiked only once? Twice? etc.


