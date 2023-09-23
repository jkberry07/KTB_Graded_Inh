%Calculates GC participation and plots Figure 3.3(c) from the dissertation

clearvars

inptlvls = 0.25:0.25:2;
klvls = 0.1:0.1:0.6; %6 levels
ti_lvls = [0 0.9:0.3:2.1]; %6 lvls

gSpikesSI = zeros(5,numel(klvls),numel(ti_lvls));
gSpikes_avg = zeros(numel(ti_lvls),length(klvls)); %avg percentage of GCs that spiked more than 1, 2,
% 3, 4, ... 16 times during the last 800 ms of the trial. 16 spikes/.8 s =
% 20Hz
gSpikes_std = zeros(numel(ti_lvls),length(klvls)); %associated standard deviations, 16 spike levels,
% 11 conductance levels

%first load null trials

inptindx = 7;
for tiindx = 1:numel(ti_lvls)
    if tiindx ==1
            tiindx2 = 0;
    else
        tiindx2 = tiindx;
    end
    for kindx = 1:numel(klvls) %input lvl
        gSpikes_temp = zeros(1,5); % 5 trials to be averaged over
        for n = 1:5  %trial n
            fname = append('LFP50_spont_SITI2_0gloms_inputlvl',num2str(inptindx),...
                    '_klvl',num2str(kindx),'_tilvl',num2str(tiindx2),'_trial',num2str(n),'.mat');
            load(fname)
            Ngc = length(gSpikes);
    
            for i=1:Ngc
                gSpikesSI(n,kindx,tiindx) = gSpikesSI(n,kindx,tiindx) + length(gSpikes{i}(gSpikes{i}>2000));
                for k=1 %for each level of spiking
                    if length(gSpikes{i}(gSpikes{i}>2000))>= k %only counting after 200 ms
                        gSpikes_temp(n) = gSpikes_temp(n) + 1;
                    end
                end
            end
        end
        gSpikes_avg(tiindx,kindx) = mean(gSpikes_temp);
        gSpikes_std(tiindx,kindx) = std(gSpikes_temp);
    end
end

%% 

% Plot gSpikes_avg at all levels
gSpikesSIavg = reshape(mean(gSpikesSI),[numel(klvls),numel(ti_lvls)]);
gSpikesSIstd = reshape(std(gSpikesSI),[numel(klvls),numel(ti_lvls)]);

figure
hold on
lbls = cell(1,numel(ti_lvls));
figname = sprintf('gSpikes_spont_SITI2_inputlvl%g_avg.png',inptindx);
for tiindx = 1:numel(ti_lvls)
    errorbar(klvls,gSpikesSIavg(:,tiindx)/Ngc,gSpikesSIstd(:,tiindx)/Ngc)
    lbls{tiindx} = sprintf('g_{tonic} = %g',ti_lvls(tiindx));
end
ylabel(append('Average spikes/GC (',num2str(Ngc),' GCs)'))
xlabel('\kappa')
title(append('Spontaneous Average Spikes/GC vs \kappa'))
legend(lbls,'NumColumns',2,'Location','best')
saveas(gcf,figname)
hold off


figure('Position',[0, 0,560,420])

%     subplot(2,3,tiindx)
yvalues = ti_lvls;
xvalues = klvls;
h = heatmap(xvalues,yvalues,squeeze(gSpikes_avg/Ngc));
h.Title=sprintf('GC Spike Participation');
h.YLabel='g_{Tonic}, nS';
h.XLabel='Spike-Independent Factor, \kappa';
h.YDisplayData = flip(yvalues);
h.ColorLimits = [0,1];
fontsize(18,'points')
% colormap hsv
figname = sprintf('Fig_TonicInh_gspiking_hmap.png');
saveas(gcf,figname)



% figure('Position',[0, 0,560,420])
% hold on
% lbls = cell(1,numel(klvls));
% for i = 1:numel(klvls)
%     errorbar(ti_lvls, gSpikes_avg(:,i)/Ngc, gSpikes_std(:,i)/Ngc,'LineWidth',2.1)
%     lbls{i} = sprintf('\\kappa = %g',klvls(i));
% end
% legend(lbls,'NumColumns',2,'Location','northeast')
% ylabel('GC Spike Participation')
% xlabel('g_{Tonic}')
% xlim([0 2.1])
% xticks([0 0.3 0.6 0.9 1.2 1.5 1.8 2.1])
% ylim([0,1.2])
% yticks([0 .2 .4 .6 .8 1])
% grid on
% fontsize(18,'points')
% ax = gca;
% ax.LineWidth = 1.5;
% fig = gcf;
% figname = sprintf('Fig_TonicInh_gspiking.png');
% exportgraphics(fig,figname,'Resolution',450)
% hold off

figure('Position',[0, 0,280,420])
hold on
lbls = cell(1,numel(ti_lvls));
for i = 1:numel(ti_lvls)
    errorbar(klvls, gSpikes_avg(i,:)/Ngc, gSpikes_std(i,:)/Ngc,'LineWidth',2.1)
    lbls{i} = sprintf('g_{Tonic} = %g nS',ti_lvls(i));
end
% legend(lbls,'NumColumns',1,'Location','bestoutside','Box','off')
ylabel('GC Spike Participation')
xlabel('\kappa')
xlim([0.1 .6])
xticks(.1:.1:.6)
ylim([0,1])
yticks([0 .2 .4 .6 .8 1])
grid on
fontsize(18,'points')
ax = gca;
ax.LineWidth = 1.5;
fig = gcf;
figname = sprintf('Fig_TonicInh_gspiking.png');
exportgraphics(fig,figname,'Resolution',600)
hold off

%% Plot it again, but vs g_tonic

figure
hold on
%histdata = zeros(16,19);
lbls = cell(1,numel(klvls));
for i = 1:numel(klvls)
    errorbar(ti_lvls, gSpikes_avg(:,i)/Ngc, gSpikes_std(:,i)/Ngc)
    lbls{i} = sprintf('\\kappa = %g',klvls(i));
end
%histdata(16,:) = gSpikes_avg(16,:);

legend(lbls,'NumColumns',2,'Location','best')
ylabel('% of GCs that spiked at least once')
xlabel('g_{Tonic}')
title(append('GC Spiking vs g_{Tonic}, Spontaneous Min Mean Input = ',num2str(inptlvls(inptindx))))
saveas(gcf,sprintf('gSpikes_spont_SITI2_vsTI_inptlvl%g_to20Hz.png',inptindx))
hold off


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


