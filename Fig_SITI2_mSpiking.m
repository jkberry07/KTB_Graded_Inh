%calculate total mitral cell spikes (and average per mc) vs inhibitory
%level and mc spike rate over the trial in 19 ms bins


% mcspikerate = zeros(19,8001); %will cut out the first 200 ms
% mcspikerate_std = zeros(19,8001);
% mcspikerate_1ms = zeros(19,800);
% mcspikerate_1ms_std = zeros(19,800);

clearvars

inptlvls = 0.25:0.25:2;
klvls = 0.1:0.1:0.6; %6 levels
ti_lvls = [0 0.9:0.3:2.1]; %6 lvls


load('glomeruli50.mat')

Nmc = 1793;

%first load null trials

sil_mcs = zeros(numel(klvls),numel(ti_lvls));
sil_mcs_std = zeros(numel(klvls),numel(ti_lvls));
sil_mcs_temp = zeros(5,numel(klvls),numel(ti_lvls));

mcedges = 1:4:149;
mcbartics = -1:4:147;
mspike_h_temp2 = zeros(5,numel(mcedges)-1);
mean_hists = zeros(numel(klvls),numel(ti_lvls),numel(mcedges));
std_hists = zeros(numel(klvls),numel(ti_lvls),numel(mcedges));

mSpikes_temp = zeros(5,numel(klvls),numel(ti_lvls),8001); %5 trials, 5 ti lvls, 9 k lvls, 8001 time steps 
inptindx = 7;
pop_mean_temp = zeros(5,1);
pop_std_temp = zeros(5,1);
pop_mean = zeros(numel(klvls),numel(ti_lvls));
pop_std = zeros(numel(klvls),numel(ti_lvls));
for kindx = 1:numel(klvls)
    for tiindx = 1:numel(ti_lvls)
        if tiindx ==1
            tiindx2 = 0;
        else
            tiindx2 = tiindx;
        end
        for n = 1:5  %trial n
            fname = append('LFP50_spont_SITI2_0gloms_inputlvl',num2str(inptindx),...
                '_klvl',num2str(kindx),'_tilvl',num2str(tiindx2),'_trial',num2str(n),'.mat');
            load(fname)
            mSpikes_temp(n,kindx,tiindx,:) = sum(mSpikeTrain(:,2001:end));
            sil_mcs_temp1 = sum(mSpikeTrain(:,2001:end),2);
            sil_mcs_temp(n,kindx,tiindx) = numel(sil_mcs_temp1(sil_mcs_temp1==0));
    
            %histogram stuff
            mspike_h_temp1 = sum(mSpikeTrain(:,2001:end),2)/0.8; %divide by .8 s to get Hz, spike rate for each MC
            pop_mean_temp(n) = mean(mspike_h_temp1); %average across the population
            pop_std_temp(n) = std(mspike_h_temp1); %standard deviation across the population
            mspike_h_temp2(n,:) = histogram(mspike_h_temp1,mcedges).Values;
        end
%     mcspikerate(m,:) = mean(mSpikes_temp(:,m,:)); %this is the average number of mc spikes at each time point
%     mcspikerate_std(m,:) = std(mSpikes_temp(:,m,:));
        sil_mcs(kindx,tiindx) = mean(sil_mcs_temp(:,kindx,tiindx));
        sil_mcs_std(kindx,tiindx) = std(sil_mcs_temp(:,kindx,tiindx));
    
        mean_hists(kindx,tiindx,:) = [sil_mcs(kindx,tiindx) mean(mspike_h_temp2)];
        std_hists(kindx,tiindx,:) = [sil_mcs_std(kindx,tiindx) std(mspike_h_temp2)];

        pop_mean(kindx,tiindx) = mean(pop_mean_temp);
        pop_std(kindx,tiindx) = sqrt(sum(pop_std_temp.^2)/numel(pop_std_temp)); %average the variances then take sqrt to get std again
    end
end


% for i = 1:800
%     %mcspikerate_1ms(:,i) = sum(mcspikerate(:,10*i-9:10*i),2)/0.001; %this sums the spikes in 1 ms bins
%     rel_data = sum(mSpikes_temp(:,:,10*i-9:10*i),3)/0.001;
%     mcspikerate_1ms(:,i) = mean(rel_data);
%     mcspikerate_1ms_std(:,i) = std(rel_data);
% end
%%
mspikesTI = reshape(mean(sum(mSpikes_temp,4)),[numel(klvls),numel(ti_lvls)]); %total number of MC spikes for each level of inh
mspikesTI_std = reshape(std(sum(mSpikes_temp,4)),[numel(klvls),numel(ti_lvls)]); %associated standard deviations, 19 inh lvls

%%
figure
hold on
lgnd = cell(1,numel(ti_lvls));
for tiindx = 1:numel(ti_lvls)
    errorbar(klvls,mspikesTI(:,tiindx)/(Nmc*0.8), mspikesTI_std(:,tiindx)/(Nmc*0.8),'LineWidth',2.1)
    if tiindx ==1
        lgnd{tiindx} = sprintf('g_{Tonic} = %g nS',0);
    else
        lgnd{tiindx} = sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx));
    end
end
ylabel(append('Average MC Spike Rate, Hz'))
ylim([10 18])
% yticks(9:3:24)
xlabel('\kappa')
xlim([0.1 0.6])
xticks(0.1:.1:0.6)
ax = gca;
ax.LineWidth = 1.5;
legend(lgnd,'NumColumns',1,'Location','bestoutside','Box','off')
fontsize(18,'points')
grid on
fig = gcf;
figname = sprintf('Fig_TonicInh_mspiking.png');
exportgraphics(fig,figname,'Resolution',600)
hold off

%this one takes standard dev over the population
figure('Position',[0, 0,280,420])
hold on
lgnd = cell(1,numel(ti_lvls));
for tiindx = 1:numel(ti_lvls)
    errorbar(klvls,pop_mean(:,tiindx), pop_std(:,tiindx)/sqrt(5*Nmc),'LineWidth',2.1) %5*Nmc because 5 trials, Nmc MCs
    if tiindx ==1
        lgnd{tiindx} = sprintf('g_{Tonic} = %g nS',0);
    else
        lgnd{tiindx} = sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx));
    end
end
ylabel(append('Average MC Spike Rate, Hz'))
ylim([10 18])
% yticks(9:3:24)
xlabel('\kappa')
xlim([0.1 0.6])
xticks(0.1:.1:0.6)
ax = gca;
ax.LineWidth = 1.5;
% legend(lgnd,'NumColumns',1,'Location','bestoutside','Box','off')
fontsize(18,'points')
grid on
fig = gcf;
figname = sprintf('Fig_TonicInh_mspiking_v2.png');
exportgraphics(fig,figname,'Resolution',600)
hold off

%plot number of silent MCs
figure('Position',[0, 0,540,420])
hold on
lgnd2 = cell(1,numel(ti_lvls));
for tiindx = 1:numel(ti_lvls)
    errorbar(klvls,1 - sil_mcs(:,tiindx)/Nmc,sil_mcs_std(:,tiindx)/Nmc,'LineWidth',2.1)
    if tiindx ==1
        lgnd{tiindx} = sprintf('g_{Tonic} = %g nS',0);
    else
        lgnd{tiindx} = sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx));
    end
end
ylabel(append('MC Spike Participation'))
xlabel('\kappa')
xlim([0.1 0.6])
xticks(0.1:.1:0.6)
ax = gca;
ax.LineWidth = 1.5;
legend(lgnd,'NumColumns',1,'Location','bestoutside','Box','off')
fontsize(18,'points')
grid on
fig = gcf;
figname = sprintf('Fig_TonicInh_mSilent.png');
exportgraphics(fig,figname,'Resolution',600)
hold off


%%

% for tiindx = 1:numel(ti_lvls)
%     for kindx = 1:numel(klvls)
%         figure
%         figname = append('mSpikeRate_Hist_spont_SITI2_inputlvl',num2str(inptindx),...
%             '_klvl',num2str(kindx),'_tilvl',num2str(tiindx),'.png');
%         meanhist_plot = reshape(mean_hists(kindx,tiindx,:),[1,numel(mcedges)])/Nmc;
%         stdhist_plot = reshape(std_hists(kindx,tiindx,:),[1,numel(mcedges)])/Nmc;
%         bar(mcbartics,meanhist_plot);
%         hold on
%         er = errorbar(mcbartics,meanhist_plot,stdhist_plot);
%         er.Color = [0 0 0];
%         er.LineStyle = 'none';
%         title(append('Spontaneous Spike Rate Histogram, Min Mean input = ',...
%             num2str(inptlvls(inptindx)), ', \kappa = ',num2str(klvls(kindx)),...
%             ', g_{Tonic} =  ', num2str(ti_lvls(tiindx)),' nS'))
%         ylabel(append('Fraction of MCs(',num2str(Nmc),' total MCs)'))
%         xlabel('Spiking Rate, Hz')
%         hold off
%         saveas(gcf,figname)
%     end
% end

%% Plot nsil and avg MC spike rate again, but vs g_tonic

% figure
% hold on
% lgnd = cell(1,numel(klvls));
% for kindx = 1:numel(klvls)
%     errorbar(ti_lvls,mspikesTI(kindx,:)/(Nmc*0.8), mspikesTI_std(kindx,:)/(Nmc*0.8))
%     lgnd{kindx} = sprintf('\\kappa = %g',klvls(kindx));
% end
% title(append('Spontaneous Avg MC Spike Rate vs g_{Tonic}, Min Input = ',...
%         num2str(inptlvls(inptindx))))
% ylabel(append('Spontaneous Avg MC Spike Rate, Hz (',num2str(Nmc),' MCs)'))
% xlabel('g_{Tonic}, nS')
% legend(lgnd,'Location','best')
% figname = sprintf('mSpikes_totalavg_spont_SITI2_vsTI_inputlvl%d.png',inptindx);
% saveas(gcf,figname)
% hold off
% 
% %plot number of silent MCs
% figure
% hold on
% lgnd2 = cell(1,numel(klvls));
% for kindx = 1:numel(klvls)
%     errorbar(ti_lvls,sil_mcs(kindx,:)/Nmc,sil_mcs_std(kindx,:)/Nmc)
%     lgnd{kindx} = sprintf('\\kappa = %g',klvls(kindx));
% end
% title(append('Fraction Silent MCs (Spontaneous) vs g_{Tonic}, Min Input = ',...
%         num2str(inptlvls(inptindx))))
% ylabel(append('Avg fraction silent MCs (',num2str(Nmc),' total MCs)'))
% xlabel('g_{Tonic}')
% legend(lgnd,'Location','best')
% figname = append('mSpikes_nsil_spont_SITI2_vsTI_inputlvl',num2str(inptindx),'.png');
% saveas(gcf,figname)
% hold off




%% Plot spike rate averaged over trials for each level

% t_ms = 1:800;
% for m = 1:19
%     figure
%     errorbar(t_ms,mcspikerate_1ms(m,:),mcspikerate_1ms(m,:),'Marker','o')
%     xlabel('time, ms')
%     ylabel('MC spike rate, Hz')
%     title(append('MC spike rate, tonic inh to GC, g = ',num2str(m),' nS'))
%     fname2 = append('mcSpikeRate_TonInh_lvl',num2str(m),'.png');
%     saveas(gcf,fname2)
% end

