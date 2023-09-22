%calculates power spectra from the data for tonic inhibition trials, identifies spectral peaks for theta, beta and gamma, and plots Figure 3.3(a) from the dissertation.
%Power spectra calculation is the same as that of Kersen et al., found at https://github.com/dkersen/olfactory-bulb
%(D. E. C. Kersen, G. Tavoni, and V. Balasubramanian. Connectivity and dynamics in the olfactory bulb. PLoS Comput Biol, 18(2):e1009856, 2022. ISSN 1553-7358.
%doi: 10.1371/journal.pcbi.1009856.)


TS = 0.1;
fs = 1/(TS/1000);
fcut = 200;
[b,a] = butter(6,fcut/(fs/2),'low');

freqlength = 100; 
window_size = 4000;
shift = round(window_size/2);
numTrials = 5;

glomtot=89;
gloms = 0;

inptlvls = 0.25:0.25:2;
klvls = 0.1:0.1:0.6; %6 levels
ti_lvls = [0 0.9:0.3:2.1]; %6 lvls

% save_loc = 'LFP50_TonicInh/LFPs_with_TonicInh_contribution/';

%initialize spectrum peak data, numel(lvls) levels
g_peaks = zeros(numel(klvls),numel(ti_lvls)); g_freqs = zeros(numel(klvls),numel(ti_lvls));
b_peaks = zeros(numel(klvls),numel(ti_lvls)); b_freqs = zeros(numel(klvls),numel(ti_lvls)); 
th_peaks = zeros(numel(klvls),numel(ti_lvls)); th_freqs = zeros(numel(klvls),numel(ti_lvls));
g_peaks_sem = zeros(numel(klvls),numel(ti_lvls)); b_peaks_sem = zeros(numel(klvls),numel(ti_lvls)); 
th_peaks_sem = zeros(numel(klvls),numel(ti_lvls));

inptindx = 7;
for kindx = 1:numel(klvls) 
    for tiindx = 1:numel(ti_lvls)
        if tiindx ==1
            tiindx2 = 0;
        else
            tiindx2 = tiindx;
        end
        for trial = 1:numTrials
            fname = append('LFP50_spont_SITI2_0gloms_inputlvl',num2str(inptindx),...
                    '_klvl',num2str(kindx),'_tilvl',num2str(tiindx2),'_trial',num2str(trial),'.mat');
            load(fname)
            
            LFP_GABA_ton = -LFP_GABA_ton; %in experiment, the sign was flipped relative to the 
            %rest of LFP data
    
            % mVoltS = detrend(mVoltS);f
            LFP_tot = LFP_NMDA+LFP_AMPA+LFP_GABA + LFP_GABA_ton;
            y = filtfilt(b,a,LFP_tot);
        
            y1 = y(2001:10000);
            y1 = detrend(y1);
        
            [pyy, f] = pwelch(y1,window_size,shift,[],fs);
        
            
            if trial == 1
                pyy_tot = zeros(numTrials,length(pyy));
            end
            
            pyy_tot(trial,:) = pyy;
        end
    
        f = f';
        LFP_mean = mean(pyy_tot);
        LFP_SEM = std(pyy_tot)/sqrt(numTrials);
        
        %find peak value and frequency for theta, beta, and gamma bands, as
        %well as overall max
        
        g_peak = max(LFP_mean(f>35 & f<80)); %35 - 80 Hz
        gindx = find(LFP_mean==g_peak);
        while g_peak<LFP_mean(gindx-1) %make sure it's not the shoulder of beta peak
            g_peak = max(LFP_mean(gindx+1:33));
            gindx = find(LFP_mean==g_peak);
        end
        g_peak_freq = f(LFP_mean==g_peak);
    
        b_peak = max(LFP_mean(f>14 & f<35)); %14 - 35 Hz
        bindx = find(LFP_mean==b_peak);
        while b_peak<LFP_mean(bindx+1) %make sure it's not shoulder of gamma peak
            b_peak = max(LFP_mean(7:bindx-1));
            bindx = find(LFP_mean==b_peak);
        end
        while b_peak<LFP_mean(bindx-1) %make sure it's not shoulder of theta peak
            b_peak = max(LFP_mean(bindx+1:15));
            bindx = find(LFP_mean==b_peak);
        end
        b_peak_freq = f(LFP_mean==b_peak); 
    
        th_peak = max(LFP_mean(f>2 & f<14)); %2 - 14 Hz
        thindx = find(LFP_mean==th_peak);
        counter=0;
        while (~isempty(th_peak)&&th_peak<LFP_mean(thindx+1)) %make sure it's not shoulder of beta peak
            th_peak = max(LFP_mean(1:thindx-1));
            counter=counter+1;
            if ~isempty(th_peak)
                thindx = find(LFP_mean==th_peak);
            end
            
        end
        if isempty(th_peak)
            th_peak_freq = 0;
            th_peak = 0;
        else
            th_peak_freq = f(LFP_mean==th_peak);
        end
        
    %     max_peak = max(LFP_mean);
    %     max_peak_freq = f(LFP_mean==max_peak);
    
        % errorbar(f(1:freqlength),LFP_mean(1:freqlength),LFP_SEM(1:freqlength))
        % hold on
        % errorbar(f(1:freqlength),LFP_mean_2(1:freqlength),LFP_SEM_2(1:freqlength))
        % xlim([0,100])
        
        figname = sprintf('LFP50_spont_TI_%gglom_inputlvl%g_klvl%g_tilvl%g.png',...
            gloms,inptindx,kindx,tiindx);
    
        errorbar(f(1:freqlength),LFP_mean(1:freqlength),LFP_SEM(1:freqlength))
        xlim([0,100])
        ylim([0,max(LFP_mean(1:freqlength))])
        title(append('Spont Firing, \kappa = ',num2str(klvls(kindx)),...
            ', g_{Tonic} = ',num2str(ti_lvls(tiindx)),...
            ' nS, min input = ',num2str(inptlvls(inptindx))));
        
        peak_vals = {['Theta: ' num2str(th_peak) ' at ' num2str(th_peak_freq) ' Hz'], ...
            ['Beta: ' num2str(b_peak) ' at ' num2str(b_peak_freq) ' Hz'], ...
            ['Gamma: ' num2str(g_peak) ' at ' num2str(g_peak_freq) ' Hz']}; %, ...
    %         ['Max: ' num2str(max_peak) ' at ' num2str(max_peak_freq) ' Hz']};
        text(36,0.85*max(LFP_mean),peak_vals, 'FontSize',12.5);
    
    %     saveas(gcf,append(save_loc,figname))
        saveas(gcf,figname)
    
        g_peaks(kindx,tiindx) = g_peak; g_freqs(kindx,tiindx) = g_peak_freq; 
        b_peaks(kindx,tiindx) = b_peak; b_freqs(kindx,tiindx) = b_peak_freq; 
        th_peaks(kindx,tiindx) = th_peak; th_freqs(kindx,tiindx) = th_peak_freq;
    
        g_peaks_sem(kindx,tiindx) = LFP_SEM(gindx); 
        b_peaks_sem(kindx,tiindx) = LFP_SEM(bindx);
        th_peaks_sem(kindx,tiindx) = LFP_SEM(thindx);
    end
end
%% 


figure('Position',[0, 0, 560, 420])
g_b = g_peaks - b_peaks - g_peaks_sem - b_peaks_sem; th_g = th_peaks - g_peaks; th_b = th_peaks - b_peaks;

yvalues = ti_lvls;
xvalues = klvls;
h = heatmap(xvalues,yvalues,transpose(g_b), Colormap=parula(2));
h.Title=sprintf('Gamma > Beta');
h.YLabel='g_{Tonic} (nS)';
h.XLabel='Spike-Independent Factor, \kappa';
h.YDisplayData = flip(yvalues);
h.ColorLimits = [min(g_b,[],'all'), -min(g_b,[],'all')];
h.CellLabelColor = 'none';
h.ColorbarVisible = 'off';

fontsize(18,'points')

fig = gcf;
figname = sprintf('Fig_TonicInh_g_b.png');
exportgraphics(fig,figname,'Resolution','600')


% for tiindx = 1:numel(ti_lvls)
%     figure
%     hold on
%     errorbar(klvls,g_peaks(:,tiindx),g_peaks_sem(:,tiindx))
%     errorbar(klvls,b_peaks(:,tiindx),b_peaks_sem(:,tiindx))
%     errorbar(klvls,th_peaks(:,tiindx),th_peaks_sem(:,tiindx))
%     legend('gamma', 'beta', 'theta')
%     title(sprintf('Spont, Peak Power vs \\kappa, Min Mean Input = %g, g_{Tonic} = %g nS',...
%         inptlvls(inptindx),ti_lvls(tiindx)))
%     xlabel('\kappa')
%     ylabel('Peak Power, V^2/Hz')
%     ylim([0,1.5e-4])
%     if tiindx ==1
%         tiindx2 = 0;
%     else
%         tiindx2 = tiindx;
%     end
%     saveas(gcf,sprintf('Spont_SITI2_Peaks_inputlvl%g_tilvl%g.png',inptindx,tiindx2))
%     hold off
% end
% 
% for tiindx = 1:numel(ti_lvls)
%     figure
%     plot(klvls,g_freqs(:,tiindx),klvls,b_freqs(:,tiindx),klvls,th_freqs(:,tiindx));
%     legend('gamma', 'beta', 'theta')
%     title(sprintf('Spont, Peak Freqs vs Tonic Inh, Min Mean Input = %g, g_{Tonic} = %g nS',...
%         inptlvls(inptindx),ti_lvls(tiindx)))
%     xlabel('\kappa')
%     ylabel('Peak Freq, Hz')
%     if tiindx ==1
%         tiindx2 = 0;
%     else
%         tiindx2 = tiindx;
%     end
%     saveas(gcf,sprintf('Spont_SITI2_Peak_freqs_inputlvl%g_tilvl%g.png',inptindx,tiindx2))
% end


%% Plot vs TI as well (just the peaks)
% for kindx = 1:numel(klvls)
%     figure
%     hold on
%     errorbar(ti_lvls,g_peaks(kindx,:),g_peaks_sem(kindx,:))
%     errorbar(ti_lvls,b_peaks(kindx,:),b_peaks_sem(kindx,:))
%     errorbar(ti_lvls,th_peaks(kindx,:),th_peaks_sem(kindx,:))
%     legend('gamma', 'beta', 'theta')
%     title(sprintf('Spont, Peak Power vs g_{Tonic}, Min Mean Input = %g, \\kappa = %g',...
%         inptlvls(inptindx),klvls(kindx)))
%     xlabel('g_{Tonic}, nS')
%     ylabel('Peak Power, V^2/Hz')
%     ylim([0,1.5e-4])
%     saveas(gcf,sprintf('Spont_SITI2_Peaks_inputlvl%g_klvl%g.png',inptindx,kindx))
%     hold off
% end

%% Thing I Tried with plotting the peaks: Didn't work out, but it was instructive 
%for using legend in some new ways.
% colors = ["#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30" "#A2142F"];
% ebs = cell(1,3*numel(ti_lvls));
% lgndlbls = strings(1,3*numel(ti_lvls));
% ebsindx = 0;
% figure
% hold on
% for tiindx = 1:numel(ti_lvls)
%     ebsindx = ebsindx + 1;
%     ebs{ebsindx} = errorbar(klvls,g_peaks(:,tiindx),g_peaks_sem(:,tiindx),...
%         'Color',colors(tiindx),'LineStyle',"-");
%     lgndlbls(ebsindx) = sprintf("Gamma, g_{T} = %g",ti_lvls(tiindx));
%     ebsindx = ebsindx + 1;
%     ebs{ebsindx} = errorbar(klvls,b_peaks(:,tiindx),b_peaks_sem(:,tiindx),...
%         'Color',colors(tiindx),'LineStyle',"--");
%     lgndlbls(ebsindx) = sprintf("Beta, g_{T} = %g",ti_lvls(tiindx));
%     ebsindx = ebsindx + 1;
%     ebs{ebsindx} = errorbar(klvls,th_peaks(:,tiindx),th_peaks_sem(:,tiindx),...
%         'Color',colors(tiindx),'LineStyle',":");
%     lgndlbls(ebsindx) = sprintf("Theta, g_{T} = %g",ti_lvls(tiindx));
% end
% legend([ebs{1:3} ebs{4} ebs{7} ebs{10} ebs{13}],{lgndlbls(1), lgndlbls(2),...
%     lgndlbls(3),lgndlbls(4),lgndlbls(7), lgndlbls(10), lgndlbls(13)},"Location","best")
% title(sprintf('Spont, Peak Power vs \\kappa, Min Mean Input = %g',inptlvls(inptindx)))
% xlabel('\kappa')
% ylabel('Peak Power, V^2/Hz')
% ylim([0,1.5e-4])
% saveas(gcf,sprintf('Spont_SITI2_Peaks_inputlvl%g.png',inptindx))
% hold off
