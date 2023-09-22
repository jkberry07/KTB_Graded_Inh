%calculates power spectra from the data, identifies spectral peaks, and plots Figures 3.9 - 3.11, B.1 and B.2 from the dissertation.
%Power spectra calculation is the same as that of Kersen et al., found at https://github.com/dkersen/olfactory-bulb
%(D. E. C. Kersen, G. Tavoni, and V. Balasubramanian. Connectivity and dynamics in the olfactory bulb. PLoS Comput Biol, 18(2):e1009856, 2022. ISSN 1553-7358.
%doi: 10.1371/journal.pcbi.1009856.)


%IMPORTANT: Here, gamma variable describes 19 - 35 Hz, beta as 9 - 18 Hz, theta as 4 - 8 Hz (so they are not the gamma, beta, 
%and theta bands, they are the three oscillatory regions observed in these experiments)

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

inptlvls = 0.25:0.25:4.75;
inptindxs = 7:2:19;
ti_lvls = [0 0.9:0.3:2.1]; % only looking at indexes 1, 2, and 5 (0, .9, and 1.8)
ti_indxs = [1,2,5];
mGABA_lvls = 0.18:0.05:1.13; %20 levels, only looking at 10 of them
mGABAindxs = 2:2:20;

% save_loc = 'LFP50_TonicInh/LFPs_with_TonicInh_contribution/';

%initialize spectrum peak data, numel(lvls) levels
g_peaks = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs)); 
g_powers = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_low_powers = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_high_powers = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_freqs = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
b_peaks = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
b_powers = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
b_freqs = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs)); 
th_peaks = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
th_powers = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
th_freqs = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_peaks_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_powers_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_low_powers_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
g_high_powers_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
b_peaks_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs)); 
b_powers_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
th_peaks_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
th_powers_sem = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs));
secondary_peaks = zeros(size(g_peaks));
secondary_freqs = zeros(size(g_peaks));
secondary_peaks_sem = zeros(size(g_peaks));

LFP_means = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs),2049);
LFP_SEMs = zeros(numel(inptindxs),numel(ti_lvls),numel(mGABAindxs),2049);

mGABAcounter = 0;
for mGABAindx = mGABAindxs
    mGABAcounter = mGABAcounter+1;
    inptcounter = 0;
    for inptindx = inptindxs
        inptcounter = inptcounter+1;
        ticounter = 0;
        for tiindx = ti_indxs
            ticounter = ticounter + 1;
            if tiindx ==1
                tiindx2 = 0;
            else
                tiindx2 = tiindx;
            end
%             if tiindx ==5 && inptindx~=7
%                 continue
%             end
            for trial = 1:numTrials
                fname = append('LFP50_GPFI_mGABAlvl',num2str(mGABAindx),'_0gloms_inputlvl',num2str(inptindx),...
                        '_tilvl',num2str(tiindx),'_trial',num2str(trial),'.mat');
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
            
            g_peak = max(LFP_mean(f>19 & f<35)); %35 - 80 Hz
            gindx = find(LFP_mean==g_peak);
            while g_peak<LFP_mean(gindx-1) %make sure it's not the shoulder of beta peak
                g_peak = max(LFP_mean(gindx+1:15));
                if isempty(g_peak) %this means there was no distinct beta peak
                    g_peak = 0;
                    break
                end
                gindx = find(LFP_mean==g_peak);
            end
            if g_peak==0
                g_peak_freq = 0;
            else
                g_peak_freq = f(LFP_mean==g_peak);
            end 
            g_power = mean(sum(pyy_tot(:,f>19 & f<35),2));
            g_power_std = std(sum(pyy_tot(:,f>19 & f<35),2));%find the integrated power density for gamma 
            % g_low_power = mean(sum(pyy_tot(:,f>35 & f<70),2));
            % g_low_power_std = std(sum(pyy_tot(:,f>35 & f<70),2)); %for low gamma
            % g_high_power = mean(sum(pyy_tot(:,f>=70 & f<100),2));
            % g_high_power_std = std(sum(pyy_tot(:,f>=70 & f<100),2)); %for high gamma
        
            b_peak = max(LFP_mean(f>9 & f<=18)); %9-18 Hz
            bindx = find(LFP_mean==b_peak);
            while b_peak<LFP_mean(bindx+1) %make sure it's not shoulder of gamma peak
                b_peak = max(LFP_mean(5:bindx-1));
                if isempty(b_peak) %this means there was no distinct beta peak
                    b_peak = 0;
                    break
                end
                bindx = find(LFP_mean==b_peak);
            end
            while b_peak<LFP_mean(bindx-1) %make sure it's not shoulder of theta peak
                b_peak = max(LFP_mean(bindx+1:8));
                if isempty(b_peak) %this means there was no distinct beta peak
                    b_peak = 0;
                    break
                end
                bindx = find(LFP_mean==b_peak);
            end
            if b_peak==0
                b_peak_freq = 0;
            else
                b_peak_freq = f(LFP_mean==b_peak);
            end 
            b_power = mean(sum(pyy_tot(:,f>9 & f<18),2));
            b_power_std = std(sum(pyy_tot(:,f>9 & f<18),2));%find the integrated power density for beta
        
            th_peak = max(LFP_mean(f>2 & f<=8)); %2 - 8 Hz
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
            th_power = mean(sum(pyy_tot(:,f>2 & f<=8),2));
            th_power_std = std(sum(pyy_tot(:,f>2 & f<=8),2));%find the integrated power density for theta
            
            %repeat for oscillations between 10-80 Hz

            sec_peak = max(LFP_mean(f>9 & f<80));
            sec_indx = find(LFP_mean==sec_peak);
            while sec_peak<LFP_mean(sec_indx-1) %make sure it's not the shoulder of theta peak
                sec_peak = max(LFP_mean(sec_indx+1:33));
                sec_indx = find(LFP_mean==sec_peak);
            end
            sec_peak_freq = f(LFP_mean==sec_peak);
            % sec_power = mean(sum(pyy_tot(:,f>35 & f<100),2));
            % sec_power_std = std(sum(pyy_tot(:,f>35 & f<100),2));%find the integrated power density for gamma 

        %     max_peak = max(LFP_mean);
        %     max_peak_freq = f(LFP_mean==max_peak);
        
            % errorbar(f(1:freqlength),LFP_mean(1:freqlength),LFP_SEM(1:freqlength))
            % hold on
            % errorbar(f(1:freqlength),LFP_mean_2(1:freqlength),LFP_SEM_2(1:freqlength))
            % xlim([0,100])
            
            figname = sprintf('LFP50_GPFI_mGABAlvl%g_%gglom_inputlvl%g_tilvl%g.png',...
                mGABAindx,gloms,inptindx,tiindx);
        
            errorbar(f(1:freqlength),LFP_mean(1:freqlength),LFP_SEM(1:freqlength))
            xlim([0,100])
%             ylim([0,max(LFP_mean(1:freqlength))])
            ylim([0,1e-3])
            title(append('Spont GPFI, mGABA = ',num2str(mGABA_lvls(mGABAindx)),...
                ' nS, g_{Tonic} = ',num2str(ti_lvls(tiindx)),...
                ' nS, min input = ',num2str(inptlvls(inptindx))));
            
            peak_vals = {['Theta: ' num2str(th_peak) ' at ' num2str(th_peak_freq) ' Hz'], ...
                ['Beta: ' num2str(b_peak) ' at ' num2str(b_peak_freq) ' Hz'], ...
                ['Gamma: ' num2str(g_peak) ' at ' num2str(g_peak_freq) ' Hz']}; %, ...
        %         ['Max: ' num2str(max_peak) ' at ' num2str(max_peak_freq) ' Hz']};
            text(36,0.85*1e-3,peak_vals, 'FontSize',12.5);
        
        %     saveas(gcf,append(save_loc,figname))
            saveas(gcf,figname)
        
            g_peaks(inptcounter,ticounter,mGABAcounter) = g_peak; 
            g_freqs(inptcounter,ticounter,mGABAcounter) = g_peak_freq;
            g_powers(inptcounter,ticounter,mGABAcounter) = g_power;
            g_low_powers(inptcounter,ticounter,mGABAcounter) = g_low_power;
            g_high_powers(inptcounter,ticounter,mGABAcounter) = g_high_power;
            b_peaks(inptcounter,ticounter,mGABAcounter) = b_peak; 
            b_freqs(inptcounter,ticounter,mGABAcounter) = b_peak_freq; 
            b_powers(inptcounter,ticounter,mGABAcounter) = b_power;
            th_peaks(inptcounter,ticounter,mGABAcounter) = th_peak; 
            th_freqs(inptcounter,ticounter,mGABAcounter) = th_peak_freq;
            th_powers(inptcounter,ticounter,mGABAcounter) = th_power;
        
            secondary_peaks(inptcounter,ticounter,mGABAcounter) = sec_peak; 
            secondary_freqs(inptcounter,ticounter,mGABAcounter) = sec_peak_freq;

            g_peaks_sem(inptcounter,ticounter,mGABAcounter) = LFP_SEM(gindx); 
            b_peaks_sem(inptcounter,ticounter,mGABAcounter) = LFP_SEM(bindx);
            th_peaks_sem(inptcounter,ticounter,mGABAcounter) = LFP_SEM(thindx);
            secondary_peaks_sem(inptcounter,ticounter,mGABAcounter) = LFP_SEM(sec_indx);

            g_powers_sem(inptcounter,ticounter,mGABAcounter) = g_power_std/sqrt(numTrials);
            g_low_powers_sem(inptcounter,ticounter,mGABAcounter) = g_low_power_std/sqrt(numTrials);
            g_high_powers_sem(inptcounter,ticounter,mGABAcounter) = g_high_power_std/sqrt(numTrials);
            b_powers_sem(inptcounter,ticounter,mGABAcounter) = b_power_std/sqrt(numTrials);
            th_powers_sem(inptcounter,ticounter,mGABAcounter) = th_power_std/sqrt(numTrials);

            LFP_means(inptcounter,ticounter,mGABAcounter,:) = LFP_mean;
            LFP_SEMs(inptcounter,ticounter,mGABAcounter,:) = LFP_SEM;
        end
    end
end
save('GPFI_peaksfull.mat','g_peaks','g_freqs','g_powers',"g_low_powers","g_high_powers",...
    'b_peaks','b_freqs',"b_powers",'th_peaks','th_freqs',"th_powers",...
    'g_peaks_sem',"g_powers_sem","g_low_powers_sem","g_high_powers_sem",...
    "b_powers_sem",'b_peaks_sem','th_peaks_sem',"th_powers_sem",...
    "secondary_peaks","secondary_peaks_sem","secondary_freqs","LFP_means",...
    "LFP_SEMs","f","freqlength",'mGABA_lvls',...
    'ti_lvls','inptlvls','mGABAindxs','inptindxs','-v7.3')
%% 
load('GPFI_peaksfull.mat')

g_freqs(g_peaks<5e-5) = 0; %set frequency of insignificant peaks to zero
b_freqs(b_peaks<5e-5) = 0;

% g_larger = find(g_peaks>b_peaks);
% b_larger = find(b_peaks>g_peaks);
% secondary_peaks = zeros(size(g_peaks));
% secondary_freqs = zeros(size(g_peaks));
% secondary_peaks(g_larger) = g_peaks(g_larger);
% secondary_freqs(g_larger) = g_freqs(g_larger);
% secondary_peaks(b_larger) = b_peaks(b_larger);
% secondary_freqs(b_larger) = b_freqs(b_larger);

figure('Position',[0, -50, 1120, 840])
tiled = tiledlayout(4,5,"TileSpacing","compact");
tiindx = 5;
tiindx2 = 3;
for inptcounter = fliplr(2:5)
    for mGABAcounter = 3:7
        nexttile
        errorbar(f(1:freqlength),squeeze(LFP_means(inptcounter,tiindx2,mGABAcounter,1:freqlength)),...
            squeeze(LFP_SEMs(inptcounter,tiindx2,mGABAcounter,1:freqlength)))
        xlim([0,80])
        % ylim([0,max(LFP_means(:,tiindx,:,1:freqlength),[],'all')])
        ylim([0,8e-3])
        if mGABAcounter ==3
            ylabel(sprintf('%g', inptlvls(inptindxs(inptcounter))))
        end
        if inptcounter ==2
            xlabel(['\rm ' sprintf('%g nS',mGABA_lvls(mGABAindxs(mGABAcounter)))])
        end 
        if ~(inptcounter==5 && mGABAcounter==7)
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
        end
    end
end
title(tiled,sprintf('Power Spectra (mV^2/Hz vs. Hz), g_{Tonic} = %g nS',ti_lvls(tiindx)),'fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_spectra.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')


% figure('Position',[0, 0,1200,600])
% ticounter = 0;
% for tiindx = ti_indxs
%     ticounter = ticounter+1;
%     subplot(1,3,ticounter)
%     yvalues = inptlvls(inptindxs);
%     xvalues = mGABA_lvls(mGABAindxs);
%     h = heatmap(xvalues,yvalues,squeeze(th_peaks(:,ticounter,:)));
%     h.Title=sprintf('GPFI Peak Theta Power, g_{Tonic} = %g nS',ti_lvls(tiindx));
%     h.YLabel='Min Mean Input';
%     h.XLabel='GABA Conductance, nS';
%     h.YDisplayData = flip(yvalues);
%     h.ColorLimits = [min(th_peaks,[],'all')*1.1,...
%         max(th_peaks,[],"all")*1.1];
%     colormap parula
% end
% figname = sprintf('GPFI_th_peaks.png');
% saveas(gcf,figname)
% 
% 
% figure('Position',[0, 0,1200,600])
% ticounter = 0;
% for tiindx = ti_indxs
%     ticounter = ticounter+1;
%     subplot(1,3,ticounter)
%     yvalues = inptlvls(inptindxs);
%     xvalues = mGABA_lvls(mGABAindxs);
%     h = heatmap(xvalues,yvalues,squeeze(th_powers(:,ticounter,:)));
%     h.Title=sprintf('GPFI Integrated Theta PSD, g_{Tonic} = %g nS',ti_lvls(tiindx));
%     h.YLabel='Min Mean Input';
%     h.XLabel='GABA Conductance, nS';
%     h.YDisplayData = flip(yvalues);
%     h.ColorLimits = [min(th_powers,[],'all')*1.1,...
%         max(th_powers,[],"all")*1.1];
%     colormap parula
% end
% figname = sprintf('GPFI_th_powers.png');
% saveas(gcf,figname)
% 
% 
figure('Position',[0, 0,1120,420])
ticounter = 0;
tiled = tiledlayout(1,3,"TileSpacing","compact");
ti_indxs = [1,2,5];
for tiindx = ti_indxs
    ticounter = ticounter+1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(g_peaks(:,ticounter,:)));
    h.Title=['\rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(g_peaks,[],'all'),...
        max(g_peaks,[],"all")];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap turbo
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Peak Power (19-35 Hz), Graded and Firing Inhibition','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_highpeaks.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')


figure('Position',[0, 0,1120,420])
ticounter = 0;
tiled = tiledlayout(1,3,"TileSpacing","compact");
ti_indxs = [1,2,5];
for tiindx = ti_indxs
    ticounter = ticounter+1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(b_peaks(:,ticounter,:)));
    h.Title=['\rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(b_peaks,[],'all'),...
        max(b_peaks,[],"all")];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap turbo
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Peak Power (9-18 Hz), Graded and Firing Inhibition','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_midpeaks.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')


figure('Position',[0, 0,1120,420])
ticounter = 0;
tiled = tiledlayout(1,3,"TileSpacing","compact");
ti_indxs = [1,2,5];
for tiindx = ti_indxs
    ticounter = ticounter+1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(th_peaks(:,ticounter,:)));
    h.Title=['\rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(th_peaks,[],'all'),...
        max(th_peaks,[],"all")];
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    colormap turbo
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Peak Power (2-8 Hz), Graded and Firing Inhibition','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_lowpeaks.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')




figure('Position',[0, 0,1120,420])
ticounter = 0;
tiled = tiledlayout(1,3,"TileSpacing","compact");
% secondary_freqs([1,2],3,1) = nan;
for tiindx = ti_indxs
    ticounter = ticounter+1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(g_freqs(:,ticounter,:)),...
        Colormap=turbo,MissingDataColor='r');
    h.Title=['\rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    h.ColorLimits = [0,35];
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Dominant Frequency (19-35 Hz), Graded and Firing Inhibition','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_highpeakfreqs.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')

figure('Position',[0, 0,1120,420])
ticounter = 0;
tiled = tiledlayout(1,3,"TileSpacing","compact");
for tiindx = ti_indxs
    ticounter = ticounter+1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(b_freqs(:,ticounter,:)),...
        Colormap=turbo,MissingDataColor='r');
    h.Title=['\rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    h.ColorLimits = [0,35];
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Dominant Frequency (9-18 Hz), Graded and Firing Inhibition','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_midpeakfreqs.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')


figure('Position',[0, 0,1120,420])
ticounter = 0;
tiled = tiledlayout(1,3,"TileSpacing","compact");
for tiindx = ti_indxs
    ticounter = ticounter+1;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(th_freqs(:,ticounter,:)),...
        Colormap=turbo,MissingDataColor='r');
    h.Title=['\rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))];
    h.YDisplayData = flip(yvalues);
    if tiindx ~=1
        h.YDisplayLabels = nan(size(h.YDisplayData));
    end
    h.ColorbarVisible = 'off';
    h.CellLabelColor = 'none';
    h.ColorLimits = [0,8];
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Dominant Frequency (2-8 Hz), Graded and Firing Inhibition','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_lowpeakfreqs.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')

%% Plot the two frequency bands on the same figure, frequency then power

figure('Position',[0, 0,1120,420])
tiled = tiledlayout(1,2);
% secondary_freqs([1,2],3,1) = nan;
for tiindx = 5
    ticounter = 3;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(b_freqs(:,ticounter,:)),...
        Colormap=turbo,MissingDataColor='r');
    h.Title='\rm Dominant Frequency (9-18 Hz)';
    h.YDisplayData = flip(yvalues);
    h.CellLabelColor = 'none';
    h.ColorLimits = [0,35];
end
h.ColorbarVisible = 'on'; %only on for the last map

for tiindx = 5
    ticounter = 3;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(g_freqs(:,ticounter,:)),...
        Colormap=turbo,MissingDataColor='r');
    h.Title = "\rm Dominant Frequency (19-35 Hz)";
    h.YDisplayData = flip(yvalues);
    h.CellLabelColor = 'none';
    h.ColorLimits = [0,35];
end
h.ColorbarVisible = 'on'; %only on for the last map
% title(tiled,['Graded and Firing Inhibition, \rm ' sprintf('g_{Tonic} = %g nS',ti_lvls(tiindx))],'fontweight','bold')
title(tiled,'Graded and Firing Inhibition, g_{Tonic} = 1.8 nS','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_twobands_freq.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')



figure('Position',[0, 0,1120,420])
tiled = tiledlayout(1,2);
for tiindx = 5
    ticounter = 3;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(b_peaks(:,ticounter,:)));
    h.Title='\rm Peak Power (9-18 Hz)';
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(b_peaks,[],'all'),...
        max(b_peaks,[],"all")];
    h.CellLabelColor = 'none';
    colormap turbo
end
h.ColorbarVisible = 'on'; %only on for the last map

for tiindx = 5
    ticounter = 3;
    nexttile
    yvalues = inptlvls(inptindxs);
    xvalues = mGABA_lvls(mGABAindxs);
    h = heatmap(xvalues,yvalues,squeeze(g_peaks(:,ticounter,:)));
    h.Title='\rm Peak Power (19-35 Hz)';
    h.YDisplayData = flip(yvalues);
    h.ColorLimits = [min(g_peaks,[],'all'),...
        max(g_peaks,[],"all")];
    h.CellLabelColor = 'none';
    colormap turbo
end
h.ColorbarVisible = 'on'; %only on for the last map
title(tiled,'Graded and Firing Inhibition, g_{Tonic} = 1.8 nS','fontweight','bold')
xlabel(tiled,'MC GABA Conductance, nS')
ylabel(tiled,'Minimum \mu_r')
fontsize(18,'points')
figname = sprintf('Fig_GPFI_twobands_peaks.png');
fig = gcf;
exportgraphics(fig,figname,'Resolution','600')




%% Plot vs TI as well (just the peaks)
% for inptindx = [7,9,11]
%     figure
%     hold on
%     errorbar(ti_lvls,g_peaks(inptindx,:),g_peaks_sem(inptindx,:))
%     errorbar(ti_lvls,b_peaks(inptindx,:),b_peaks_sem(inptindx,:))
%     errorbar(ti_lvls,th_peaks(inptindx,:),th_peaks_sem(inptindx,:))
%     legend('gamma', 'beta', 'theta')
%     title(sprintf('Spont, Peak Power vs g_{Tonic}, Min Mean Input = %g',...
%         inptlvls(inptindx)))
%     xlabel('g_{Tonic}, nS')
%     ylabel('Peak Power, V^2/Hz')
%     ylim([0,1.5e-4])
%     figname = sprintf('VDSI1_Peaks_inputlvl%g.png',inptindx);
%     saveas(gcf,figname)
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
