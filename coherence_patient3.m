%% analysis
%
close all; clear all
addpath('/home/rajat/Dropbox/Teresa_Rajat/AP_EEG/fieldtrip-20170618/')

prefix = 'P3_run_';
seg_type = 'movie'; %'camera'; % choose btn ['movie', 'camera']
config = 'nofilter_nobaseline';

%--- First viewing

runID = 'A1';
load([prefix runID '_seg_' seg_type '_' config '.mat']);
dataA1 = data;

runID = 'B1';
load([prefix runID '_seg_' seg_type  '_' config '.mat']);
dataB1 = data;

runID = 'C1';
load([prefix runID '_seg_' seg_type  '_' config '.mat']);
dataC1 = data;

%--- second viewing

runID = 'A2';
load([prefix runID '_seg_' seg_type  '_' config '.mat']);
dataA2 = data;

runID = 'B2';
load([prefix runID '_seg_' seg_type  '_' config '.mat']);
dataB2 = data;

runID = 'C2';
load([prefix runID '_seg_' seg_type  '_' config '.mat']);
dataC2 = data;


data = ft_appenddata([], dataA1, dataB1, dataC1, dataA2, dataB2, dataC2);

%------------- Channels of interests -------------------------------------%
Parietal = [42,33,41];
PM = [30,31];
Visual = [45,38,55,59,37,65,60,67,36,39];
Temporal = 12;

selectedChannels = [Parietal PM Visual Temporal];
nSelectedChannels = length(selectedChannels);



%------------- Channels of interests -------------------------------------%


channelcmb = cell(nSelectedChannels*nSelectedChannels, 2);
for i=1:nSelectedChannels
    for j = 1:nSelectedChannels
        channelcmb{(i-1)*nSelectedChannels + j, 1} = data.label{i};
        channelcmb{(i-1)*nSelectedChannels + j, 2} = data.label{j};
    end
end



%%
ntrials = length(data.trialinfo);

%--------- Bands of interest -------------------------------------------
alpha = 8:2:12;
low_beta = 12:2:20;
high_beta = 20:2:28;
low_gamma = [28:2:46 54:2:64];
high_gamma = [64:2:98 102:2:148 152:2:196];

% freqs of interest
foi = unique([alpha low_beta high_beta low_gamma high_gamma]);

foi = 5:5:180;


nChannelcmb = length(channelcmb);
n_rpt = 12;
nfreq = length(foi);

nchannels = length(data.label);
coh_alltrials = zeros(ntrials, nchannels, nchannels, nfreq);

%%
data_intact = data;
data_intact.trial = data.trial(data.trialinfo <= 20);
data_intact.time  = data.time(data.trialinfo <= 20);
data_intact.sampleinfo  = data.sampleinfo(data.trialinfo <= 20,:);
data_intact.trialinfo = data.trialinfo(data.trialinfo <= 20);

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foi        = foi;
cfg.toi        = 0:0.2:20; %-0.3:0.1:0.8;
%cfg.pad        = 'nextpow2'; % efficient FFT
cfg.t_ftimwin  = 7./cfg.foi;
cfg.tapsmofrq = 4;
%cfg.trials     = trial_i;
%cfg.channelcmb = channelcmb; %{'A49-E-A48-E' 'A52-E-A51-E'; 'A49-E-A48-E' 'A52-E-A51-E';}; % {'A60-E-A59-E' 'A3-E-A2-E';};

freq          = ft_freqanalysis(cfg, data_intact);
cfg           = [];
cfg.method    = 'coh';
coh_intact           = ft_connectivityanalysis(cfg, freq);

%%
data_scrambled = data;
data_scrambled.trial = data.trial(data.trialinfo > 20);
data_scrambled.time  = data.time(data.trialinfo > 20);
data_scrambled.sampleinfo  = data.sampleinfo(data.trialinfo <= 20,:);
data_scrambled.trialinfo = data.trialinfo(data.trialinfo > 20);

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foi        = foi;
cfg.toi        = 0:0.2:20; %-0.3:0.1:0.8;
%cfg.pad        = 'nextpow2'; % efficient FFT
cfg.t_ftimwin  = 7./cfg.foi;
cfg.tapsmofrq = 4;
%cfg.trials     = trial_i;
%cfg.channelcmb = channelcmb; %{'A49-E-A48-E' 'A52-E-A51-E'; 'A49-E-A48-E' 'A52-E-A51-E';}; % {'A60-E-A59-E' 'A3-E-A2-E';};

freq          = ft_freqanalysis(cfg, data_scrambled);
cfg           = [];
cfg.method    = 'coh';
coh_scrambled           = ft_connectivityanalysis(cfg, freq);

%%
for trial_i = 1:ntrials
    
    cfg           = [];
    cfg.method    = 'mtmfft';
    cfg.taper     = 'dpss';
    cfg.output    = 'fourier';
    cfg.foi        = foi;
    cfg.toi        = 0:0.1:20; %-0.3:0.1:0.8;
    %cfg.pad        = 'nextpow2'; % efficient FFT
    cfg.t_ftimwin  = 4./cfg.foi;
    cfg.tapsmofrq = 1;
    %cfg.trials     = trial_i; 
    cfg.channelcmb = channelcmb; %{'A49-E-A48-E' 'A52-E-A51-E'; 'A49-E-A48-E' 'A52-E-A51-E';}; % {'A60-E-A59-E' 'A3-E-A2-E';};
    
    freq          = ft_freqanalysis(cfg, data);
    cfg           = [];
    cfg.method    = 'coh';
    coh           = ft_connectivityanalysis(cfg, freq);
    
    coh_alltrials(trial_i, :, :, :) = coh.cohspctrm;
    
end

frequencies = coh.freq;

data.trialinfo(end-7) = 12;
intact_trials = (data.trialinfo <= 20 & data.trialinfo ~= -100);
coh_intact = coh_alltrials(intact_trials, :, :, :);

scram_trials = (data.trialinfo > 20 & data.trialinfo ~= 100);
coh_scram = coh_alltrials(scram_trials, :, :, :);


index_alpha      = coh.freq > (alpha(1)-1) & (coh.freq < (alpha(end)+1));
index_low_beta   = coh.freq > (low_beta(1)-1) & (coh.freq < (low_beta(end)+1));
index_high_beta  = coh.freq > (high_beta(1)-1) & (coh.freq < (high_beta(end)+1));
index_low_gamma  = coh.freq > (low_gamma(1)-1) & (coh.freq < (low_gamma(end)+1));
index_high_gamma = coh.freq > (high_gamma(1)-1) & (coh.freq < (high_gamma(end)+1));


coh_intact_alpha      = mean(coh_intact(:, :, :, index_alpha), 4); % Average across frequencies
coh_intact_low_beta   = mean(coh_intact(:, :, :, index_low_beta), 4); 
coh_intact_high_beta  = mean(coh_intact(:, :, :, index_high_beta), 4); 
coh_intact_low_gamma  = mean(coh_intact(:, :, :, index_low_gamma), 4); 
coh_intact_high_gamma = mean(coh_intact(:, :, :, index_high_gamma), 4); 


coh_scram_alpha      = mean(coh_scram(:, :, :, index_alpha), 4); % Average across frequencies
coh_scram_low_beta   = mean(coh_scram(:, :, :, index_low_beta), 4); 
coh_scram_high_beta  = mean(coh_scram(:, :, :, index_high_beta), 4); 
coh_scram_low_gamma  = mean(coh_scram(:, :, :, index_low_gamma), 4); 
coh_scram_high_gamma = mean(coh_scram(:, :, :, index_high_gamma), 4); 

save('coh_measures.mat', 'coh_intact_alpha', 'coh_intact_low_beta', 'coh_intact_high_beta', 'coh_intact_low_gamma', 'coh_intact_high_gamma', ...
    'coh_scram_alpha', 'coh_scram_low_beta', 'coh_scram_high_beta', 'coh_scram_low_gamma', 'coh_scram_high_gamma','frequencies');

%% T-test


intact_trial_list = data.trialinfo(intact_trials);
scram_trial_list = data.trialinfo(scram_trials);

[~, sorted_intact] = sort(intact_trial_list);
[~, sorted_scram] = sort(scram_trial_list);

intact_bands = {'coh_intact_alpha', 'coh_intact_low_beta', 'coh_intact_high_beta', 'coh_intact_low_gamma', 'coh_intact_high_gamma'};
scram_bands = {'coh_scram_alpha', 'coh_scram_low_beta', 'coh_scram_high_beta', 'coh_scram_low_gamma', 'coh_scram_high_gamma'};
band_name = {'coh_alpha', 'coh_low_beta', 'coh_high_beta', 'coh_low_gamma', 'coh_high_gamma'};


for i=1:length(band_name)

    figure
    coh_band_intact = eval(intact_bands{i});
    coh_band_scram = eval(scram_bands{i});
    [h, p, ci, stats] = ttest(coh_band_intact(sorted_intact, :, :), coh_band_scram(sorted_scram, :, :));
    T_thresh = 3.0;
    imagesc(squeeze(stats.tstat .* (abs(stats.tstat) > T_thresh)))
    caxis([-5 5])

    % Sample the channels that show high coherence
    [chx, chy] = find(squeeze(abs(stats.tstat) > T_thresh));
    save([band_name{i} '_channels.mat'], 'chx', 'chy');
end





%% 

avg_coh_intact = squeeze(mean(coh_intact,1));
avg_coh_intact = reshape(avg_coh_intact, [size(avg_coh_intact,1)*size(avg_coh_intact,2) size(avg_coh_intact,3)]);
plot(coh.freq, avg_coh_intact');

%plot(coh.freq, mean(avg_coh_intact,1));


















%%
% % %     cfg            = [];
% % %     cfg.output     = 'powandcsd';
% % %     cfg.method     = 'mtmconvol'; %'mtmfft';
% % % 
% % % 
% % % %fries_low = 7:1:17;
% % % %fries_high = 42:4:93;
% % % 
% % % 
% % %     cfg.foi        = unique([ alpha low_beta high_beta low_gamma high_gamma]);
% % %     cfg.toi        = 0:0.1:20; %-0.3:0.1:0.8;
% % %     cfg.t_ftimwin  = 4./cfg.foi;
% % %     %cfg.tapsmofrq  = 1;
% % %     cfg.taper = 'hanning';
% % %     %cfg.keeptrials = 'yes';
% % % 
% % %     cfg.trials     = 12; %[1 2]; % trial_i; 
% % %     cfg.channelcmb = channelcmb; %{'A49-E-A48-E' 'A52-E-A51-E'; 'A49-E-A48-E' 'A52-E-A51-E';}; % {'A60-E-A59-E' 'A3-E-A2-E';};
% % %     
% % %     if (data.trialinfo(trial_i) <= 20 & data.trialinfo(trial_i) ~= -100)
% % %         freq_intact = ft_freqanalysis(cfg, data);
% % %     end
% % %     
% % %     if(data.trialinfo(trial_i) > 20 & data.trialinfo(trial_i) ~= 100)
% % %         freq_scrambled = ft_freqanalysis(cfg, data);
% % %     end
% % %     
% % %     %cfg.channel    = 'A*'; %'all'; %{'A60-E-A59-E' 'A3-E-A2-E';};   %'all';
% % %     
% % %     %cfg.keeptrials = 'yes';
% % %     
% % %     
% % % 
% % % %     cfg.trials     = (data.trialinfo == -100 | data.trialinfo == 100);
% % % %     freq_baseline = ft_freqanalysis(cfg, data);
% % % 
% % % 
% % %     cfg            = [];
% % %     cfg.method     = 'coh';
% % %     %cfg.channel = 'A*';
% % %     %cfg.channelcmb = {'A49-E-A48-E' 'A52-E-A51-E'};
% % % 
% % %     
% % %      
% % % 
% % %     if (data.trialinfo(trial_i) <= 20 & data.trialinfo(trial_i) ~= -100)
% % %         fd_intact             = ft_connectivityanalysis(cfg, freq_intact);
% % %         coh_intact(intact_trial,:,:) = nanmean(fd_intact.cohspctrm,3);
% % %         
% % %         
% % %         %index_alpha = fd_intact.freq > (alpha(1)-1) & (fd_intact.freq < (alpha(end)+1));
% % %         index_low_beta = fd_intact.freq > (low_beta(1)-1) & (fd_intact.freq < (low_beta(end)+1));
% % %         index_high_beta = fd_intact.freq > (high_beta(1)-1) & (fd_intact.freq < (high_beta(end)+1));
% % %         index_low_gamma = fd_intact.freq > (low_gamma(1)-1) & (fd_intact.freq < (low_gamma(end)+1));
% % %         index_high_gamma = fd_intact.freq > (high_gamma(1)-1) & (fd_intact.freq < (high_gamma(end)+1));
% % %     
% % %         %coh_intact_alpha(intact_trial,:) = mean(coh_intact(intact_trial,:,index_alpha),3);
% % %         coh_intact_low_beta(intact_trial,:) = mean(coh_intact(intact_trial,:,index_low_beta),3);
% % %         coh_intact_high_beta(intact_trial,:) = mean(coh_intact(intact_trial,:,index_high_beta),3);
% % %         coh_intact_low_gamma(intact_trial,:) = mean(coh_intact(intact_trial,:,index_low_gamma),3);
% % %         coh_intact_high_gamma(intact_trial,:) = mean(coh_intact(intact_trial,:,index_high_gamma),3);
% % %         intact_trial = intact_trial + 1;
% % %     end
% % % 
% % %     if(data.trialinfo(trial_i) > 20 & data.trialinfo(trial_i) ~= 100)
% % %         fd_scrambled          = ft_connectivityanalysis(cfg, freq_scrambled);
% % %         coh_scram(scram_trial,:,:) = nanmean(fd_scrambled.cohspctrm,3);
% % %         
% % %         %index_alpha = fd_scrambled.freq > (alpha(1)-1) & (fd_scrambled.freq < (alpha(end)+1));
% % %         index_low_beta = fd_scrambled.freq > (low_beta(1)-1) & (fd_scrambled.freq < (low_beta(end)+1));
% % %         index_high_beta = fd_scrambled.freq > (high_beta(1)-1) & (fd_scrambled.freq < (high_beta(end)+1));
% % %         index_low_gamma = fd_scrambled.freq > (low_gamma(1)-1) & (fd_scrambled.freq < (low_gamma(end)+1));
% % %         index_high_gamma = fd_scrambled.freq > (high_gamma(1)-1) & (fd_scrambled.freq < (high_gamma(end)+1));
% % %         
% % %         %coh_scram_alpha(scram_trial,:) = mean(coh_scram(scram_trial,:,index_alpha),3);
% % %         coh_scram_low_beta(scram_trial,:) = mean(coh_scram(scram_trial,:,index_low_beta),3);
% % %         coh_scram_high_beta(scram_trial,:) = mean(coh_scram(scram_trial,:,index_high_beta),3);
% % %         coh_scram_low_gamma(scram_trial,:) = mean(coh_scram(scram_trial,:,index_low_gamma),3);
% % %         coh_scram_high_gamma(scram_trial,:) = mean(coh_scram(scram_trial,:,index_high_gamma),3);
% % %         scram_trial = scram_trial + 1;
% % %     end
% % % %fd_baseline          = ft_connectivityanalysis(cfg, freq_baseline);
% % % %fdfourier      = ft_connectivityanalysis(cfg, freqfourier);
% % % %-----------------
% % % % % cfg                  = [];
% % % % % cfg.parameter        = 'cohspctrm';
% % % % % %cfg.xlim             = [5 80];
% % % % % cfg.refchannel       = 'A49-E-A48-E';
% % % % % cfg.showlabels       = 'yes';
% % % % % %cfg.channel = 'A*';
% % % % % %figure; ft_singleplotER(cfg, fd_intact, fd_scrambled);
% % % % % 
% % % % % coh_intact = squeeze(fd_intact.cohspctrm);
% % % % % coh_scram = squeeze(fd_scrambled.cohspctrm);
% % % % % coh_baseline = squeeze(fd_baseline.cohspctrm);
% % % 
% % % % for i = 1:size(coh_intact,2)
% % % %     coh_intact(:,i) = smooth(coh_intact(:,i),5 );
% % % %     coh_scram(:,i) = smooth(coh_scram(:,i), 5);
% % % %     coh_baseline(:,i) = smooth(coh_baseline(:,i), 5);
% % % % end
% % % %
% % % 
% % %     

    

    
end



%%




%%
ci = nanmean(coh_intact,3);
cs = nanmean(coh_scram,3);
cb = nanmean(coh_baseline,3);



for i=1:size(cb,1)
 if strcmp(fd_intact.labelcmb(i,1), 'A49-E-A48-E') || strcmp(fd_intact.labelcmb(i,2), 'A49-E-A48-E')
    %fd_intact.labelcmb(i,1)
    %'fd_intact.labelcmb(i,2)
    plot(fd_intact.freq, ci(i,:))
    hold on
    plot(fd_scrambled.freq, cs(i,:))
    plot(fd_baseline.freq, cb(i,:))
    
    title([fd_intact.labelcmb(i,1) '-' fd_intact.labelcmb(i,2)])
    pause
    close all
 end

end

%%
ci = nanmean(coh_intact,3);
cs = nanmean(coh_scram,3);
cb = nanmean(coh_baseline,3);


for i=1:size(cb,1)
 %if strcmp(fd_intact.labelcmb(i,1), 'A49-E-A48-E') || strcmp(fd_intact.labelcmb(i,2), 'A49-E-A48-E')
    %fd_intact.labelcmb(i,1)
    %'fd_intact.labelcmb(i,2)
    %plot(fd_intact.freq, ci(i,:))
    hold on
    plot(fd_scrambled.freq, ci(i,:)-cs(i,:))
    %plot(fd_baseline.freq, cb(i,:))
    
    %title([fd_intact.labelcmb(i,1) '-' fd_intact.labelcmb(i,2)])
    %pause
    %close all
 %end

end

%%
times=-0.3:0.1:0.8;
for t_index=1:length(times)

t = times(t_index);

plot(fd_intact.freq, coh_intact(:,t_index))
hold on
plot(fd_scrambled.freq, coh_scram(:,t_index))
title(['Coherence spectrum at time ' num2str(t)])
pause
close all
end