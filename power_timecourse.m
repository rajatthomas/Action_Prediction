% Granger analysis

%% load the data 
close all; clear all
addpath('/home/rajat/Dropbox/Teresa_Rajat/AP_EEG/fieldtrip-master/')

runID = 'B';
seg_type = 'camera'; %'camera'; % choose btn ['movie', 'camera']
load(['run_' runID '_seg_' seg_type '.mat']);


%% 

disp('Performing ERP on intact and scrambled')

cfg = [];
cfg.keeptrials       = 'yes'; % keep trials for statistics
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq   = 30;    % smooth ERP with low-pass filter
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 1;     % reduce slow drifts
cfg.preproc.detrend  = 'yes';


cfg.trials = find(data.trialinfo <= 20 & data.trialinfo ~= -100); % select only 'object' trials (all M movies have label <= 20)
ERP_intact = ft_timelockanalysis(cfg, data);
 
cfg.trials = find(data.trialinfo > 20 & data.trialinfo ~= 100); % (all S movies are > 20)
ERP_scrambled = ft_timelockanalysis(cfg, data);
 

%% HGP

low_freq = 5;
high_freq = 150;
freq_width = 2;
% time-frequency analysis
cfg            = [];
cfg.method     = 'tfr';
cfg.keeptrials = 'yes';
cfg.toi        = data.time{1}; % keep full temporal resolution
cfg.foi        = low_freq:freq_width:high_freq;

if low_freq < 50
    rm_index = round((50-low_freq)/freq_width) + 1;
    rm_index = [rm_index rm_index+10 rm_index+20];
else
    rm_index = round((100-low_freq)/freq_width) + 1;
    rm_index = [rm_index rm_index+10 rm_index+20];
end

cfg.foi(rm_index) = []; % leave out harmonics of 50 Hz
cfg.width      = 10 * ones(1,length(cfg.foi));

cfg.trials = find(data.trialinfo <= 20 & data.trialinfo ~= -100); % select 'object' trials
TFR_intact = ft_freqanalysis(cfg,data);

cfg.trials = find(data.trialinfo > 20 & data.trialinfo ~= 100); % select 'face' trials
TFR_scrambled   = ft_freqanalysis(cfg,data);

save(['TFR_' runID '_seg_' seg_type '.mat'],'TFR_intact', 'TFR_scrambled', '-v7.3')

%%
% create HGP as empty timelock structure with same dimensions as ERP, values will be filled in in the next steps
HGP_intact = rmfield(ERP_intact,{'trial','avg','var'});
HGP_scrambled   = rmfield(ERP_scrambled,{'trial','avg','var'});

% correct for the 1/f dropoff
freqcorr = reshape(TFR_intact.freq.^2,[1 1 length(TFR_intact.freq)]); %this vector accounts for the 1/f dropoff
% use repmat to create a matrix the same size as the TFR data
freqcorr = repmat(freqcorr,[size(TFR_intact.powspctrm,1) size(TFR_intact.powspctrm,2) 1 length(TFR_intact.time)]);

% multiply data with freqcorr matrix and average over frequencies
%HGP_intact.trial = squeeze(nanmean(TFR_intact.powspctrm(:,:,:,:) .* freqcorr_intact,3));
%HGP_scrambled.trial   = squeeze(nanmean(TFR_scrambled.powspctrm(:,:,:,:) .* freqcorr_scrambled,3));

% multiply data with freqcorr matrix and average over frequencies
%%

cfg            = [];
cfg.method     = 'tfr';
cfg.keeptrials = 'yes';
cfg.toi        = data.time{1}; % keep full temporal resolution
cfg.foi        = low_freq:freq_width:high_freq;

freq_l = 80;
freq_h = 120;

closest = abs(cfg.foi - freq_l);
index_l = find(closest == min(closest));
index_l = index_l(1);

closest = abs(cfg.foi - freq_h);
index_h = find(closest == min(closest));
index_h = index_h(1);


HGP_intact.trial = squeeze(nanmean(TFR_intact.powspctrm(:,:,index_l:index_h,:) .* freqcorr_intact(:,:,index_l:index_h,:),3));
HGP_scrambled.trial   = squeeze(nanmean(TFR_scrambled.powspctrm(:,:,index_l:index_h,:) .* freqcorr_scrambled(:,:,index_l:index_h,:),3));


% calculate mean and variance
HGP_intact.avg = squeeze(nanmean(HGP_intact.trial,1));
HGP_intact.var = squeeze(nanvar(HGP_intact.trial,1));
HGP_scrambled.avg   = squeeze(nanmean(HGP_scrambled.trial,1));
HGP_scrambled.var   = squeeze(nanvar(HGP_scrambled.trial,1));
 
% baseline correction
cfg          = [];
cfg.baseline = [-.3 -.05];
    
HGP_intact_bl = ft_timelockbaseline(cfg,HGP_intact);
HGP_scrambled_bl   = ft_timelockbaseline(cfg,HGP_scrambled);

%%
cfg           = [];
cfg.parameter = 'avg';
cfg.xlim      = [-.3 .7];
cfg.channel   = 38; 
figure, ft_singleplotER(cfg,HGP_intact_bl,HGP_scrambled_bl)
hold on
plot([0, 0],[1e7,10e7],'k')

cfg.channel   = 42; 
figure, ft_singleplotER(cfg,HGP_intact_bl,HGP_scrambled_bl)
hold on
plot([0, 0],[1e7,10e7],'k')

cfg.channel   = 43; 
figure, ft_singleplotER(cfg,HGP_intact_bl,HGP_scrambled_bl)
hold on
plot([0, 0],[1e7,10e7],'k')

cfg.channel   = 46; 
figure, ft_singleplotER(cfg,HGP_intact_bl,HGP_scrambled_bl)
hold on
plot([0, 0],[1e7,10e7],'k')

cfg.channel   = 47; 
figure, ft_singleplotER(cfg,HGP_intact_bl,HGP_scrambled_bl)
hold on
plot([0, 0],[1e7,10e7],'k')

%%

close all

ch=34;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 hold on
 
 ch=38;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 hold on
 
 ch=42;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 
 ch=43;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 
 ch=46;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 
 ch=47;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 
 %%
 close all
 
 ch=38;
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_intact.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 hold on
 
 freq = TFR_intact.freq;
 P=nanmean(squeeze(nanmean(TFR_scrambled.powspctrm(:,ch,:,:),1)),2);
 P = log10(P);
 plot(freq,P)
 hold on
 
 