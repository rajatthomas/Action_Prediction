function freq_distrib = get_freq(data)

% Set the bounds of the frequency here.
low_freq = 4;
high_freq = 160;
freq_width = 2;

cfg = [];
cfg.keeptrials       = 'yes'; % keep trials for statistics
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq   = 30;    % smooth ERP with low-pass filter
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 1;     % reduce slow drifts
cfg.preproc.detrend  = 'yes';

ERP = ft_timelockanalysis(cfg, data);

% HGP

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

cfg.trials = find(data.trialinfo <= 20 & data.trialinfo ~= -100); % select 'intact' trials
TFR_intact = ft_freqanalysis(cfg,data);

cfg.trials = find(data.trialinfo > 20 & data.trialinfo ~= 100); % select 'scrambled' trials
TFR_scrambled   = ft_freqanalysis(cfg,data);

% create HGP as empty timelock structure with same dimensions as ERP, values will be filled in in the next steps
HGP_intact = rmfield(ERP,{'trial','avg','var'});
HGP_scrambled   = rmfield(ERP_scrambled,{'trial','avg','var'});

% correct for the 1/f dropoff
freqcorr = reshape(TFR_intact.freq.^2,[1 1 length(TFR_intact.freq)]); %this vector accounts for the 1/f dropoff
% use repmat to create a matrix the same size as the TFR data
freqcorr_intact = repmat(freqcorr,[size(TFR_intact.powspctrm,1) size(TFR_intact.powspctrm,2) 1 length(TFR_intact.time)]);
freqcorr_scrambled   = repmat(freqcorr,[size(TFR_scrambled.powspctrm,1) size(TFR_scrambled.powspctrm,2) 1 length(TFR_scrambled.time)]);
 
% multiply data with freqcorr matrix and average over frequencies
HGP_intact.trial = squeeze(nanmean(TFR_intact.powspctrm(:,:,:,:) .* freqcorr_intact,3));
HGP_scrambled.trial   = squeeze(nanmean(TFR_scrambled.powspctrm(:,:,:,:) .* freqcorr_scrambled,3));
 
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


end