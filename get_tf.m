function component = get_tf(data, timeOrFreq)

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

% time-frequency analysis
cfg            = [];
cfg.method     = 'tfr';
cfg.keeptrials = 'yes';
cfg.toi        = data.time{1}; % keep full temporal resolution
cfg.foi        = low_freq:freq_width:high_freq;
cfg.width      = 10 * ones(1,length(cfg.foi));
TFR = ft_freqanalysis(cfg,data);

% create POW as empty timelock structure with same dimensions as ERP, values will be filled in in the next steps
component = rmfield(ERP,{'trial','avg','var'});

% correct for the 1/f dropoff
freqcorr = reshape(TFR.freq.^2,[1 1 length(TFR.freq)]); %this vector accounts for the 1/f dropoff
% use repmat to create a matrix the same size as the TFR data
freqcorr = repmat(freqcorr,[size(TFR.powspctrm,1) size(TFR.powspctrm,2) 1 length(TFR.time)]);

if strcmp(timeOrFreq, 'power_timecourse')
    % multiply data with freqcorr matrix and average over frequencies
    component.trial = squeeze(nanmean(TFR.powspctrm(:,:,:,:) .* freqcorr,3));
    component.dimord = 'rpt_chan_time';
end

if strcmp(timeOrFreq, 'power_spectrum')
    component.trial = squeeze(nanmean(TFR.powspctrm(:,:,:,:) .* freqcorr,4));
    component.dimord = 'rpt_chan_freq';
end

% calculate mean and variance
component.avg = squeeze(nanmean(component.trial,1));
component.var = squeeze(nanvar(component.trial,1));
component.freq = cfg.foi;

end