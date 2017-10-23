% analysis

%%
addpath('/home/rajat/Dropbox/Teresa_Rajat/AP_EEG/fieldtrip-master/')

runID = 'A';
seg_type = 'camera'; %'camera'; % choose btn ['movie', 'camera']
load(['run_' runID '_seg_' seg_type '.mat']);


parietalA= [58 59;59 60;60 61;61 62;53 54;54 55;55 56;48 49;41 42];
visualA= [51 52; 45 46;46 47;35 36; 30 32;31 32;25 26;26 27];
parietalB=[58 59;59 60;23 24;24 25;53 43];
visualB=[61 62;55 56;49 50;29 30; 30 31;33 34;34 35;35 36; 36 37;38 39;39 40;40 41; 41 42];

label_parietalA = cell(1, length(parietalA));
for i=1:length(parietalA)
    label_parietalA{i} = sprintf('A%d-E-A%d-E', parietalA(i,1), parietalA(i,2));
end

label_parietalB = cell(1, length(parietalB));
for i=1:length(parietalB)
    label_parietalB{i} = sprintf('B%d-E-B%d-E', parietalB(i,1), parietalB(i,2));
end

label_visualA = cell(1, length(visualA));
for i=1:length(parietalA)
    label_visualA{i} = sprintf('A%d-E-A%d-E', visualA(i,1), visualA(i,2));
end


label_visualB = cell(1, length(visualB));
for i=1:length(parietalB)
    label_visualB{i} = sprintf('B%d-E-B%d-E', visualB(i,1), visualB(i,2));
end


% cfg            = [];
% cfg.output     = 'powandcsd';
% cfg.method     = 'mtmfft';
% cfg.foilim     = [1 200];
% cfg.tapsmofrq  = 5;
% cfg.keeptrials = 'yes';
% cfg.channel    = 'all';
% freq           = ft_freqanalysis(cfg, data);

% cfg            = [];
% cfg.output     = 'fourier';
% cfg.method     = 'mtmfft';
% cfg.foilim     = [5 200];
% cfg.tapsmofrq  = 5;
% cfg.keeptrials = 'yes';
% cfg.channel    = 'all';
% freqfourier    = ft_freqanalysis(cfg, data);

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
 
% baseline correction
%cfg = [];
%cfg.baseline = [-.3 -.05];
 
%ERP_intact_bl = ft_timelockbaseline(cfg,ERP_intact);
%ERP_scrambled_bl = ft_timelockbaseline(cfg,ERP_scrambled);

%% figure tigo test trigger

figure
cfg = [];
cfg.parameter = 'avg';
cfg.channel = 40;

if strcmp(seg_type, 'movie')
    cfg.xlim = [-5 20]; %[-.5 1.0];
    cfg.baseline = [-3 -0.2];
else
    cfg.xlim = [-.5 1.0];
    cfg.baseline = [-0.5 -0.2];
end

%%
cfg          = [];
cfg.baseline = [-.2 -.05];
cfg.channel = 50;
 
ERP_intact_bl = ft_timelockbaseline(cfg,ERP_intact);
ERP_scrambled_bl   = ft_timelockbaseline(cfg,ERP_scrambled);

ft_singleplotER(cfg,ERP_intact_bl,ERP_scrambled_bl)
%% HGP

low_freq = 5;
high_freq = 100;
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

cfg.trials = find(data.trialinfo <= 20 & data.trialinfo ~= -100); % select 'intact' trials
TFR_intact = ft_freqanalysis(cfg,data);

cfg.trials = find(data.trialinfo > 20 & data.trialinfo ~= 100); % select 'scrambled' trials
TFR_scrambled   = ft_freqanalysis(cfg,data);

% create HGP as empty timelock structure with same dimensions as ERP, values will be filled in in the next steps
HGP_intact = rmfield(ERP_intact,{'trial','avg','var'});
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
 
%%
close all
fcorr = squeeze(freqcorr_intact(10,24,:,200));
fpow = squeeze(TFR_intact.powspctrm(10,24,:,200));
fn = fcorr.*fpow;
f = TFR_intact.freq;
plot(f, log10(fn))
hold on
plot(f, log10(fcorr), 'r')
plot(f, log10(fpow), 'g')

%% clear TFR*
% plot TFR
cfg = [];
cfg.baseline = [-0.5 -0.1]; 
cfg.baselinetype = 'absolute'; 	
cfg.logplot = 'yes';
cfg.zlim = [-6 6]; %[-400000 300000];	
cfg.channel   = 10; % top figure
figure;ft_singleplotTFR(cfg, TFR_intact);
figure;ft_singleplotTFR(cfg, TFR_scrambled);

%%
cfg           = [];
cfg.parameter = 'avg';
cfg.xlim      = [-.3 .7];
cfg.channel   = 61; 
figure, ft_singleplotER(cfg,HGP_intact_bl,HGP_scrambled_bl)