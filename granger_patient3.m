% Granger analysis

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


% Patient-3
% Selection based on max(ISCi,ISCs)>0.1, max 10 per category

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



% %--------- Bands of interest -------------------------------------------
% alpha = 8:2:12;
% low_beta = 12:2:20;
% high_beta = 20:2:28;
% low_gamma = 28:4:64;
% high_gamma = 64:4:200;

% freqs of interest
%foi = unique([alpha low_beta high_beta low_gamma high_gamma]);

%n_rpt = 12;
% 
% 
% nchannels = length(data.label);
% nchannelcmb = length(channelcmb);


data_intact = data;
data_intact.trial = data.trial(data.trialinfo <= 20);
data_intact.time  = data.time(data.trialinfo <= 20);
data_intact.sampleinfo  = data.sampleinfo(data.trialinfo <= 20,:);
data_intact.trialinfo = data.trialinfo(data.trialinfo <= 20);

cfg         = [];
cfg.order   = 20;
cfg.toolbox = 'bsmart';
cfg.channel = selectedChannels;
mdata       = ft_mvaranalysis(cfg, data_intact); %_trial);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

cfg           = [];
cfg.method    = 'granger';
granger_intact       = ft_connectivityanalysis(cfg, mfreq);


data_scrambled = data;
data_scrambled.trial = data.trial(data.trialinfo > 20);
data_scrambled.time  = data.time(data.trialinfo > 20);
data_scrambled.sampleinfo  = data.sampleinfo(data.trialinfo <= 20,:);
data_scrambled.trialinfo = data.trialinfo(data.trialinfo > 20);

cfg         = [];
cfg.order   = 20;
cfg.toolbox = 'bsmart';
cfg.channel = selectedChannels;
mdata       = ft_mvaranalysis(cfg, data_scrambled); %_trial);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

cfg           = [];
cfg.method    = 'granger';
granger_scrambled       = ft_connectivityanalysis(cfg, mfreq);

ntrials = length(data.trialinfo);
nfreq = 201; %0:1:Fs/2
granger_alltrials = zeros(ntrials, nSelectedChannels, nSelectedChannels, nfreq);


% Do the granger per trial
for trial_i = 1:ntrials
    disp('-----------------------------------------------');
    disp(['Trial ' num2str(trial_i) ' being computed --->'])
    disp('-----------------------------------------------');
        
    tic()
    
    cfg = [];
    cfg.trials = trial_i;
    data_trial = ft_selectdata(cfg, data);

    cfg         = [];
    cfg.order   = 20;
    cfg.toolbox = 'bsmart';
    %cfg.channelcmb = channelcmb;
    cfg.channel = selectedChannels;
    mdata       = ft_mvaranalysis(cfg, data_trial);
    
    
    cfg        = [];
    cfg.method = 'mvar';
    mfreq      = ft_freqanalysis(cfg, mdata);
    
   
    cfg           = [];
    cfg.method    = 'granger';
    granger_out       = ft_connectivityanalysis(cfg, mfreq);
    
    granger_alltrials(trial_i, :, :, :) = granger_out.grangerspctrm;
    
    toc()
end

save(['P3_granger_values' config '.mat'], 'granger_alltrials', 'granger_intact', 'granger_scrambled', 'granger_out');

%%
k=1;
for i=1:16
    for j=1:16
        subplot(16,16,k)
        
        plot(granger_intact.freq(8:90), log(squeeze(granger_intact.grangerspctrm(i,j,8:90))))
        %axis([1 90 0 0.0005])
        hold on
        plot(granger_intact.freq(8:90), log(squeeze(granger_scrambled.grangerspctrm(i,j,8:90))))
        
        k = k + 1;
    end
end


%%

lf = 40;
hf = 90;
k=1;
for i=1:3
    for j=6:15
        subplot(3,10,k)
        plot(granger_intact.freq(lf:hf), squeeze(granger_intact.grangerspctrm(i,j,lf:hf)))
        axis([40 90 0 1e-3])
        hold on
        plot(granger_intact.freq(lf:hf), squeeze(granger_scrambled.grangerspctrm(i,j,lf:hf)))
        k=k+1;
    end
end

figure

lf = 2;
hf = 40;
k=1;
for i=6:15
    for j=1:3
        subplot(10,3,k)
        plot(granger_intact.freq(lf:hf), squeeze(granger_intact.grangerspctrm(i,j,lf:hf)))
        axis([2 40 0 1e-3])
        hold on
        plot(granger_intact.freq(lf:hf), squeeze(granger_scrambled.grangerspctrm(i,j,lf:hf)))
        k=k+1;
    end
end


%%

intact = granger_alltrials(data.trialinfo <=20, :, :, :);
scrambled = granger_alltrials(data.trialinfo >20, :, :, :);
mIntact = mean(intact(:,:,:,40:90),4);
mScram = mean(scrambled(:,:,:,40:90),4);

i=1;
j=10;

tMat = zeros(16,16);
pMat = zeros(16,16);
for i=1:16
    for j=1:16

X = mIntact(:,i,j);
Y = mScram(:,i,j);

[H,P,CI,STATS] = ttest2(X,Y,'tail','left');

tMat(i,j) = STATS.tstat;
pMat(i,j) = P;

    end
end

imagesc(tMat.*(pMat<0.05))
figure
imagesc(pMat.*(pMat<0.01))
%%



intact = granger_alltrials(data.trialinfo <=20, :, :, :);
scrambled = granger_alltrials(data.trialinfo >20, :, :, :);
mIntact = mean(intact(:,:,:,8:30),4);
mScram = mean(scrambled(:,:,:,8:30),4);

i=1;
j=10;

tMat = zeros(16,16);
pMat = zeros(16,16);
for i=1:16
    for j=1:16

X = mIntact(:,i,j);
Y = mScram(:,i,j);

[H,P,CI,STATS] = ttest2(X,Y,'tail','both');

tMat(i,j) = STATS.tstat;
pMat(i,j) = P;

    end
end

imagesc(tMat.*(pMat<0.05))














