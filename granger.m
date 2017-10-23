% Granger analysis

close all; clear all
addpath('/home/rajat/Dropbox/Teresa_Rajat/AP_EEG/fieldtrip-20170618/')

runID = 'A';
seg_type = 'movie'; %'camera'; % choose btn ['movie', 'camera']
config = 'filter_nobaseline';
load(['run_' runID '_seg_' seg_type '_' config '.mat']);
dataA = data;

runID = 'B';
seg_type = 'movie'; %'camera'; % choose btn ['movie', 'camera']
load(['run_' runID '_seg_' seg_type  '_' config '.mat']);
dataB = data;

runID = 'C';
seg_type = 'movie'; %'camera'; % choose btn ['movie', 'camera']
load(['run_' runID '_seg_' seg_type  '_' config '.mat']);
dataC = data;

data = ft_appenddata([], dataA, dataB, dataC);


parietalA= [58 59;59 60;60 61;61 62;53 54;54 55;55 56;48 49;41 42];
visualA= [51 52; 45 46;46 47;35 36; 30 31; 31 32;31 32;25 26;26 27];
parietalB=[58 59;59 60;23 24;24 25;53 54];
visualB=[61 62;55 56;49 50;29 30; 30 31;33 34;34 35;35 36; 36 37;38 39;39 40;40 41; 41 42];

nParietalA = length(parietalA);
nParietalB = length(parietalB);
nVisualA = length(visualA);
nVisualB = length(visualB);

label_parietalA = cell(1, nParietalA);
for i=1:nParietalA
    label_parietalA{i} = sprintf('A%d-E-A%d-E', parietalA(i,2), parietalA(i,1));
end

label_parietalB = cell(1, nParietalB);
for i=1:nParietalB
    label_parietalB{i} = sprintf('B%d-E-B%d-E', parietalB(i,2), parietalB(i,1));
end

label_visualA = cell(1, nVisualA);
for i=1:nVisualA
    label_visualA{i} = sprintf('A%d-E-A%d-E', visualA(i,2), visualA(i,1));
end


label_visualB = cell(1, nVisualB);
for i=1:nVisualB
    label_visualB{i} = sprintf('B%d-E-B%d-E', visualB(i,2), visualB(i,1));
end

%------------- Channels of interests -------------------------------------%


label_parietal = [label_parietalA label_parietalB];
label_visual = [label_visualA label_visualB];

nParietal = length(label_parietal);
nVisual = length(label_visual);


channelcmb = cell(nParietal*nVisual, 2);
for i=1:nParietal
    for j = 1:nVisual
        channelcmb{(i-1)*nVisual + j, 1} = label_parietal{i};
        channelcmb{(i-1)*nVisual + j, 2} = label_visual{j};
    end
end


%%

ntrials = length(data.trialinfo);

% %--------- Bands of interest -------------------------------------------
% alpha = 8:2:12;
% low_beta = 12:2:20;
% high_beta = 20:2:28;
% low_gamma = 28:4:64;
% high_gamma = 64:4:200;

% freqs of interest
%foi = unique([alpha low_beta high_beta low_gamma high_gamma]);

n_rpt = 12;
nfreq = 201;

nchannels = length(data.label);
nchannelcmb = length(channelcmb) * 4; % A-> B, B _> A, A -> A, B -> B for every pair
granger_alltrials = zeros(ntrials, nchannelcmb, nfreq);



for trial_i = 1:ntrials
    disp('-----------------------------------------------');
    disp(['Trial ' num2str(trial_i) ' being computed --->'])
    disp('-----------------------------------------------');
        
    tic()
    
    cfg = [];
    cfg.trials = trial_i;
    data_trial = ft_selectdata(cfg, data);

    cfg         = [];
    cfg.order   = 250;
    cfg.toolbox = 'bsmart';
    cfg.channelcmb = channelcmb;
    mdata       = ft_mvaranalysis(cfg, data_trial);
    mdata.label = data_trial.label;
    
    cfg        = [];
    cfg.method = 'mvar';
    mfreq      = ft_freqanalysis(cfg, mdata);
    
   
    cfg           = [];
    cfg.method    = 'granger';
    granger_out       = ft_connectivityanalysis(cfg, mfreq);
    
    granger_alltrials(trial_i, :, :) = granger_out.grangerspctrm;
    
    toc()
end

save(['granger_values' config '.mat'], 'granger_alltrials', 'granger_out');

