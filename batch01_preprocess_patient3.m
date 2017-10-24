% Preprocessing pipeline for ECoG data of Patient-3 (no triggers this time)
%
% Step-1: Create the bi-polar montage
% Step-2: Preprocess the data (removing the 1/f noise) 
% Step-3: Segment into movies 
clear all; close all
addpath('/home/rajat/Dropbox/Teresa_Rajat/AP_EEG/fieldtrip-20170618/')

% run
patientNum = 3;
runID{1} = 'A1'; runID{2} = 'B1'; runID{3} = 'C1';
runID{4} = 'A2'; runID{5} = 'B2'; runID{6} = 'C2';

runID{7} = 'D1'; runID{8} = 'E1';
runID{9} = 'D2'; runID{10} = 'E2';

% segment at which level?
%seg_type{1} = 'camera';
seg_type{1} = 'movie';
%seg_type{2} = 'camera';

n_runs = length(runID);
n_segs = length(seg_type);

POST_BUFFER = 2000; % Baseline considered 2 seconds after last camera change

filtered = false;
baseline = false;
for run_i = 1:n_runs
    for seg_i = 1:n_segs
    
        disp(['P' num2str(patientNum) '_run ' runID{run_i} ' segmenting ' seg_type{seg_i}])
    %--------------------- Load the dataset ------------------------------------
 
    if filtered
        if baseline
            output_file = sprintf('P%d_run_%s_seg_%s_filter_baseline.mat', patientNum, runID{run_i}, seg_type{seg_i});
        else
            output_file = sprintf('P%d_run_%s_seg_%s_filter_nobaseline.mat', patientNum, runID{run_i}, seg_type{seg_i});
        end
        data_dir = '/home/rajat/Dropbox/JapanECoG/raw_data/Patient3/raw_LP300Hz_HP_p08Hz_noHumFilter';
        data_file = fullfile(data_dir,['sequential_action_movie_' runID{run_i} '.mat']);
        load(data_file) % loads raw (channels X timepoints) and hdr (header info)
    else
        if baseline
            output_file = sprintf('P%d_run_%s_seg_%s_nofilter_baseline.mat', patientNum, runID{run_i}, seg_type{seg_i});
        else
            output_file = sprintf('P%d_run_%s_seg_%s_nofilter_nobaseline.mat', patientNum, runID{run_i}, seg_type{seg_i});
        end
        data_dir = '/home/rajat/Dropbox/JapanECoG/raw_data/Patient3/raw_LP300Hz_HP_p08Hz_noHumFilter';
        % trigger channel 63 from the filtered case is not present.
        data_file = fullfile(data_dir,['sequential_action_movie_' runID{run_i} '.mat']);
        load(data_file) % loads raw (channels X timepoints) and hdr (header info)
    end
    
    
    % -------------- create bipolar montage ------------------------------------
    
    % Refer to ElectrodeLocations_Patient02.pdf for Montage creation.
    % NOTE: electrodes labelled 'A' and 'B' are totally symmetric.

    
    
    bipolar_montage = diff(raw); % subtracts n+1^th channel from n^th channel

    % because there are different sets say (A1 to A4) and (A5 to A8) and so on,
    % some bipolar channels like A5-A4 or A23-A22 does not make sense. Thus, we
    % remove them.
    % For example, after diff(raw), the 4th row contains A5-A4 and 22nd row
    % A23-A22 and so on.
    
    
    to_remove_A = [8,16,24,32,40]; % For all A channels
    to_remove_B = [0,4,12,20,24] + 60; % for B channels 

    to_remove = [to_remove_A, to_remove_B];
    %to_remove = [];
    
    labels = cell(1,length(hdr.label)-1);
    chantype = cell(1,length(hdr.label)-1);
    chanunit = cell(1,length(hdr.label)-1);
    for ch_i = 1:length(hdr.label)-1
        labels{ch_i} = [hdr.label{ch_i+1} '-' hdr.label{ch_i}];
        chantype{ch_i} = 'eeg';
        chanunit = 'uV';
    end

    % remove montages that do not make sense (sensors too far apart)
    bipolar_montage(to_remove,:) = [];
    labels(to_remove) = [];
    chantype(to_remove) = [];
    

    % --------------------- convert to fieldtrip format ---------------------------
    
    % Get trial definitions
    cfg = [];
    cfg.runID = runID{run_i}; 
    cfg.patientNum = patientNum;
    cfg.Fs            = hdr.Fs;
    if (strcmp(seg_type{seg_i}, 'camera'))
        cfg.trialdef.pre = 0.5;
        cfg.trialdef.post= 0.7;
        trl = trialdefs_camera_changes(cfg);
    else % movie
        cfg.trialdef.pre = 2.0; %  secs before movie starts 
        cfg.trialdef.post= 40.0; % 40 secs of movie (shortest movie ~ 40secs)
        trl = trialdefs_movies(cfg);
    end

    if baseline
        % Baseline trials is extracted w.r.t camera changes
        trl_baseline = trialdefs_camera_changes(cfg);
        pre_baseline_trl = [1 trl_baseline(1,1) 0 -100]; % -100 is the ID for pre_baseline
        post_baseline_trl= [trl_baseline(end,2)+POST_BUFFER size(raw,2) 0 100]; %+100 for post
        trl = [trl; pre_baseline_trl; post_baseline_trl];
    end
    
    % Construct FieldTrip format datafiles
    [hdr.nChans, hdr.nSamples] = size(bipolar_montage);
    hdr.label = labels;
    hdr.chantype = chantype;
    hdr.chanunit = chanunit;

    data = [];
    data.hdr     = hdr;
    data.trial{1}= bipolar_montage;
    data.fsample = hdr.Fs;
    data.label   = hdr.label;
    data.time{1} = (1:hdr.nSamples)./hdr.Fs; % in seconds
    data.sampleinfo = [1 hdr.nSamples];
    data.cfg        = [];

    cfg = [];
    cfg.trl = trl;
    data = ft_redefinetrial(cfg,data);
    sampleinfo = data.sampleinfo;
   
    % --------------------- Preprocessing the data --------------------------------


    % Resampling
    cfg = [];
    cfg.resamplefs = 400; % resample to 400 Hz
    cfg.detrend    = 'yes';
    cfg.demean     = 'yes';
    cfg.trials     = 'all';

    data = ft_resampledata(cfg, data);
    
    % define a new data.sampleinfo
    data.sampleinfo = round(sampleinfo.*data.fsample./data.hdr.Fs);
    
    
    % Filtering
    
    if filtered
        
        disp('FILTERING-----------');
        %cfg = [];
        %cfg.dftfilter     = 'yes';
        %cfg.dftfreq       = [50 100 150];
        cfg = [];
        cfg.bsfilter      = 'yes';
        cfg.bsfreq = [48 52];
        data              = ft_preprocessing(cfg,data);
        
        cfg = [];
        cfg.bsfilter      = 'yes';
        cfg.bsfreq = [98 102];
        data              = ft_preprocessing(cfg,data);
        
        cfg = [];
        cfg.bsfilter      = 'yes';
        cfg.bsfreq = [148 152];
        data              = ft_preprocessing(cfg,data);
        
    end
    
    save(output_file, 'data');
    end
end


