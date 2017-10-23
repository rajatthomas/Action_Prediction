% Difference in power spectrum

%%
addpath('/home/rajat/Dropbox/Teresa_Rajat/AP_EEG/fieldtrip-master/')
ft_defaults
%--------- Bands of interest -------------------------------------------
alpha = 8:2:12;
low_beta = 12:2:20;
high_beta = 20:2:28;
low_gamma = [28:2:46 54:2:64];
high_gamma = [64:2:98 102:2:148 152:2:196];



%% load the super intact: (no camera changes or time warps)

run_type = {'super_intact', 'time_warped', 'intact', };

for run_i = 1:2
   for type_i=1:4 
    data = get_data(run_i, run_type{type_i});
    freq_spectrum = get_tf(data, 'power_spectrum');
   end
end
    


%% Power of super intact run1

super_intact = [1, 5, 8, 9, 10, 11, 12, 16, 20];

data = data_super_intact_run1;

Fs = 400;            % Sampling frequency 
[nchannels, ntimes] =size(data.trial{1});
ntrials = length(data.trial); %-2; % last two are baselines if present

L = 512;
f = Fs*(0:floor(L/2))/L;
fcorr = repmat(f',[1 nchannels]);

POW_super_intact_run1 = zeros(floor(L/2)+1, nchannels, ntrials);

for trial_i = 1:ntrials
    Y = fft(data.trial{trial_i}',L);
    P2 = abs(Y/L).^2;
    P1 = P2(1:floor(L/2)+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    POW_super_intact_run1(:,:,trial_i) = P1;
end


%% Power of super intact run2

super_intact = [1, 5, 8, 9, 10, 11, 12, 16, 20];

data = data_super_intact_run2;

Fs = 400;            % Sampling frequency 
[nchannels, ntimes] =size(data.trial{1});
ntrials = length(data.trial); %-2; % last two are baselines if present

L = 512;
f = Fs*(0:floor(L/2))/L;
fcorr = repmat(f',[1 nchannels]);

POW_super_intact_run2 = zeros(floor(L/2)+1, nchannels, ntrials);

for trial_i = 1:ntrials
    Y = fft(data.trial{trial_i}',L);
    P2 = abs(Y/L).^2;
    P1 = P2(1:floor(L/2)+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    POW_super_intact_run2(:,:,trial_i) = P1;
end
%% Power of scrambled run1

scrambled = [31,36,25,28,21,40,30,29,32];
scram2super = [6, 8, 2, 3, 1, 9, 5, 4, 7]; % because super_intact and scrambled movies are not in the same order

data = data_scrambled_run1;

Fs = 400;            % Sampling frequency 
[nchannels, ntimes] =size(data.trial{1});
ntrials = length(data.trial); %-2; % last two are baselines if present

L = 512;
f = Fs*(0:floor(L/2))/L;
fcorr = repmat(f',[1 nchannels]);

POW_scrambled_run1 = zeros(floor(L/2)+1, nchannels, ntrials);

for trial_i = 1:ntrials
    Y = fft(data.trial{trial_i}',L);
    P2 = abs(Y/L).^2;
    P1 = P2(1:floor(L/2)+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    POW_scrambled_run1(:,:,trial_i) = P1;
end

%%
scrambled = [31,36,25,28,21,40,30,29,32];
scram2super = [6, 8, 2, 3, 1, 9, 5, 4, 7]; % because super_intact and scrambled movies are not in the same order

data = data_scrambled_run2;

Fs = 400;            % Sampling frequency 
[nchannels, ntimes] =size(data.trial{1});
ntrials = length(data.trial); %-2; % last two are baselines if present

L = 512;
f = Fs*(0:floor(L/2))/L;
fcorr = repmat(f',[1 nchannels]);

POW_scrambled_run2 = zeros(floor(L/2)+1, nchannels, ntrials);

for trial_i = 1:ntrials
    Y = fft(data.trial{trial_i}',L);
    P2 = abs(Y/L).^2;
    P1 = P2(1:floor(L/2)+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    POW_scrambled_run2(:,:,trial_i) = P1;
end


%% plots
POW_scrambled_order_run1 = POW_scrambled_run1(:, :, scram2super);
%POW_diff = log10(POW_scrambled_order_run1) - log10(POW_super_intact_run1);


mean_POW_scram_run1 = log10(mean(POW_scrambled_order_run1,3).*fcorr);
mean_POW_super_intact_run1 = log10(mean(POW_super_intact_run1,3).*fcorr);

POW_scrambled_order_run2 = POW_scrambled_run2(:, :, scram2super);
%POW_diff = log10(POW_scrambled_order_run1) - log10(POW_super_intact_run1);


mean_POW_scram_run2 = log10(mean(POW_scrambled_order_run2,3).*fcorr);
mean_POW_super_intact_run2 = log10(mean(POW_super_intact_run2,3).*fcorr);


%mean_POW_diff = mean(POW_diff,3);
%% differences

parietalA= [58 59;59 60;60 61;61 62;53 54;54 55;55 56;48 49;41 42];
visualA= [51 52; 45 46;46 47;35 36; 30 32;31 32;25 26;26 27];
parietalB=[58 59;59 60;23 24;24 25;53 43];
visualB=[61 62;55 56;49 50;29 30; 30 31;33 34;34 35;35 36; 36 37;38 39;39 40;40 41; 41 42];

cluster = parietalA;

for i=1:length(cluster)
    
    
    channel = sprintf('A%d-B39-A%d-B39', cluster(i,2), cluster(i,1));
    
    for ch_i=1:77
        if strcmp(channel, data.label(ch_i))
            break
        end
    end
    %p = smooth(mean_POW_diff(:,ch_i),10);
    I = smooth(mean_POW_super_intact_run1(:,ch_i),5);
    S = smooth(mean_POW_scram_run1(:,ch_i),5);
   
    plot(f(5:160),I(5:160),'b')
    hold on
    plot(f(5:160), S(5:160),'r')
    
    I = smooth(mean_POW_super_intact_run2(:,ch_i),5);
    S = smooth(mean_POW_scram_run2(:,ch_i),5);
    
    plot(f(5:160),I(5:160),'b')
    
    plot(f(5:160), S(5:160),'r')
    
    title(data.label(ch_i))
    pause
    close all
end

plot



%%


for i=1:100
    loglog(f, smooth(POW_super_intact(:,i, trial_i), 10))
    hold on
    loglog(f, smooth(POW_scrambled_order_run1(:,i, trial_i), 10))
    pause
    close all
end




