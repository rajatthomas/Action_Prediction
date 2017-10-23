% Compute Fourier transform over all time

% load the data set

%%
clear all; close all
runID = 'D1'; % 'A2', similarly 'B' or 'C'

%load(['run_' runID '_seg_camera.mat']);
%load(['run_' runID '_seg_movie.mat']);
load(['P3_run_' runID '_seg_movie_nofilter_nobaseline.mat']);

%L = 512; % Points to consider for the FFT

%%
% % Fs = 1000;            % Sampling frequency
% % T = 1/Fs;             % Sampling period
% % X = delta(20,:);
% % L = length(X);             % Length of signal
% % t = (0:L-1)*T;        % Time vector
% % 
% % tot_P1 = zeros(1,L/2+1);
% % for i=1:87
% %     
% % X = delta(i,:);
% % Y=fft(X);
% % P2 = abs(Y/L);
% % P1 = P2(1:L/2+1);
% % P1(2:end-1) = 2*P1(2:end-1);
% % tot_P1 = tot_P1 + P1;
% % end
% % 
% % 
% % f = Fs*(0:(L/2))/L;
% % loglog(f,tot_P1)
% % title('Single-Sided Amplitude Spectrum of X(t)')
% % xlabel('f (Hz)')
% % ylabel('|P1(f)|')

%%

Fs = data.fsample;            % Sampling frequency 
[nchannels, ntimes] =size(data.trial{1});
ntrials = length(data.trial); %-2; % last two are baselines if present

L = 512;
f = Fs*(0:floor(L/2))/L;
fcorr = repmat(f',[1 nchannels]);

POW = zeros(floor(L/2)+1, nchannels, ntrials);

for trial_i = 1:ntrials
    Y = fft(data.trial{trial_i}(:,5000:end)',L);
    P2 = abs(Y/L).^2;
    P1 = P2(1:floor(L/2)+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    POW(:,:,trial_i) = P1;
end

trial_intact = data.trialinfo <= 20 & data.trialinfo ~= -100;
trial_scram = data.trialinfo > 20 & data.trialinfo ~= 100;

mean_pow_intact = mean(POW(:,:,trial_intact),3).*fcorr;
mean_pow_scram = mean(POW(:,:,trial_scram),3).*fcorr;
mean_trial_pow = mean(POW,3).*fcorr;


% Repeat for the baseline

baseline_pre  = data.trial{end-1};
baseline_post = data.trial{end};
baseline = [baseline_post];% baseline_post];

nChunks = length(1:L:(size(baseline,2)-L));
baseline_pow = zeros(floor(L/2)+1, nchannels, nChunks);

i_chunk=1;
for chunk_i = 1:L:(size(baseline,2)-L)
   X = baseline(:,chunk_i:(chunk_i+L-1));
   Y  = fft(X');
   P2 = abs(Y/L).^2;
   P1 = P2(1:floor(L/2)+1,:);
   P1(2:end-1,:) = 2*P1(2:end-1,:);
   baseline_pow(:,:,i_chunk) = P1;
   i_chunk = i_chunk + 1;
end

baseline_pow = mean(baseline_pow,3).*fcorr;


%%
ch_i = 24;
SMP=10;
for ch_i = 1:77
%plot(f,smooth(log10(mean_trial_pow(:,ch_i)'),SMP))
%hold on

plot(f,smooth(log10(mean_pow_intact(:,ch_i)'),SMP))
hold on
%plot(f,smooth(log10(mean_pow_scram(:,ch_i)'),SMP))
%plot(f,((abs(mean_pow_intact(:,ch_i))'-abs(mean_pow_scram(:,ch_i))')))
%hold on
%plot(f,smooth(log10(baseline_pow(:,ch_i)'),SMP))
title(['Pow. Spect. <Trial> Ch# ' data.label{ch_i}])
xlabel('f (Hz)')
ylabel('|P(f)|')
%pause
%close all
end


%%

ch_i=35;
plot(f,smooth(log10(mean_pow_scram(:,ch_i)'),SMP),'r-')
hold on
plot(f,smooth(log10(A1_mean_scram(:,ch_i)'),SMP),'r.')

plot(f,smooth(log10(mean_pow_intact(:,ch_i)'),SMP),'b-')
plot(f,smooth(log10(A1_mean_intact(:,ch_i)'),SMP),'b.')


%plot(f,smooth(log10(mean_pow_scram(:,ch_i)'),SMP))
%plot(f,((abs(mean_pow_intact(:,ch_i))'-abs(mean_pow_scram(:,ch_i))')))
%hold on
%plot(f,smooth(log10(baseline_pow(:,ch_i)'),SMP))
title(['Pow. Spect. <Trial> Ch# ' data.label{ch_i}])
xlabel('f (Hz)')
ylabel('|P(f)|')
hold on
%plot(f,smooth(log10(mean_pow_scram(:,ch_i)'),SMP))

%%

diff_freq = mean_pow_intact - mean_pow_scram;
corr_freq_mat = corrcoef(diff_freq');


ch_i = 24;
SMP=10;
for ch_i = 1:98

    
loglog(f,smooth(abs(mean_pow_scram(:,ch_i)'- mean_pow_intact(:,ch_i)'),SMP))
hold on
%plot(f,smooth(log10(baseline_pow(:,ch_i)'),SMP))
%title(['Pow. Spect. <Trial> Ch# ' data.label{ch_i}])
%xlabel('f (Hz)')
%ylabel('|P(f)|')
%pause
%close all
end
