function [tfr, time, freq] =  tf_decompose(tc, fs, nfft)

% tf_decompose.m
% Given a time course (tc) with sample freq fs, 
% calculate the time-frequency decomposition between fl and fh


[tfr, freq, time] = spectrogram(tc,100,98,nfft,fs);