%%% single speaker data %%%%
clear all; clc; close all;
startup;
% load('ARCTIC_SLT_raw_speech_128_16_60_LB_TIMIT_V1b_w128_b60_pc1.000000e-001_sigpc0.1_p0.5_pl16_plambda109_sp102_5.000000e-003_eps0.0001_epsdecay0.001_l2reg1_bs50_0151102T104016.mat')
load('C:\Users\hardik\Documents\MATLAB\Unsupervised_learning\ConvRBM_speech_filterbank\sorted_filters\W3_SLT_60_filters.mat');
fs = 16000;
N = 128;
W = W3;
 W = reshape(W(end:-1:1, :),[128,60]);
wfft = abs(fft(W));
wfft = wfft(1:64,:);
%%% 31,34,56
tx = (0:length(W)-1)./fs;
tx = tx*1000;
tf = 0.001.*linspace(1,64,64)*(fs/N);
subplot = @(m,n,p) subtightplot (m, n, p, [0.09 0.1], [0.15 0.01], [0.1 0.05]);
subplot(2,2,1)
plot(tx,W(:,31)); axis tight;ylabel('Amplitude');
subplot(2,2,2)
plot(tf,wfft(:,31));axis tight;ylabel('Amplitude');
subplot(2,2,3)
plot(tx,W(:,46));axis tight; xlabel('Time (ms)'); ylabel('Amplitude');
subplot(2,2,4)
plot(tf,wfft(:,46));axis tight; xlabel('Frequency (kHz)'); ylabel('Amplitude');

