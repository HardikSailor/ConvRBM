clear all; close all; clc;

startup;
% melfeats = load('E:\TIMIT_features_new\mfs\dr1_mcpm0_si1824_tr.wav.mfs');
load('sorted_Aurora4_Adam_dropout0p3_LRELU_40128_21052017.mat');
[sig2,fs] = audioread('E:\speech_data\ASV2017\eval\E_1008250.wav');
%       sig = pre_processing(sig,fs);
sig = bsxfun(@minus, sig,mean(sig));
sig = bsxfun(@rdivide, sig,(std(sig)));

ws = 128;
spacing = 16;
 winlen = fs*0.025; % 25 ms window length
 winhop = fs*0.010; % 10 ms window hop (overlap)
  tx = (0:length(sig)-1)./fs;

%    hbias_vec = zeros(80,1);
%  imdatatest = hardik_trim_audio_for_spacing_fixconv(sig, ws, spacing);
 [poshidexp,poshidprobs2] = tirbm_inference_fixconv_1d(imdatatest, W, hbias_vec, pars);

% [poshidexp,poshidprobs2] = perfectrecon_tirbm_inference(imdatatest, W, 0.01.*hbias_vec, pars);
 [poshidstates poshidprobs] = my_tirbm_sample_multrand_1d(poshidexp, spacing);
recon_data = tirbm_reconstruct_LB_fixconv_1d(poshidprobs, W, pars);

%% spectrograms
%   tx2 = (0:length(imdatatest)-1)./fs;
% STFT = log(abs(spectrogram(sig,winlen,winhop,[],fs)));
% STFT_recon = log(abs(spectrogram(recon_data,winlen,winhop,[],fs)));
%% visualizing reconstruction using spectrograms
% left= 1;
% bottom1=0.5;
% bottom2=0.1;
% width=1;
% height=1; % which is also bottom1-bottom2

% subplot(221)
% % axes('Position',[left bottom2 width height])
% plot(tx,sig);axis tight; xlabel('Time (ms)');ylabel('Amplitude'); title('Original Speech');
% subplot(222)
% % axes('Position',[left bottom2 width height])
% plot(tx2,recon_data);axis tight;xlabel('Time (ms)');ylabel('Amplitude');title('Reconstructed Speech');
% subplot(223)
% % axes('Position',[left bottom2 width height])
% imagesc(tx,[0:2:8],STFT); axis xy;axis tight;ylabel('Freq.(kHz)');xlabel('Time (ms)');
% subplot(224)
% % axes('Position',[left bottom2 width height])
% imagesc(tx2,[0:2:8],STFT_recon);axis xy;axis tight;ylabel('Freq.(kHz)')
% xlabel('Time (ms)');
%  tightfig;
% spaceplots([0.02 0.02 0.02 0.02], [.02 .02]);
%% reconstruction of small segment
ss1 = imdatatest(40000:40900);
ss2 = recon_data(40000:40900);
residual = ss1-ss2;
  txs = (0:length(ss1)-1)./fs;
subplot(311)
plot(txs,ss1);axis tight;ylabel('Amplitude')
subplot(312)
plot(txs,ss2);axis tight;ylabel('Amplitude')
subplot(313)
plot(txs,residual);axis tight;xlabel('Time (ms)');ylabel('Amplitude')