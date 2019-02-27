 clear all; clc;
startup;
% melfeats = load('E:\TIMIT_features_new\mfs\dr1_mcpm0_si1824_tr.wav.mfs');
%    load('ARCTIC_SLT_raw_speech_128_16_60_LB_TIMIT_V1b_w128_b60_pc1.000000e-001_sigpc0.1_p0.5_pl16_plambda109_sp102_5.000000e-003_eps0.0001_epsdecay0.001_l2reg1_bs50_0151102T104016.mat')
%  load('ARCTIC_SLT_raw_speech_128_16_60_LB_TIMIT_V1b_w128_b60_pc1.000000e-001_sigpc0.1_p0.5_pl16_plambda109_sp102_5.000000e-003_eps0.0001_epsdecay0.001_l2reg1_bs50_0151101T202405.mat')
% load('ARCTIC_BDL_128_relu_reg_W40.mat');
%  load('5040_raw_speech_128_16_60_LB_TIMIT_V1b_w128_b60_pc3.000000e-002_sigpc0.03_p10_pl16_plambda109_sp102_5.000000e-003_eps0.0001_epsdecay0.001_l2reg1_bs50_0150902T115354.mat')
%    load('TIMIT_scrambled_4250_dropout0p3_40128_19082017.mat')
%  load('ProsodyGuj2_6460.mat');
%  load('Pro_ASR_comb_64_16_40_LB_TIMIT_V1b_w64_b60_pc5.000000e-02_sigpc0.05_p1_pl16_plambda109_sp102_5.000000e-03_eps0.0001_epsdecay0.001_l2reg1_bs50_0160405T200655.mat');
%  load('TIMIT_raw_speech_128_16_40_LB_TIMIT_V1b_w128_b40_pc3.000000e-002_sigpc0.03_p10_pl16_plambda109_sp102_5.000000e-003_eps0.0001_epsdecay0.001_l2reg1_bs50_0150930T103511.mat');
% load('WSJ0_raw_speech_128_40_LB_WSJ0_V1b_w128_b40_pc1.000000e-01_sigpc0.1_p1_pl16_plambda109_sp102_5.000000e-03_eps0.0001_epsdecay0.001_l2reg1_bs50_0151105T225627.mat');
%    load('Aurora4_train_12860all_traintest_corrected.mat');
%   load('Replay_Adam_dropout_preemphasis_60160_LReLU_05092017.mat');
%  load('Replay_Adam_dropout_preemphasis_40128_LReLU_06092017.mat')
%  load('Replay_Adam_dropout_preemphasis_6096_LReLU_08092017.mat')

%    load('Aurora4_Adam_dropout0p3_40128_17122016.mat')
%  load('sorted_MSIS2018Gujarati_Adam_AD0p3_40128_LReLU_21012018.mat')
%      load('MSIS2018Telugu_WSJ1TL_Adam_AD0p3_40128_LReLU_25012018.mat')
%   load('sorted_Aurora4_Adam_dropout0p3_LRELU_40128_21052017.mat')
%   load('sorted_ReplayTrainAll_Adam_LReLU_60128_b0p5_26012017.mat')
%   load('NewReplayTrain_nopreemphasis_Adam_AD0p3_80128_LReLU_10052018.mat')
%    load('sorted_replayTrain_Adam_LReLU_preemphasis_60128_07092017.mat')
 load('E:\Codes_phd\Cleaned_ConvRBM_dropout_Adam_full\Sound_classification_learned_weights\ESC50_Adam_dropout0p5_40176_18122016.mat')
%    load('sorted_Replay_Adam_dropout_preemphasis_6096_LReLU_08092017.mat');
%              W = W3;
%       W = reshape(W(end:-1:1, :),[128,40]);
%[sig1,fs] = audioread('E:\speech_data\AVS2015\p225_026.wav');
% [sig1,fs] = audioread('E:\speech_data\ASV2017\eval\E_1008240.wav');
% [sig2,fs] = audioread('E:\speech_data\ASV2017\eval\E_1008250.wav');
%[sig1,fs] = audioread('E:\speech_data\ASV2017\eval\E_1010305.wav');
% [sig1,fs] = audioread('E:\speech_data\ASV_2017_updated_database\ASVspoof2017_V2_train\T_1000392.wav');
[sig2,fs] = audioread('E:\speech_data\ESC_all_22kHz_single\1-977-A.wav');
% [sig1,fs] = audioread('/media/hardik/PhD_work/speech_data/Aurora4_test_wavfiles/440c020e9.wav');
%   sig1 = shufflewins(sig2,round(fs*.004),round(fs*.50));
sig2 = bsxfun(@minus, sig2,mean(sig2));
sig2 = bsxfun(@rdivide, sig2,(std(sig2)+0.0001));
%       sig1 = pre_processing(sig1,fs);
% sig1 = sig2(:,1);
% sig1 = resample(sig1,160,441);
% sig = bsxfun(@minus, sig1,mean(sig1));
% sig = bsxfun(@rdivide, sig,(std(sig)+0.0001));
%  sig = sig(30000:90000);               
 winlen = round(fs*0.025); % 25 ms window length
 winhop = round(fs*0.010); % 10 ms window hop (overlap)
  tx = (0:length(sig2)-1)./fs;
%%% spectrogram, Mel-spectrum and Gammatone spectrum %%%%%%%%%
x_frames = enframe(sig2,hamming(winlen),winhop);
STFT = abs(fft(x_frames',512));
STFT = STFT(1:257,:);
% STFT = (abs(spectrogram(sig,winlen,winhop,512,fs)));
 STFT = (STFT.^2);

melbank = filtbankm(40,512,fs,100,fs/2,'m');
Melspec = STFT'*melbank';
Melspec = log(Melspec);
% Melspec = (Melspec);
%     Melspec = bsxfun(@minus,Melspec,mean(Melspec));
%    Melspec = bsxfun(@rdivide,Melspec,std(Melspec));
  Melspec = Melspec';
%  gspec = gammatonegram(sig,fs,winlen/fs,winhop/fs,40,0,fs/2);
%  gspec = (log(gspec));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PARAMETERS
% B = 96;
% fmax = fs/2;
% fmin = fmax/2^9;
% d = 16;
% cf = 39;
% ZsdD = 'ZsdD';
% 
% %% COMPUTE CQCC FEATURES
% [CQcc, LogP_absCQT, TimeVec, FreqVec, Ures_LogP_absCQT, Ures_FreqVec] = ...
%     cqcc(sig, fs, B, fmax, fmin, d, cf, ZsdD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bssis selection based in variance                
v1 = var(W);
p1 = sort(v1);
p2 = find(v1>0);
W2 = W(:,p2);
% bais selection based on L2norm
% for k=1:80
%     n1(k) = norm(W(:,k));
% end
% n2 = sort(n1);
% n3 = find(n2>0.1113);
% W2 = W(:,n3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wfft = abs(fft(W2,201));
wfft=wfft(1:100,:);
% locations = zeros(1,p2);
for i=1:length(p2)
     filter_fd = wfft(:,i)./max(wfft(:,i));
     fp = find(filter_fd==1);
     locations(i) = fp;
%     filter_td =W2(:,i);
%      y = conv(sig,filter_fd,'same');
%      res(i,:) = y;

end
[s1,I1] = sort(locations);
% [us1,I2] = unique(s1);
% [s2,I2] = sort(us1);
wfft3 = wfft(:,I1);
wfft4 = bsxfun(@rdivide,wfft3,max(wfft3));
W3 = W2(:,I1);
 W4 = W3;
%   W4 = bsxfun(@minus,W3,mean(W3));
 hbias_vec = hbias_vec(I1);
% hbias_vec = [];

%%
% res = [];
% matlabpool('open',3);
% tic;
 for j=1:length(s1)

     filter_td =W4(:,j);
%      filter_td = filter_td./norm(filter_td);
     y = conv(sig2,filter_td,'same');
     res(j,:) = y;
     
     
 end
%  toc;
% pars = [];
% hbias_vec = zeros(size(hbias_vec));
% dropoutrate = 0;
%     [poshidactvalues,poshidsample,DropMask] = test_ConvRBM_sample_multrand_1d(sig, W, hbias_vec, pars,dropoutrate);
%     negdata = ConvRBM_reconstruct_LB_fixconv_1d(poshidactvalues, W, vbias_vec, pars);
%           negdata = bsxfun(@minus, negdata,mean(negdata));
%                 negdata = bsxfun(@rdivide, negdata,(std(negdata)+0.0001));
%   [feats_s] = ConvRBM_feat_extract(sig,fs,W3,winlen,winhop,hbias_vec);
  [feats_N] = ConvRBM_feat_extract(sig2,fs,W3,winlen,winhop);
% [Mel_s_feats] = Mel_spec_temp(sig,fs,winlen,winhop,40);
[Mel_N_feats] = Mel_spec_temp(sig2,fs,winlen,winhop,40);

% [feats1] = AMFM_ConvRBM_feat_extract(sig,fs,W3,winlen,winhop,hbias_vec);

%  matlabpool('close');
%  res = res';
%   res = bsxfun(@minus, res, mean(res));
%    res = bsxfun(@rdivide, res, std(res));
% %  res = log(res+0.0001);
% res = res';
%  Y = max_pooled((res),winlen,winhop);
%          Y = Y.^(1/15);
%            Y = log(Y+0.0001);
 %% plotting
%  TEO_feats = TEO_feats;
 figure;
 subplot(311)
plot(tx,sig2);axis('tight');ylabel('Amplitude'); 
% subplot(412)
% imagesc(tx,[0:2:8],Melspec); axis xy;axis('tight');ylabel('Freq.(kHz)');
 subplot(312)
 imagesc(tx,[0:2:8],feats_N);axis xy;axis('tight');ylabel('Freq.(kHz)');
 subplot(313)
 imagesc(tx,[0:2:8],Mel_N_feats);axis xy;axis('tight');ylabel('Freq.(kHz)');
xlabel('Time (s)');
%  figure;
%%%%%%%%%% plot time scambling %%%%%%%%%%%%%%%%
%%% Utt: Spring Street is straight ahead. for wav: dr3_fgrw0_sx72_tr.wav
%  figure;
%  subplot(321)
% plot(tx,sig2);axis('tight');ylabel('Amplitude'); 
%  subplot(322)
% plot(tx,sig);axis('tight');
% % subplot(412)
% % imagesc(tx,[0:2:8],Melspec); axis xy;axis('tight');ylabel('Freq.(kHz)');
%  subplot(323)
%  imagesc(tx,[0:2:8],Mel_N_feats);axis xy;axis('tight');ylabel('Freq.(kHz)');
%  subplot(324)
%  imagesc(tx,[0:2:8],Mel_s_feats);axis xy;axis('tight');
%  subplot(325)
%  imagesc(tx,[0:2:8],feats_N);axis xy;axis('tight');ylabel('Freq.(kHz)');
% xlabel('Time (s)');
%  subplot(326)
%  imagesc(tx,[0:2:8],feats_s);axis xy;axis('tight');
% xlabel('Time (s)');
% M=60; N=60; L = 22;
%   CC11 = dct(feats2);
% %   CC = CC(1:13,:);
%    CC = cep_coeff(feats1,M,N,L);
%    D1 = DeltasFunction(CC,4);
%    D2 = DeltasFunction(D1,4);
% CC2 = deltas(CC);
% CC3 = deltas(CC2);
% ConvRBMCC = [CC;CC2;CC3];
% ConvRBMCC2 = [CC;D1;D2];
% CCmvn = CC';
% DC1 = bsxfun(@minus, CCmvn,mean(CCmvn));
% % DC2 = bsxfun(@rdivide, DC1,std(DC1));
% hh1 = DeltasFunction(DC1,4);
% hh2 = DeltasFunction(hh1,4);
% featers = [DC1';hh1';hh2'];
% % imagesc(featers)
% % plot(ans)
% % W3 = W(,);
%      set(gca,'LooseInset',get(gca,'TightInset'));
%  set(gca,'Visible','off');
%  set(gca,'XTickLabel','');set(gca,'YTickLabel','');