%%% ConvRBM feature extraction %%%
clear all; clc; close all;
startup;
% melfeats = load('E:\TIMIT_features_new\mfs\dr1_mcpm0_si1824_tr.wav.mfs');
% load('5040_raw_speech_128_32_60_LB_TIMIT_V1b_w128_b80_pc3.000000e-002_sigpc0.03_p10_pl16_plambda109_sp102_5.000000e-003_eps0.0001_epsdecay0.001_l2reg1_bs50_0150708T131409.mat');
 [sig,fs] = wavread('F:\Database_collection\timit_na_SA\dr1_mcpm0_si1824_tr.wav');
%       sig = pre_processing(sig,fs);
               sig = bsxfun(@minus, sig,mean(sig));
                sig = bsxfun(@rdivide, sig,(std(sig)));
                
 winlen = fs*0.025; % 25 ms window length
 winhop = fs*0.010; % 10 ms window hop (overlap)
  tx = (0:length(sig)-1)./fs;
  
  load('selectedW.mat');
  [l,b]=size(W4);
%   W4 = bsxfun(@minus,W4,mean(W4));
%   W4 = bsxfun(@rdivide,W4,max(W4));

res = [];
 for j=1:b

     filter_td =W4(:,j);
     y = conv(sig,filter_td,'same');
     res(j,:) = y;
     
     
 end
% res = res';
%  res = bsxfun(@minus, res, mean(res));
%  res = bsxfun(@rdivide, res, std(res));
% res = (res);
 Y = max_pooled((res),winlen,winhop);
 Y1 = log(Y);
  Y2 = log(Y+0.001);
  imagesc(tx,[0:2:8],Y1); axis xy;axis('tight');ylabel('Freq.(kHz)');figure;
  imagesc(tx,[0:2:8],Y2); axis xy;axis('tight');ylabel('Freq.(kHz)');
