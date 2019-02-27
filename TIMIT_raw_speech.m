%%%% ASR feature extractor %%%
%%%% ConvRBM learned weights
clear all; close all; clc;
% parpool(4);
load('sorted_MSIS2018Gujarati_Adam_AD0p3_40128_LReLU_21012018.mat');

fs = 16000;
winlen = fs*0.025; % 25 ms window length
winhop = fs*0.010; % 10 ms window hop (overlap)
nfilter=40;
filelist = dir('/home/hardik/Documents/database/gu-in-Measurement/Audios/*.wav');
mkdir('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiBlindTest_40128fbankTEOAM_htk');


parfor index=1:length(filelist)
    fprintf('.\n');
    fname = strcat('/home/hardik/Documents/database/gu-in-Measurement/Audios/',filelist(index).name);
    a1 = filelist(index).name;
    a2 = a1(1:length(a1)-4);
    htkname2 = strcat('/home/hardik/Documents/Features/MSIS2018_feats/Gujarati/GujaratiBlindTest_40128fbankTEOAM_htk/',a2,'.rawf');
    
    [sig,fs] = audioread(fname);
%            sig = pre_processing(sig,fs);
    sig = bsxfun(@minus, sig,mean(sig));
    sig = bsxfun(@rdivide, sig,(std(sig))+0.0001);


  %[featsmean] = AMFM_ConvRBM_feat_extract(sig,fs,W3,winlen,winhop,hbias_vec);

   [feats] = ConvRBM_feat_extract(sig,fs,W3,winlen,winhop);
    %          CC = cep_coeff(feats,40,13,22); %% to get DCT-based features
    %       [feats] = Mel_spec_temp(sig,fs,winlen,winhop,nfilter);
    %       Mel-filterbank features
    %     gspec = gammatonegram(sig,fs,winlen/fs,winhop/fs,nfilter,0,fs/2,0);
    %     gspec = (log(gspec));
    %      ConvRBMCC = dct(feats);
    %      ConvRBMCC = ConvRBMCC(1:13,:);
    % %     gdct = dct(gspec);
    % %     gdct = gdct(1:13,:);
    %      ConvRBMCC = ConvRBMCC';
    %     gdct = gdct';
%     featsmax = featsmax';
    feats = feats';
%         featsmean = bsxfun(@minus, featsmean, mean(featsmean));
%     featsmean = bsxfun(@rdivide, featsmean, std(featsmean));
%     writehtk(htkname1,featsmax,0.01,9);
    writehtk(htkname2,feats,0.01,9); % write it in HTK format to use for ASR
%     writehtk(htkname3,ConvRBMCCmax,0.01,9);
%     writehtk(htkname4,ConvRBMCCmean,0.01,9);
    %           TEO_feats = TEO_feats';
    %          STFT = log(STFT');
    %         writehtk(htkname1,feats,0.01,9);
    %       writehtk(htkname2,ConvRBMCC,0.01,9);
    %      writehtk(htkname3,CC,0.01,9);
    % save(htkname2,'STFT','-mat');
end

