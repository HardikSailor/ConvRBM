%%%%% CRBM feature extraction %%%%%%%%%

clear all; clc; close all;

% ws = 6;
% spacing = 3;
load('raw_speech_50_10_40_LB_TIMIT_V1b_w50_b80_pc1.000000e-002_sigpc0.01_p0.9_pl10_plambda109_sp102_1.000000e-001_eps0.01_epsdecay0.01_l2reg1_bs50_0150503T140657.mat');

filelist = dir('E:\TIMIT_features_new\mfs\*.wav.mfs');
CD_mode = 'mf';
bias_mode = 'simple';
sigmaPC = 3;
ws = 50;
spacing = 10;
num_bases = 80;
pbias=0.01;
pbias_lb=0.01;
pbias_lambda=0.1;
epsilon=0.1;
l2reg=0.01;
epsdecay=0.01;
for index=1:1
    tic;
    fprintf('Processing %s\n', filelist(index).name);
    fname=strcat('E:\TIMIT_features_new\mfs\',filelist(index).name);
    name1=strcat('E:\TIMIT_features_new\CRBM2_MFS\',filelist(index).name);
    htkname1=strcat('E:\TIMIT_features_new\CRBM2_MFS_htk\',filelist(index).name);
    %     name2=strcat('E:\TIMIT_features_new\MFS_W\',filelist(index).name);
    %     htkname2=strcat('E:\TIMIT_features_new\MFS_W_htk\',filelist(index).name);
    data = load(fname);
    imdata = Ewhiten*log(data);
    imdata = my_trim_audio_for_spacing_fixconv(imdata, ws, spacing);
    
    %         MFS = log(data);
    imdatatr = imdata';
    %     mfs_feats = MFS';
    imdatatr = reshape(imdatatr, [size(imdatatr,1), 1, size(imdatatr,2)]);
    
    recon_data  = tirbm_inference_fixconv_1d(imdatatr, W, hbias_vec, pars);
    %     [ferr dW_total dh_total dv_total poshidprobs poshidstates negdata stat] = my_tirbm_CD_LB_sparse_audio(imdatatr, W, hbias_vec, vbias_vec, pars, 'mf', bias_mode, spacing, l2reg);
    
    CRBM_feats=reshape(recon_data,[size( recon_data ,1) size( recon_data ,3)]);
    %%%%% making same length %%%
    [h1,h2]=size(CRBM_feats);
    [o1,o2]=size(data');
    a_new =CRBM_feats;
    if h1<o1
        extra = o1-h1;
        for i=1:extra
            a_new = [a_new;CRBM_feats(end,:)];
        end
    end
%     size(a_new)
%     size(data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     CRBM_feats = CRBM_feats';
    
    save('-ascii',name1,'a_new');
    writehtk(htkname1,a_new,0.01,9);
    %     save('-ascii',name2,'mfs_feats');
    %     writehtk(htkname2,mfs_white,0.01,9);
    clear data;
    clear imdatatr;
    clear imdata;
    clear a_new;
    %     clear recon_data;
    toc;
end
