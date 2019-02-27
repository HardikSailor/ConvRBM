function ReLU_RAW_CRBM_hardik_main(ws, num_bases, spacing, pbias, pbias_lb, pbias_lambda, epsilon, l2reg, epsdecay)
% clear all; close all; clc;
%%% NOTE:
%%% 1. visible bias is NOT used. 
% addpath ../voicebox
% addpath F:\Deep_Learning\Tools\crbm_audio\voicebox;
dataname = 'TIMIT';
batch_size = 1;
filelist = dir('H:\TIDIGITS_16khz\');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the ConvRBM parameters
sigma_start = 0.2;
sigma_stop = 0.2;
%sparseCost=0.5;
%sparseTarget=0.3;
%dsh=0;
%meanActiv=0;
%lambda = 0.9;
CD_mode = 'mf';
bias_mode = 'none'; %%% sparsity regularization

% Etc parameters
K_CD = 1; %%% Single step constrastive divergence (CD)

% Initialization
W = [];
vbias_vec = [];
hbias_vec = [];
pars = [];

C_sigm = 1;

num_trials = 50;  %%% No. of training epochs

% Initialize variables
if ~exist('pars', 'var') || isempty(pars)
    pars=[];
end

if ~isfield(pars, 'ws'), pars.ws = ws; end
if ~isfield(pars, 'num_bases'), pars.num_bases = num_bases; end
if ~isfield(pars, 'spacing'), pars.spacing = spacing; end

if ~isfield(pars, 'pbias'), pars.pbias = pbias; end
if ~isfield(pars, 'pbias_lb'), pars.pbias_lb = pbias_lb; end
if ~isfield(pars, 'pbias_lambda'), pars.pbias_lambda = pbias_lambda; end
if ~isfield(pars, 'bias_mode'), pars.bias_mode = bias_mode; end

if ~isfield(pars, 'std_gaussian'), pars.std_gaussian = 1; end
if ~isfield(pars, 'sigma_start'), pars.sigma_start = sigma_start; end
if ~isfield(pars, 'sigma_stop'), pars.sigma_stop = sigma_stop; end

if ~isfield(pars, 'K_CD'), pars.K_CD = K_CD; end
if ~isfield(pars, 'CD_mode'), pars.CD_mode = CD_mode; end
if ~isfield(pars, 'C_sigm'), pars.C_sigm = C_sigm; end

if ~isfield(pars, 'num_trials'), pars.num_trials = num_trials; end
if ~isfield(pars, 'epsilon'), pars.epsilon = epsilon; end

disp(pars)

if ~exist('W', 'var') || isempty(W)
    %     W = 0.001*randn(pars.ws, numchannels, pars.num_bases);
    W = 0.01*randn(pars.ws, pars.num_bases);
    %          W = -1+2*rand(pars.ws, pars.num_bases);
    %       W = sqrt(2/pars.ws*pars.spacing).*randn(pars.ws, pars.num_bases);
end

if ~exist('vbias_vec', 'var') || isempty(vbias_vec)
    %      vbias_vec = 1*ones(pars.num_bases,1);
    vbias_vec = zeros(1,1);
end

if ~exist('hbias_vec', 'var') || isempty(hbias_vec)
    %           hbias_vec = -0.01*ones(pars.num_bases,1);
    hbias_vec = zeros(pars.num_bases,1);
    %      hbias_vec = -1*rand(pars.num_bases,1);
end

fname_prefix = sprintf('TIDIGITS_raw_speech_128_16_40_LB_%s_V1b_w%d_b%02d_pc%d_sigpc%g_p%g_pl%g_plambda%g_sp%d_%s_eps%g_epsdecay%g_l2reg%g_bs%02d_%s', dataname, ws, num_bases, pbias, pbias_lb, pbias_lambda, spacing, CD_mode, epsilon, epsdecay, l2reg, batch_size, datestr(now, 30));
fname_save = sprintf('%s', fname_prefix);
fname_mat  = sprintf('%s.mat', fname_save);
fname_out = fname_mat;
mkdir(fileparts(fname_save));
fname_out

%%%% Setting momentum
initialmomentum  = 0.5;
finalmomentum    = 0.9;

error_history = [];
sparsity_history = [];

dWnorm_history_CD= [];
dWnorm_history_l2= [];

Winc=0;
vbiasinc=0;
hbiasinc=0;
filelist_temp = filelist; %%% take full filelist or less
for t=1:num_trials
    veclist = randperm(numel(1:floor(length(filelist_temp))));
    if(t>10)  %%%% fixed learning rate (LR) for first 10 epochs
        epsilon = pars.epsilon/(1+epsdecay*t);
        %     epsilon_h = epsilon_h/(1+epsdecay*t); %%% ON if different LR for
        %     hidden units
    end
    
    tic;
    ferr_current_iter = [];
    sparsity_curr_iter = [];
    tt1 = 0;
    for j=veclist;
        tt1 = tt1+1;
        %         disp(j)
        %          rand_j =
        %          floor((length(filelist))*rand(1));randi(numel(filelist));
        %         ON for Randomly selecting files
        %          
        fprintf('.');
        %            fprintf('Processing %s\n', filelist(rand_j).name);
        fname=strcat('H:\TIDIGITS_16khz\',filelist(j).name);
        [sig,fs]=wavread(fname);
        %%% Zero man and unit variance pre-processing
        sig = bsxfun(@minus, sig,mean(sig));
        sig = bsxfun(@rdivide, sig,(std(sig)));

        imdata = hardik_trim_audio_for_spacing_fixconv(sig, ws, spacing); %%% Trim speech file for valid length convolution
        imdatatr = imdata;
        %%% Momentum is fixed 0.5 for first five training epochs and then
        %%% set to 0.9
        if t<5,
            momentum = initialmomentum;
        else
            momentum = finalmomentum;
        end
        [ferr, dW, dh, dv, poshidprobs, poshidstates, negdata, stat]= fobj_tirbm_CD_LB_sparse_audio(imdatatr, W, hbias_vec, vbias_vec, pars, CD_mode, bias_mode, spacing, l2reg,momentum);
        ferr_current_iter = [ferr_current_iter, ferr];
        sparsity_curr_iter = [sparsity_curr_iter, mean(poshidprobs(:))];
        
        dWnorm_history_CD(j,t) = stat.dWnorm_CD;
        dWnorm_history_l2(j,t) = stat.dWnorm_l2;
        
        % update parameters
        Winc = momentum*Winc + epsilon*dW;
        W = W + Winc;
        %%% use below code for RBM debuging for weights
        % if(mod(tt1,5)==0)
        %  subplot(221);plot(Winc(:,1)); subplot(222);plot(Winc(:,20));subplot(223);plot(Winc(:,30));subplot(224);plot(Winc(:,10));figure;
        %   subplot(221);plot(W(:,1)); subplot(222);plot(W(:,20));subplot(223);plot(W(:,30));subplot(224);plot(W(:,10));figure;
        %  end
        vbiasinc = momentum*vbiasinc + epsilon*dv;
        vbias_vec = vbias_vec + vbiasinc;
        
        hbiasinc = momentum*hbiasinc + epsilon*dh;
        hbias_vec = hbias_vec + hbiasinc;
        % spectrogram(negdata,400,160);figure;
        % subplot(221);plotit(W);subplot(222);plotit(hbias_vec);subplot(223);plotit(Winc);subplot(224);plotit(hbiasinc);figure;
    end
    mean_err = mean(ferr_current_iter);
    mean_sparsity = mean(sparsity_curr_iter);
    %%% below code NOT required for NOW !
    %     if (pars.std_gaussian < sigma_stop) % stop decaying after some point
    %         pars.std_gaussian = pars.std_gaussian*0.95;
    %     end
    error_history(t) = mean(ferr_current_iter);
    sparsity_history(t) = mean(sparsity_curr_iter);
    
    toc
    %%% ON below code to save files at each epochs
    %     if(mod(t,10)==0)
    % fname_prefix = sprintf('500_raw_speech_128_32_80_LB_%s_V1b_w%d_b%02d_pc%d_sigpc%g_p%g_pl%g_plambda%g_sp%d_%s_eps%g_epsdecay%g_l2reg%g_bs%02d_%s', dataname, ws, num_bases, pbias, pbias_lb, pbias_lambda, spacing, CD_mode, epsilon, epsdecay, l2reg, batch_size, datestr(now, 30));
    % fname_save = sprintf('%s', fname_prefix);
    % fname_mat  = sprintf('%s.mat', fname_save);
    % fname_out = fname_mat;
    % mkdir(fileparts(fname_save));
    % fname_out
    
    fprintf('epoch %d error = %g \tsparsity_hid = %g\n', t, mean(ferr_current_iter), mean(sparsity_curr_iter));
    save(fname_mat, 'W', 'pars', 't', 'vbias_vec', 'hbias_vec', 'error_history', 'sparsity_history', 'dWnorm_history_CD', 'dWnorm_history_l2');
    %     end
    %         save(fname_mat, 'W', 'pars', 't', 'vbias_vec', 'hbias_vec', 'error_history', 'sparsity_history', 'E', 'S', 'sigmaPC', 'Ewhiten', 'Eunwhiten', 'dWnorm_history_CD', 'dWnorm_history_l2');
      %%% ON below code to save files at each epochs
    %     if 0 %mod(t, 10)==0
    %         fname_mat_timestamp  = sprintf('%s_%04dEPOCHS.mat', fname_save, t);
    %         save(fname_mat_timestamp, 'W', 'pars', 't', 'vbias_vec', 'hbias_vec', 'error_history', 'sparsity_history', 'E', 'S', 'sigmaPC', 'Ewhiten', 'Eunwhiten', 'dWnorm_history_CD', 'dWnorm_history_l2');
    %     end
    
    disp(sprintf('results saved as %s\n', fname_mat));
    %          subplot(221);plotit(W);subplot(221);plotit(hbias_vec);subplot(223);plotit(Winc);subplot(224);plotit(hbiasinc);
    %  plot(negdata(:));
    %  subplot(221);plot(W(:,1)); subplot(222);plot(W(:,20));subplot(223);plot(W(:,30));subplot(224);plot(W(:,10));figure;
    %     imagesc(poshidprobs);figure;
end

return


function [ferr dW_total dh_total dv_total poshidprobs poshidstates negdata stat] = ...
    fobj_tirbm_CD_LB_sparse_audio(imdata, W, hbias_vec, vbias_vec, pars, CD_mode, bias_mode, spacing, l2reg,momentum)

ws = size(W,1);

poshidexp= ConvRBM_inference_fixconv_1d(imdata, W, hbias_vec, pars);
%   imagesc(poshidexp);figure;
[poshidstates poshidprobs] = ConvRBM_sample_multrand_1d(poshidexp, spacing);
%  if strcmp(CD_mode, 'mf'), poshidstates = poshidprobs; end
%     imagesc(poshidstates);colormap('jet');figure;
posprods = tirbm_vishidprod_fixconv_1d(imdata, poshidprobs, ws);
poshidact = squeeze(sum(poshidprobs,2));
posvisact = squeeze(sum(imdata,1));

% neghidprobs = poshidprobs;
% poshidstates = poshidprobs ;
for j=1:pars.K_CD
    negdata = tirbm_reconstruct_LB_fixconv_1d(poshidstates, W, pars);
    %      negdata  = bsxfun(@minus, negdata ,mean(negdata ));
    %      negdata  = bsxfun(@rdivide, negdata ,(std(negdata )));
    [neghidexp,neghidprobs2] = tirbm_inference_fixconv_1d(negdata, W, hbias_vec, pars);
    [neghidstates neghidprobs] = my_tirbm_sample_multrand_1d(neghidexp, spacing);
    %     if strcmp(CD_mode, 'mf'), neghidstates = neghidprobs; end
end
negprods = tirbm_vishidprod_fixconv_1d(negdata, neghidprobs, ws);
neghidact = squeeze(sum(neghidprobs,2));
negvisact = squeeze(sum(negdata,1));
%     plot(negdata(:));figure;
%    spectrogram(negdata(:),400,160);figure;
ferr = mean( (imdata(:)-negdata(:)).^2 );

if strcmp(bias_mode, 'none')
    dhbias = 0;
    dvbias = 0;
    dW = 0;
elseif strcmp(bias_mode, 'simple')
    %     plot(squeeze(mean(poshidprobs,2)));figure;
    dhbias = squeeze(mean(poshidprobs,2)) - pars.pbias;
    %       hist(mean(poshidprobs,2));figure;
    dvbias = 0;
    dW = 0;
else
    error('wrong adjust_bias mode!');
end
% meanActiv = lembda*meanActiv+(1-lembda)*(mean(mean(poshidprobs)));
% dsh = meanActiv-sparseTarget;
% numcases1 = size(poshidprobs,1)*size(poshidprobs,2);
% numcases1 = 1;
% numcases2=1;
numcases1 = size(poshidprobs,2);
%
numcases2 = size(imdata,1)*size(imdata,2);

dW_total1 = (posprods-negprods)/numcases1;
dW_total2 = - l2reg*W;
% dW_total2 = - 0.1*sign(W);
%  dW_total3 = - pars.pbias_lambda*dW;
%  dW_total4 = - pars.pbias_lambda*dhbias;
% dW_total3 = (1-momentum)*dW;
% disp(dW_total3);
dW_total = dW_total1 + dW_total2;

stat = [];
stat.dWnorm_CD = norm(vec(dW_total1));
stat.dWnorm_l2 = norm(vec(dW_total2));


dh_total = (poshidact-neghidact)/numcases1 - pars.pbias_lambda*dhbias;
% dh_total = (poshidact-neghidact)/numcases1 - - sparseCost * dsh;

dv_total = ((posvisact-negvisact)/numcases2);

% fprintf('||W||=%g, ||dWprod|| = %g, ||dWl2|| = %g, ||dWsparse|| = %g\n', sqrt(sum(W(:).^2)), sqrt(sum(dW_total1(:).^2)), sqrt(sum(dW_total2(:).^2)), sqrt(sum(dW_total3(:).^2)));

return


% function im2 = trim_audio_for_spacing_fixconv(im2, ws, spacing)
% if mod(size(im2,2)-ws+1, spacing)~=0
%     n = mod(size(im2,2)-ws+1, spacing);
%     im2(:, 1:floor(n/2), :) = [];
%     im2(:, end-ceil(n/2)+1:end, :) = [];
% end
% return
%
%
%
% function flist = get_speech_filenames2()
%
% flist = [];
% fid = fopen('F:\Deep_Learning\Tools\crbm_audio\code\flist_train.txt', 'r');
% while 1
%     tline = fgetl(fid);
%     if ~ischar(tline), break, end
%     % disp(tline);
%     flist{length(flist)+1} = tline;
% end
% fclose(fid);
%
% return



