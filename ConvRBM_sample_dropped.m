function [poshidact,poshidsample] = ConvRBM_sample_dropped(imdata, W, hbias_vec, pars)
% function [poshidactv,poshidprobs2] = ConvRBM_inference_fixconv_1d(imdata, W, hbias_vec, pars)
% cutoff = log(realmin('double'));
% k=1;
ws = size(W,1);
numbases = size(W,2);
dropoutrate = 0.7;
poshidexp = zeros(numbases, size(imdata,1)-ws+1);
%   disp(size(poshidexp2));
%  poshidexp2 = zeros(size(imdata,1)-ws+1, numbases);
% poshidact=zeros(numbases, size(imdata,1)-ws+1);
poshidsample = zeros(numbases, size(imdata,1)-ws+1);
% gen_noise = zeros(1, size(imdata,1)-ws+1);
H = reshape(W(end:-1:1, :),[ws,numbases]);
 DropMask = rand(size(poshidexp,2),1)>dropoutrate;
%     poshidexp = conv2_mult(imdata, H, 'valid');
% xact = zeros(1, size(imdata,1)-ws+1);
% tic;
for b=1:numbases
%     gen_noise=zeros(size(poshidexp(b,:)));
    poshidexp(b,:) = conv(imdata, H(:,b), 'valid');
    poshidexp(b,:) = pars.C_sigm/(pars.std_gaussian^2).*(poshidexp(b,:) + hbias_vec(b));
    poshidexp(b,:) = stable_response(poshidexp(b,:));
    if nargout>1
%     poshidexp(poshidexp*k>-cutoff) = -cutoff/k;
%     poshidexp(poshidexp*k<cutoff) = cutoff/k;  %%% this lines is for preventing NaN values
    gen_noise = (sigmoid(poshidexp(b,:)).^2).*randn(size(poshidexp(b,:))); %%% zero mean sigmoid(x)as variance
    poshidsample(b,:) = DropMask.*max(0,poshidexp(b,:)+gen_noise);
    poshidact(b,:) = DropMask.*max(0,poshidexp(b,:)); 
%     clear gen_noise;
    end
end
% toc;
% poshidact =  max(0,poshidexp);
return
