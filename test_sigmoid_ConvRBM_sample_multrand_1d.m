function [poshidact,poshidsample] = test_sigmoid_ConvRBM_sample_multrand_1d(imdata, W, hbias_vec, pars)
% function [poshidactv,poshidprobs2] = ConvRBM_inference_fixconv_1d(imdata, W, hbias_vec, pars)
% cutoff = log(realmin('double'));
% k=1;
ws = size(W,1);
numbases = size(W,2);

poshidexp = zeros(numbases, size(imdata,1)-ws+1);
%   disp(size(poshidexp2));
%  poshidexp2 = zeros(size(imdata,1)-ws+1, numbases);
 poshidact=zeros(numbases, size(imdata,1)-ws+1);
poshidsample = zeros(numbases, size(imdata,1)-ws+1);
% gen_noise = zeros(1, size(imdata,1)-ws+1);
H = reshape(W(end:-1:1, :),[ws,numbases]);

%     poshidexp = conv2_mult(imdata, H, 'valid');
% xact = zeros(1, size(imdata,1)-ws+1);
% tic;
for b=1:numbases
%     gen_noise=zeros(size(poshidexp(b,:)));
    poshidexp(b,:) = conv(imdata, H(:,b), 'valid');
    poshidexp(b,:) = pars.C_sigm/(pars.std_gaussian^2).*(poshidexp(b,:) + hbias_vec(b));
    poshidexp(b,:) = stable_response(poshidexp(b,:));
%     poshidact(b,:) = sigmoid(poshidexp(b,:));
%     poshidsample(b,:) = double(rand(size(poshidact(b,:)))<poshidact(b,:));
end
% toc;
 [poshidact,poshidsample] = softmax(poshidexp);
return
