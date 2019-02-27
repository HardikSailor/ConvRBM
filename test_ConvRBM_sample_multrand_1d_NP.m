function [poshidact] = test_ConvRBM_sample_multrand_1d_NP(imdata, W, hbias_vec,dropoutrate)
%%%% Function to get hidden unit activation values in negative phase %%%
ws = size(W,1);
numbases = size(W,2);

poshidexp = zeros(numbases, size(imdata,1)+ws-1);

H = reshape(W(end:-1:1, :),[ws,numbases]);


% tic;
for b=1:numbases
    poshidexp(b,:) = conv(imdata, H(:,b), 'full');
    poshidexp(b,:) = poshidexp(b,:) + hbias_vec(b);

    poshidexp(b,:) = stable_response(poshidexp(b,:));
          DropMask(b,:) = rand(size(poshidexp(b,:)))>dropoutrate; %%% K x 1 matrix
      poshidexp(b,:) = poshidexp(b,:).*DropMask(b,:);

end
% toc;
%  poshidact =  max(0,poshidexp); %%% ReLU activation
 poshidact =  max(0,poshidexp)+0.01.*min(0,poshidexp); %%% Leaky ReLU activation

return
