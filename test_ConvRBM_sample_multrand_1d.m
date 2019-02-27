function [poshidact,poshidsample,DropMask] = test_ConvRBM_sample_multrand_1d(imdata, W, hbias_vec,dropoutrate)

ws = size(W,1);
numbases = size(W,2);

poshidexp = zeros(numbases, size(imdata,1)+ws-1);

poshidsample = zeros(numbases, size(imdata,1)+ws-1);
DropMask = zeros(numbases, size(imdata,1)+ws-1);

H = reshape(W(end:-1:1, :),[ws,numbases]);


for b=1:numbases
    poshidexp(b,:) = conv(imdata, H(:,b), 'full');
    poshidexp(b,:) = poshidexp(b,:) + hbias_vec(b);
    
    poshidexp(b,:) = stable_response(poshidexp(b,:));
    DropMask(b,:) = rand(size(poshidexp(b,:)))>dropoutrate; 
    poshidexp(b,:) = poshidexp(b,:).*DropMask(b,:);
    if nargout>1
        gen_noise = (sigmoid(poshidexp(b,:)).^2).*randn(size(poshidexp(b,:))); %%% zero mean sigmoid(x)as variance
        %      poshidsample(b,:) = max(0,poshidexp(b,:)+gen_noise); %%% NReLU
        poshidsample(b,:) = max(0,poshidexp(b,:)+gen_noise)+0.01.*min(0,poshidexp(b,:)+gen_noise); %% Leaky  NReLU
        
    end
end
% toc;
%  poshidact =  max(0,poshidexp); %%% ReLU activation
poshidact =  max(0,poshidexp)+0.01.*min(0,poshidexp); %%% Leaky ReLU activation

return
