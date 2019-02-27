%%%% function to extract features %%%%
function [feats] = ConvRBM_feat_extract(sig,fs,W3,winlen,winhop)
Y = [];
[l,b]=size(W3);

W3 = reshape(W3(end:-1:1, :),[l,b]);

for j=1:b
    filter_td =W3(:,j);
    y = (conv(sig,filter_td,'same'));
    
    %   y = max(0,y); %%% ReLU nonlinearity at test time
    
    y = abs(y); %%% absolute nonlinearity at test time
    
    out = enframe(y,hamming(winlen),winhop);
    %
    %   out3 = max(out,[],2); %%% max pooling
    % %  out = out';
    out3 = mean(out,2); %%% average pooling
    % % out2 = sqrt(mean(out.^2,2)); %%% L2 pooling
    %
    %   Y1(j,:) = out2;
    Y2(j,:) = out3;
    
    %        Y(j,:) = y;
end

% feats1 = log(Y1+0.0001);
feats = log(Y2+0.0001);
%   feats = Y2.^(1/15); %%% power-law nonlinearity
%      feats = feats';
%     feats = bsxfun(@minus,feats,mean(feats));
%     feats = bsxfun(@rdivide,feats,std(feats));
%     feats = feats';
%    feats = Y;
end