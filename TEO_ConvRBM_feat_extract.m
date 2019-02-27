%%%% function to extract features %%%%
function [TEO_feats] = TEO_ConvRBM_feat_extract(sig,fs,W4,winlen,winhop,hbias_vec)
% Y1 = [];
 Y2 = [];
[l,b]=size(W4);
%   W4 = bsxfun(@minus,W4,mean(W4));
%   W4 = bsxfun(@rdivide,W4,max(W4));
W4 = reshape(W4(end:-1:1, :),[l,b]);
% res = zeros(b,length(sig));
[b1,a1] = butter(1,0.125);

for j=1:b
    %          filter_td =W4(:,j)./norm(W4(:,j));
    filter_td =W4(:,j);
    y = (conv(sig,filter_td,'same'));
%     y = abs(y);
%     ymax = max(0,y)+0.01.*min(0,y);
    ylp =  filter(b1,a1,y);
%     ybank =  filter(b1,a1,ymax);

    y3 = abs(teager(ylp));
    [AM, FM]=AMFM_all_update(y3,fs);
    AM = abs(AM);
%     out1 = enframe(ylp,hamming(winlen),winhop);
    out2 = enframe(AM,hamming(winlen),winhop);
    
    %                out2 = max(out,[],2);
    %  out = out';
%     st1 = mean(out1,2);
        st2 = mean(out2,2);

%            st1 = sqrt(mean(out1.^2,2));
    
%           Y1(j,:) = st1;
    Y2(j,:) = st2;
end
%   Y = Y';
%  Y = bsxfun(@minus,Y,mean(Y));
%   Y = bsxfun(@rdivide,Y,std(Y));
%   Y = Y';
%  Y = max(Y,0);
%    res = max(res,0);
% res = res';
%  res = bsxfun(@minus, res, mean(res));
%  res = bsxfun(@rdivide, res, std(res));
% res = (res);
%  Y = max_pooled((res),winlen,winhop);
% feats = log(Y1+0.001);
TEO_feats = log(Y2+0.0001);

%         feats = Y.^(1/15);
%     feats = feats';
%    feats = bsxfun(@minus,feats,mean(feats));
%   feats = bsxfun(@rdivide,feats,std(feats));
%    feats = feats';

%  feats = Y;
end