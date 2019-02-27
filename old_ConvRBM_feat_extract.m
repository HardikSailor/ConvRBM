%%%% function to extract features %%%%
function feats = ConvRBM_feat_extract(sig,W4,winlen,winhop)

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
 feats = log(Y+0.001);


end