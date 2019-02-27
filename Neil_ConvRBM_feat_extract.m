%%%% function to extract features %%%%
function [feats] = Neil_ConvRBM_feat_extract(sig,fs,W3,winlen,winhop)
     sig = pre_processing(sig,fs);
    sig = bsxfun(@minus, sig,mean(sig));
    sig = bsxfun(@rdivide, sig,(std(sig))+0.0001);
 Y = [];
[l,b]=size(W3);
%   W4 = bsxfun(@minus,W4,mean(W4));
%   W4 = bsxfun(@rdivide,W4,max(W4));
        W3 = reshape(W3(end:-1:1, :),[l,b]);
% res = zeros(b,length(sig));
% [b1,a1] = butter(1,0.125);

for j=1:b
%          filter_td =W4(:,j)./norm(W4(:,j));
     filter_td =W3(:,j);
    y = (conv(sig,filter_td,'same'));%+hbias_vec(j);

%   y = bsxfun(@minus, y, mean(y));
%                    y = max(0,y);
%                  x_teager  = cal_freqweighted_energy(y,fs,'abs_teager');
%      y = bsxfun(@rdivide, y, std(y));
          y = abs(y);
%                   ylp =  filter(b1,a1,y);

%     y = hilbert(y);
% instfreq = fs/(2*pi)*diff(unwrap(angle(y)));
% %     res(j,:) = y;
%      out = enframe(instfreq,winlen,winhop);
         out = enframe(y,hamming(winlen),winhop);
% 
%                      out3 = max(out,[],2);
% %  out = out';
                            out3 = mean(out,2);
% %       out2 = sqrt(mean(out.^2,2));
% 
     %   Y1(j,:) = out2;
                 Y2(j,:) = out3;

%        Y(j,:) = y;
end
%  Y = Y';
%  display(size(Y));
% WIF = IF_AMweighted(Y,fs);
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
        % feats1 = log(Y1+0.0001);
                     feats = log(Y2+0.0001);

%            feats = Y2.^(1/15);
%      feats = feats';
%     feats = bsxfun(@minus,feats,mean(feats));
%     feats = bsxfun(@rdivide,feats,std(feats));
%     feats = feats';

%    feats = Y;
end