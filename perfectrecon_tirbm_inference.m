function [poshidactv,poshidprobs2] = perfectrecon_tirbm_inference(imdata, W, hbias_vec, pars)

ws = size(W,1);
numbases = size(W,2);

  poshidexp2 = zeros(numbases, size(imdata,1)-ws+1);
%   disp(size(poshidexp2));
%  poshidexp2 = zeros(size(imdata,1)-ws+1, numbases);
poshidprobs2=zeros(numbases, size(imdata,1)-ws+1);
% poshidsample = zeros(numbases, size(imdata,1)-ws+1);
 if nargout>1
     poshidprobs2 = zeros(size(poshidexp2));
 end

    H = reshape(W(end:-1:1, :),[ws,numbases]);

    poshidexp2 = conv2_mult(imdata, H, 'same');
% xact = zeros(1, size(imdata,1)-ws+1);
for b=1:numbases
    poshidexp2(b,:) = pars.C_sigm/(pars.std_gaussian^2).*(poshidexp2(b,:) + hbias_vec(b));
     if nargout>1
           poshidprobs2(b,:) = 1./(1 + exp(-poshidexp2(b,:)));
%           xact = poshidexp2(b,:)+(poshidprobs2(b,:)).*randn(size(poshidprobs2(b,:)));
%           if xact>0
%           poshidsample(b,:) = xact;
%           else
%             poshidsample(b,:) = 0.01.*xact;
%           end
%             poshidsample(b,:) = max(0,xact);
%            xsample(b,:) = max(0,poshidexp2(b,:));

     end
end
poshidactv =  poshidexp2;
return