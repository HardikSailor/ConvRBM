function vishidprod2 = ConvRBM_vishidprod_fixconv_1d(imdata, H, ws)


numbases = size(H,1);
selidx1 = size(H,2):-1:1;
vishidprod2 = zeros(ws,numbases);
for b=1:numbases
    vishidprod2(:,b) = conv2_mult2(imdata, H(b,selidx1), 'valid'); 
end

return
