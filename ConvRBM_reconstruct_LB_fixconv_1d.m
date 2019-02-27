function negdata = ConvRBM_reconstruct_LB_fixconv_1d(S, W,vbias_vec)

ws = size(W,1);
patch_M = size(S,2);

numbases = size(W,2);

S2 = S;
negdata2 = zeros(1,patch_M-ws+1);

for b = 1:numbases,
    
    negdata2 = negdata2 + conv(S2(b,:), W(:,b), 'valid');
    
end

 negdata = negdata2'+vbias_vec;

return