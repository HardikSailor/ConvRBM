function [H HP Hc HPc] = ConvRBM_sample_multrand_1d(poshidexp, spacing)
%%% looking inside ConvRBM_sample_multrand_1d function %%%%%%%%%%%%
% poshidexp is 2d array
% imagesc(poshidexp);figure;
 cutoff = log(realmin('double'));
 k=1;
 poshidexp(poshidexp*k>-cutoff) = -cutoff/k;
 poshidexp(poshidexp*k<cutoff) = cutoff/k;  %%% this three lines is for preventing NaN values
%%
poshidprobs = poshidexp;
% poshidprobs_mult = zeros(spacing+1, size(poshidprobs,1)*size(poshidprobs,2)*size(poshidprobs,3)/spacing);
% poshidprobs_mult = zeros(spacing+1, size(poshidprobs,1)*size(poshidprobs,2)/spacing);
% poshidprobs_mult(end,:) = 1; %%% changed and solved error
% TODO: replace this with more realistic activation, bases..
for r=1:spacing
    %     temp = poshidprobs(r:spacing:end, :, :);
    x = poshidprobs(:,r:spacing:end);
    temp1 = max(0,x);  %%%% ReLU activation 
    temp2 = max(0,x+((sigmoid(x).^2).*randn(size(x)))); %%%% hidden unit sampling using noisy ReLU
    poshidprobs_mult_act(r,:) = temp1(:);
    poshidprobs_mult_sample(r,:) = temp2(:);
end
P = poshidprobs_mult_act;
S = poshidprobs_mult_sample;
% [S1 P1] = multrand2(poshidprobs_mult');
% S = S1';
% P = P1';
% clear S1 P1
% P = max(0,P11);
% subplot(211);imagesc(P);subplot(212);imagesc(P11);figure;
% convert back to original sized matrix
H = zeros(size(poshidexp));
HP = zeros(size(poshidexp));
for r=1:spacing
    H(:,r:spacing:end) = reshape(S(r,:), [size(H,1), size(H,2)/spacing]);
    HP(:,r:spacing:end) = reshape(P(r,:), [size(H,1), size(H,2)/spacing]);
end

if nargout >2
    Sc = sum(S(:,1:end-1));
    Pc = sum(P(:,1:end-1));
    Hc = reshape(Sc, [size(poshidexp,1), size(poshidexp,2)/spacing]);
    HPc = reshape(Pc, [size(poshidexp,1), size(poshidexp,2)/spacing]);
end

return
