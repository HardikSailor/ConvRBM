function [H HP Hc HPc] = tirbm_sample_multrand_1d(poshidexp, spacing)
%%% looking inside tirbm_sample_multrand_1d function %%%%%%%%%%%%
% poshidexp is 3d array
% imagesc(poshidexp);figure;
%       poshidexp = max(min(poshidexp,30),-30); % DEBUG: This is for preventing NAN values
%       imagesc(poshidexp);figure;
 cutoff = log(realmin('single'));
 k=1;
 poshidexp(poshidexp*k>-cutoff) = -cutoff/k;
 poshidexp(poshidexp*k<cutoff) = cutoff/k;
%%
poshidprobs = exp(poshidexp);
% poshidprobs_mult = zeros(spacing+1, size(poshidprobs,1)*size(poshidprobs,2)*size(poshidprobs,3)/spacing);
poshidprobs_mult = zeros(spacing+1, size(poshidprobs,1)*size(poshidprobs,2)/spacing);

poshidprobs_mult(end,:) = 1; %%% changed and solved error
% TODO: replace this with more realistic activation, bases..
for r=1:spacing
    %     temp = poshidprobs(r:spacing:end, :, :);
    temp = poshidprobs(:,r:spacing:end);
    
    poshidprobs_mult(r,:) = temp(:);
end
[S1 P1] = multrand2(poshidprobs_mult');
S = S1';
P = P1';
clear S1 P1
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
