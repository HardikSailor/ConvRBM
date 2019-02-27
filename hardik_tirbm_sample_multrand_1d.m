function [H] = hardik_tirbm_sample_multrand_1d(poshidexp, spacing)
%%% looking inside tirbm_sample_multrand_1d function %%%%%%%%%%%%
% poshidexp is 3d array
% imagesc(poshidexp);figure;
    poshidexp = max(min(poshidexp,30),-30); % DEBUG: This is for preventing NAN values
%   imagesc(poshidexp);figure;  
% cutoff = log(realmin('double'));
% k=1;
% poshidexp(poshidexp*k>-cutoff) = -cutoff/k;
% poshidexp(poshidexp*k<cutoff) = cutoff/k;
%%
 poshidprobs2 = exp(poshidexp);
% poshidprobs_mult = zeros(spacing+1, size(poshidprobs,1)*size(poshidprobs,2)*size(poshidprobs,3)/spacing);
poshidprobs_mult = zeros(spacing, size(poshidprobs2,1)*size(poshidprobs2,2)/spacing);

% poshidprobs_mult(end,:) = 1; %%% changed and solved error
% TODO: replace this with more realistic activation, bases..
for r=1:spacing
    %     temp = poshidprobs(r:spacing:end, :, :);
    temp = poshidprobs2(:,r:spacing:end);
    
    poshidprobs_mult(r,:) = temp(:);
end
S1 = hardik_multrand2(poshidprobs_mult');
S = S1';
clear S1
% P = max(0,P11);
% subplot(211);imagesc(P);subplot(212);imagesc(P11);figure;
% convert back to original sized matrix
H = zeros(size(poshidprobs2));
% HP = zeros(size(poshidprobs2));
for r=1:spacing
    H(:,r:spacing:end) = reshape(S(r,:), [size(H,1), size(H,2)/spacing]);
%     HP(:,r:spacing:end) = reshape(P(r,:), [size(H,1), size(H,2)/spacing]);
end

% if nargout >2
%     Sc = sum(S(:,1:end-1));
%     Pc = sum(P(:,1:end-1));
%     Hc = reshape(Sc, [size(poshidprobs2,1), size(poshidprobs2,2)/spacing]);
%     HPc = reshape(Pc, [size(poshidprobs2,1), size(poshidprobs2,2)/spacing]);
% end

return
