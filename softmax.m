function [P S] = softmax(poshidexp)
% P is 2-d matrix: 2nd dimension is # of choices
P = exp(poshidexp);
% sumP = row_sum(P); 
sumP = sum(P,1);
P = P./repmat(sumP, [size(P,1),1]);

% cumP = cumsum(P,2);
% rand(size(P));
unifrnd = rand(size(P,1),1);
S = P > repmat(unifrnd,[1,size(P,2)]);
% Sindx = diff(temp,1,2);
% S = zeros(size(P));
% S(:,1) = 1-sum(Sindx,2);
% S(:,2:end) = Sindx;
 S = double(S);
% size(P)
% size(S)
end