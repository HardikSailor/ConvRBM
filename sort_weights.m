function [Wsorted,I1] = sort_weights(W,num)

[a,b] = size(W);

if num == a
    a = num;
    W = W';
    nfft = 2*b+1;
else
    b = num;
    nfft = 2*a+1;

end
wfft = abs(fft(W,nfft));
wfft=wfft(1:(floor(nfft/2)),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 locations = zeros(1,num);
for i=1:num
     filter_fd = wfft(:,i)./max(wfft(:,i));
     fp = find(filter_fd==1);
     locations(i) = fp;
end
[s1,I1] = sort(locations);

Wsorted = W(:,I1);


end