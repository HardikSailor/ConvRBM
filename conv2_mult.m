function y = conv2_mult(a, B, convopt)
y = [];
for i=1:size(B,2)
%     display(size(a));
%     display(size(B(:,:,i)));
    y(i,:) = conv(a, B(:,i), convopt);
%      display(size(y))
end
return
