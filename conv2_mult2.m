function y = conv2_mult2(A, b, convopt)
% y = [];
% for i=1:size(A,1)
    y = conv2(b,A', convopt);
% end
return
