function [thetainc,vt] = SSD_optimize(dW,dv,dh,M,J,epsilon)
vbiasinc = (sum(abs(dv)))/J.*sign(dv);
hbiasinc = (sum(abs(dh)))/M.*sign(dh);
[A,lambda,B] = svd(W);
Winc = epsilon*(norm(lambda,1)/(M*J))*A*B';
end