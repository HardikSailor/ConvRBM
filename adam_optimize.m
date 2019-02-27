function [thetainc,mt,vt] = adam_optimize(dtheta,mt,vt,b1,b2,epsilon,t)
mt  = b1*mt + (1-b1)*dtheta;
vt  = b2*vt + (1-b2)*(dtheta.^2);
mcapt = mt/(1-(b1.^t));
vcapt = vt/(1-(b2.^t));
thetainc = epsilon*(mcapt./(sqrt(vcapt))+10.^-8);
end