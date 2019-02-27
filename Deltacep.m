function [d]=Deltacep(x,N)
%% x is stored such that every column represents a feture vector 

[n1,n2]=size(x); % n1= dimension , n2 = number of frames 

%% Do frame repeatation part 
x1=[];
for k=1:N
x1=[x1,x(:,1)];
end
x1=[x1,x];
for k=1:N
x1=[x1,x(:,end)];
end
% ,x,x(:,end-N+1:end)];
for j=1:n1
for k=1+N:size(x1,2)-N % for every frame 
%     sm1=zeros(n1,1); % zero vector of feature dim 
sm1=0;
    sm2=0; % denome
    for k1=1:N
        sm1=sm1+(k1*(x1(j,k+k1)-x1(j,k-k1)));
        sm2=sm2+k1^2;
    end
    sm2=2*sm2; 
    y(j,k)=sm1/sm2;
end
end
%% Do frame removing
d=y(:,N+1:size(x1,2)-N);
