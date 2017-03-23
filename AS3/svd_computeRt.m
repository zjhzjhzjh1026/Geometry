function [ Rt , flag ] = svd_computeRt( A , B )

miu_x=mean(A,2);
miu_y=mean(B,2);
sigma_x=sum(var(A,1,2));
sigma_y=sum(var(B,1,2));
Sigma_xy=zeros(3,3);
for m=1:3
    Sigma_xy =Sigma_xy + (B(:,m)-miu_y)*(A(:,m)-miu_x)';
end;
Sigma_xy=Sigma_xy/3;
[U,D,V] = svd(Sigma_xy);
S=eye(3);
if(det(U)*det(V)==-1)
   S(3,3)=-1;
end;
R=U*S*V';
c=trace(D*S)/sigma_x;
t=miu_y-c*R*miu_x;
Rt=[R,t];
flag=0;
if(rank(Sigma_xy)<(size(A,1)-1))
    flag=1;
end;

end

