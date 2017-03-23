function [ x_corerction ] = sampson_correction( x_1 , x_2 , H )

delta=[];
for i=1:size(x_1,1)
    Ah=[-(x_1(i,1)*H(2,1)+x_1(i,2)*H(2,2)+H(2,3))+x_2(i,2)*(x_1(i,1)*H(3,1)+x_1(i,2)*H(3,2)+H(3,3));...
        x_1(i,1)*H(1,1)+x_1(i,2)*H(1,2)+H(1,3)-x_2(i,1)*(x_1(i,1)*H(3,1)+x_1(i,2)*H(3,2)+H(3,3))];
    J=[-H(2,1)+x_2(i,2)*H(3,1) -H(2,2)+x_2(i,2)*H(3,2) 0 x_1(i,1)*H(3,1)+x_1(i,2)*H(3,2)+H(3,3);...
        H(1,1)-x_2(i,1)*H(3,1) H(1,2)-x_2(i,1)*H(3,2) -(x_1(i,1)*H(3,1)+x_1(i,2)*H(3,2)+H(3,3)) 0];
    lambda=inv(J*J')*(-Ah);
    tmp=J'*lambda;
    delta=[delta;tmp(1) tmp(2)];
end;
x_corerction = x_1 + delta;

end

