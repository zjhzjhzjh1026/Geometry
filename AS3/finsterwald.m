function [ p_out ] = finsterwald( label , X_in , x_normalized )
num1=label(1);
num2=label(2);
num3=label(3);

p1=X_in(num1,:);
p2=X_in(num2,:);
p3=X_in(num3,:);

u1=x_normalized(num1,1);
v1=x_normalized(num1,2);
u2=x_normalized(num2,1);
v2=x_normalized(num2,2);
u3=x_normalized(num3,1);
v3=x_normalized(num3,2);

%% Finsterwalder starter
a=sqrt(sum((p2-p3).^2));
b=sqrt(sum((p3-p1).^2));
c=sqrt(sum((p2-p1).^2));

f=1;
j1=[u1,v1,f]'/sqrt(u1^2+v1^2+f^2);
j2=[u2,v2,f]'/sqrt(u2^2+v2^2+f^2);
j3=[u3,v3,f]'/sqrt(u3^2+v3^2+f^2);

cos_a=j2'*j3;
cos_b=j1'*j3;
cos_c=j1'*j2;

%% root lambda
sin_a2=1-cos_a^2;
sin_b2=1-cos_b^2;
sin_c2=1-cos_c^2;
a2=a^2;
b2=b^2;
c2=c^2;

G=c2*(c2*sin_b2-b2*sin_c2);
H=b2*(b2-a2)*sin_c2+c2*(c2+2*a2)*sin_b2+2*b2*c2*(-1+cos_a*cos_b*cos_c);
I=b2*(b2-c2)*sin_a2+a2*(a2+2*c2)*sin_b2+2*a2*b2*(-1+cos_a*cos_b*cos_c);
J=a2*(a2*sin_b2-b2*sin_a2);
fac=[G H I J];
lambda=roots(fac);
for i = 1:size(lambda,1)        %%pick any real lambda
    if(isreal(lambda(i)))
        lambda_0=lambda(i);
        break;
    end;
end;

%% A-F
A=1+lambda_0;
B=-cos_a;
C=(b2-a2)/b2-lambda_0*c2/b2;
D=-lambda_0*cos_c;
E=(a2/b2+lambda_0*c2/b2)*cos_b;
F=-a2/b2+lambda_0*((b2-c2)/b2);

%% u
p=sqrt(B^2-A*C);
q=sign(B*E-C*D)*sqrt(E^2-C*F);

m1=(-B+p)/C;
n1=(-(E-q))/C;
A1=b2-m1*m1*c2;
B1=c2*(cos_b-n1)*m1-b2*cos_c;
C1=-c2*n1^2+2*c2*n1*cos_b+b2-c2;
u1=-sign(B1)/A1*(abs(B1)+sqrt(B1^2-A1*C1));
u2=C1/A1/u1;
p_out=[];
if(isreal(u1))
    pp1 = compute_pi( u1,m1,n1,a2,cos_a,j1,j2,j3 );
    pp2 = compute_pi( u2,m1,n1,a2,cos_a,j1,j2,j3 );
    p_out= cat(3,p_out,pp1);
    p_out= cat(3,p_out,pp2);
end;

m2=(-B-p)/C;
n2=(-(E+q))/C;
A2=b2-m2*m2*c2;
B2=c2*(cos_b-n2)*m2-b2*cos_c;
C2=-c2*n2^2+2*c2*n2*cos_b+b2-c2;
u3=-sign(B2)/A2*(abs(B2)+sqrt(B2^2-A2*C2));
u4=C2/A2/u3;                    %% pick real
if(isreal(u3))
    pp1 = compute_pi( u3,m2,n2,a2,cos_a,j1,j2,j3 );
    pp2 = compute_pi( u4,m2,n2,a2,cos_a,j1,j2,j3 );
    p_out= cat(3,p_out,pp1);
    p_out= cat(3,p_out,pp2);
end;

end

