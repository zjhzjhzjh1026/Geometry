function [ F ] = seven_points_algo( label , x1_ho , x2_ho  , poly_func)

x_seven1 = zeros(7,3);
x_seven2 = zeros(7,3);
a = zeros(7,9);
for i=1:7
    x_seven1(i,:) = x1_ho(label(i),:);
    x_seven2(i,:) = x2_ho(label(i),:);
end;
miux1=mean(x_seven1,1);
varx1=var(x_seven1);
sx1=sqrt(2/sum(varx1));

miux2=mean(x_seven2,1);
varx2=var(x_seven2);
sx2=sqrt(2/sum(varx2));

T1=[sx1,0,-miux1(1)*sx1;0,sx1,-miux1(2)*sx1;0,0,1];
T2=[sx2,0,-miux2(1)*sx2;0,sx2,-miux2(2)*sx2;0,0,1];

x1_nor=x_seven1*T1';      %%left multiply ->right multiply
x2_nor=x_seven2*T2';

for i = 1:7
    a(i,:) = kron(x2_nor(i,:),x1_nor(i,:));
end;

[~,~,V] = svd(a);
F1 = V(:,9);
F2 = V(:,8);

root_poly = poly_func(F1(1),F1(2),F1(3),F1(4),F1(5),F1(6),F1(7),F1(8),...
    F1(9),F2(1),F2(2),F2(3),F2(4),F2(5),F2(6),F2(7),F2(8),F2(9));
root_poly_reverse = root_poly(end:-1:1);
alp = roots(root_poly_reverse);
F=[];
F_1 = reshape(F1,3,3);
F_2 = reshape(F2,3,3);
F_1 = F_1';
F_2 = F_2';
for i=1:3
    if(isreal(alp(i)))
        tmp = alp(i) * F_1 + F_2;
        tmp = T2' * tmp * T1;
        F = cat(3,F,tmp);
    end;
end;


end

