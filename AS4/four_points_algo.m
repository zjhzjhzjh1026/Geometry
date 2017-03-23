function [ H_inv ] = four_points_algo( label , x )

x1=x(label(1),:);
x2=x(label(2),:);
x3=x(label(3),:);
x4=x(label(4),:);

x=[x1' x2' x3'];
lambda=x\x4';
H_inv=[lambda(1)*x1' lambda(2)*x2' lambda(3)*x3'];

end

