function func_B1 = cal_B1
syms x y
x_hat=[x y];
x_bar=deparameterization(x_hat).';
x1_hat_inhomo=x_bar(1:2,:)./x_bar(3,:);
func_B1([x y])=jacobian(x1_hat_inhomo.',x_hat);
end