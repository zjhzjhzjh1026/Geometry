function func_B2 = cal_B2
syms h1 h2 h3 h4 h5 h6 h7 h8
h_para=[h1 h2 h3 h4 h5 h6 h7 h8];
H=reshape(deparameterization(h_para).',3,3).';
syms x y
x_hat=[x y];
x_bar=deparameterization(x_hat).';
x1_hat_inhomo=x_bar(1:2,:)./x_bar(3,:);
x2_hat=H*x_bar;
x2_hat_inhomo=x2_hat(1:2,:)./x2_hat(3,:);
func_B2([h1 h2 h3 h4 h5 h6 h7 h8 x y])=jacobian(x2_hat_inhomo.',x_hat);
end