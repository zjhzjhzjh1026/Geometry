function [ J ] = compute_matrix_deriv( x_hat , p_hat_norm , X_nor)
n=12;
nor=norm(p_hat_norm);  %%compute norm

p_mat=deparameterization(p_hat_norm);  %%compute w
p_mat=reshape(p_mat,4,3)';
w=p_mat(3,:) * X_nor';

xd=[X_nor,zeros(1,4),-x_hat(1)*X_nor;zeros(1,4),X_nor,-x_hat(2)*X_nor]/w;
a=sinc_override(nor/2)/(-4)*p_hat_norm;
if(nor==0)
    b=eye(n-1)/2;
else
    b=sinc_override(nor/2)/2*eye(n-1)+sinc_deri(nor/2)/4/nor*p_hat_norm'*p_hat_norm;
end;
J=xd*[a;b];
end

