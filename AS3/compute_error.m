function [ error ] = compute_error( P , X_all , x_in)

x_hat_homo=P*X_all';
x_hat_homo=(x_hat_homo./x_hat_homo(3,:))';
x_hat_inhomo=x_hat_homo(:,1:2);
error=sum((x_hat_inhomo-x_in).^2,2);

end

