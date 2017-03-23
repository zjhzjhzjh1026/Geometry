function [ J ] = compute_Jacobian( x_homo_hat , X_inhomo , w_hat )

J=[];
for i = 1 : size(x_homo_hat,1)
    p1=[1/x_homo_hat(i,3),0,-x_homo_hat(i,1)/(x_homo_hat(i,3)^2);...
        0,1/x_homo_hat(i,3),-x_homo_hat(i,2)/(x_homo_hat(i,3)^2)];
    p2=eye(3);
    p3=compute_angel_axis_derivative(w_hat , X_inhomo(i,:));
    J1=p1*p2*p3;
    J2=p1*eye(3);
    J=[J;J1,J2];
end

end

