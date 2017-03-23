function [ p ] = compute_angel_axis_derivative(w_hat , X_inhomo)

theta = norm(w_hat);
s = (1 - cos(theta))/(theta^2);
a = X_inhomo;
if(abs(theta)<0.000001)
    p = skew(-a);
else
    p = sinc_override(theta)*skew(-a)+cross(w_hat,a)'*sinc_deri(theta)...
        *(w_hat'/norm(w_hat))+cross(w_hat,cross(w_hat,a))'*...
        (sin(theta)*(theta^2)-2*theta+2*theta*cos(theta))/(theta^4)...
        *(w_hat'/norm(w_hat))+s*eye(3)*(skew(w_hat)*skew(-a)+...
        skew(-cross(w_hat,a)));
end;

end

