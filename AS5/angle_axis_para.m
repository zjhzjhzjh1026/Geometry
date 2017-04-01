function [ w ] = angle_axis_para( R )

[~,~,V] = svd(R-eye(3));
v=V(:,end);
v_hat = [R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
sin_theta = v'*v_hat/2;
cos_theta = (trace(R)-1)/2;
theta=atan2(sin_theta,cos_theta);
w=theta*v/norm(v);

theta = norm(w);  %%normalization
if(theta>pi)
    w = w * (1 - 2*pi/theta*ceil((theta-pi)/(2*pi)));
end;

end

