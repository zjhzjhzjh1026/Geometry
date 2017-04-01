function [ de_r_w ] = jacobian_angle_axis(w)  %%3*1

theta = norm(w);
de_skew_w=[0,0,0;0,0,-1;0,1,0;0,0,1;0,0,0;-1,0,0;0,-1,0;1,0,0;0,0,0];
if(theta>1e-5)
    de_theta_w = w'/theta;
    s=(1-cos(theta))/(theta^2);
    M=w*w';
    de_m_w=kron(w,eye(3))+kron(eye(3),w);
    de_s_w=((theta*sin(theta)-2*(1-cos(theta)))/theta^3)*de_theta_w;
    vec_I=eye(3);
    vec_I=vec_I(:);
    vec_skew=skew(w)';
    vec_skew=vec_skew(:);
    vec_M=M';
    vec_M=vec_M(:);
    de_r_w=-vec_I*sin(theta)*de_theta_w+sinc_override(theta)*de_skew_w+...
        vec_skew*sinc_deri(theta)*de_theta_w+s*de_m_w+vec_M*de_s_w;
else
    de_r_w=de_skew_w;
end;    

end

