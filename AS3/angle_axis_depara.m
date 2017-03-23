function [ R ] = angle_axis_depara( w )

theta=norm(w);
R=cos(theta)*eye(3)+sinc_override(theta)*skew(w)...
    +(1-cos(theta))/(theta^2)*w*w';

end

