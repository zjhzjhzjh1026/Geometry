function [ w ] = angle_axis_norm( w )

theta = norm(w);
if(theta>pi)
    w = w * (1 - 2*pi/theta*ceil((theta-pi)/(2*pi)));
end;

end

