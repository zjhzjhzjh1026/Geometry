function [ P1 , P2 ] = F2P_new( w_u , w_v , s )

U = angle_axis_depara( w_u );
V = angle_axis_depara( w_v );
sigma = deparameterization (s);

P1 = [eye(3) zeros(3,1)];

sig = [sigma(1),sigma(2),(sigma(1)+sigma(2))/2];
D_new = diag(sig);
Z = [0,-1,0;1,0,0;0,0,1];
M = U * Z * D_new * V'; 
P2 = [M -U(:,3)];

end

