function [ F ] = fundamental_deparameterization( w_u , w_v , s )

U = angle_axis_depara( w_u );
V = angle_axis_depara( w_v );
sigma = deparameterization (s);
F = sigma(1)*U(:,1)*V(:,1)'+sigma(2)*U(:,2)*V(:,2)';

end

