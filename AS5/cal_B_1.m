function [ xd ] = cal_B_1( x_hat , P , X_para )

X_depara =  deparameterization( X_para );
w=P(3,:) * X_depara';
xd=[P(1,:)-x_hat(1)*P(3,:);P(2,:)-x_hat(2)*P(3,:)]/w;

end

