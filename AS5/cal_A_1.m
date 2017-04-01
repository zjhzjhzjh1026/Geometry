function [ xd ] = cal_A_1( x1_hat , P , X_para )

X_depara =  deparameterization( X_para );
w=P(3,:) * X_depara';
xd=[X_depara,zeros(1,4),-x1_hat(1)*X_depara;zeros(1,4),X_depara,-x1_hat(2)*X_depara]/w;

end

