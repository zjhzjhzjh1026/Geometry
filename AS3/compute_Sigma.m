function [ Sigma_x ] = compute_Sigma( X_inlier , miu_x)

Sigma_x=zeros(3,3);
for m=1:size(X_inlier,2)
    Sigma_x =Sigma_x + (X_inlier(:,m)-miu_x)*(X_inlier(:,m)-miu_x)';
end;
Sigma_x=Sigma_x/size(X_inlier,2);

end

