function [ error ] =sampson_error_F( x1_ho , x2_ho , F )
%For fundamental matrix
error = [];
for i = 1 : size(x1_ho , 1)
    numerator = x2_ho(i,:) * F * x1_ho(i,:)';
    denominator=(x2_ho(i,:)*F(:,1))^2 + (x2_ho(i,:)*F(:,2))^2 +...
        (F(1,:)*x1_ho(i,:)')^2 + (F(2,:)*x1_ho(i,:)')^2;
    temp = numerator^2/denominator;
    error = [error;temp];
end;

end

