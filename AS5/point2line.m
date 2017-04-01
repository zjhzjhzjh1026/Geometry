function [ lines ] = point2line( rand_num , num , x1_outlier , F)

x1 = zeros(num,3);

for i=1:num
    x1(i,:) = [x1_outlier(rand_num(i),:) 1];
end;

lines = F * x1';

end

