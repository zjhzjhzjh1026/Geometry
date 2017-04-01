function [ cost ] = cal_cost (f_1,f_2,a,b,c,d,t)

cost_1 = t^2/(1+t^2*f_1^2);
cost_2 = (c*t+d)^2/((a*t+b)^2+f_2^2*(c*t+d)^2);
cost = cost_1 + cost_2;

end

