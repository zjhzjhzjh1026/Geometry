function func_root = cal_t
syms f_1 f_2 a b c d t
g=t*((a*t+b)^2+(f_2^2)*(c*t+d)^2)^2-(a*d-b*c)*((1+(f_1^2)*(t^2))^2)*(a*t+b)*(c*t+d);
equ = collect(g,t);
func_root([f_1 f_2 a b c d])=coeffs(equ,t);
end