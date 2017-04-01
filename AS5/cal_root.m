function func_root = cal_root 
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 b1 b2 b3 b4 b5 b6 b7 b8 b9 alp 
F1_sym=[a1,a2,a3;a4,a5,a6;a7,a8,a9];
F2_sym=[b1,b2,b3;b4,b5,b6;b7,b8,b9];
F = alp * F1_sym + F2_sym ;
detF = det(F);
equ = collect(detF,alp);
func_root([a1 a2 a3 a4 a5 a6 a7 a8 a9 b1 b2 b3 b4 b5 b6 b7 b8 b9])=coeffs(equ,alp);
end