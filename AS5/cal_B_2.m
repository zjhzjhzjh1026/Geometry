function func_B = cal_B_2

syms X Y Z
X_3D=[X Y Z];
X_depa = deparameterization( X_3D ).';
func_B([X Y Z])=jacobian(X_depa,X_3D);

end