function [ w_u , w_v , s ] = fundamental_parameterization (F_DLT)

[U , D , V] = svd(F_DLT);
sigma = [D(1,1),D(2,2)];
if (det(U)<0)
    U = -U;
end;
if(det(V)<0)
    V = -V;
end;
w_u = angle_axis_para(U);
w_v = angle_axis_para(V);
s = parameterization(sigma);

end

