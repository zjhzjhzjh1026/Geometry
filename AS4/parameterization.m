function [ v ] = parameterization( v_overline )
% Parameterize a vector
v_overline= v_overline/norm(v_overline,'fro');      %%normalize to 1
len=size(v_overline,2);

a=v_overline(1);
b=v_overline(2:len);

x_sinc=acos(a);
tmp=sinc_override( x_sinc );

v=2/tmp*b;
end

