function [ v_overline ] = deparameterization( v )
a=cos(norm(v)/2);

x_sinc=norm(v)/2;
tmp=sinc_override( x_sinc );
b=tmp/2*v;

v_overline=[a,b];
end

