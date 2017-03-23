function [ v_nor ] = normalize( v )

nor=norm(v);
if(nor>pi)
    v_nor=(1-2*pi/nor*ceil((nor-pi)/pi/2))*v;
else
    v_nor=v;
end;

end

