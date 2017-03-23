function [ y ] = sinc_override( x )

if (x==0)     %%sinc function
    y=1;
else
    y=sin(x)/x;
end;

end

