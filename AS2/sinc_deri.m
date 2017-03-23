function [ y ] = sinc_deri( x )

if(x==0)
    y=0;
else
    y=cos(x)/x-sin(x)/x/x;

end

