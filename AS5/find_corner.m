function [x_corner y_corner]=find_corner( im , x , y ,window_size);
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
len=size(x,1);
kernel=[-1,8,0,-8,1]/12;
win_half=fix(window_size/2);
Ix=double(imfilter(im,kernel,'symmetric'));
Iy=double(imfilter(im,kernel','symmetric'));
x_corner=zeros(len,1);
y_corner=zeros(len,1);
for i = 1:len  
    A=zeros(2,2);
    B=zeros(2,1);
    window_x=double(Ix(x(i)-win_half:x(i)+win_half,y(i)-win_half:y(i)+win_half));
    window_y=double(Iy(x(i)-win_half:x(i)+win_half,y(i)-win_half:y(i)+win_half));
    A(1,1)=sum(sum(window_x.^2));
    A(1,2)=sum(sum(window_x.*window_y));
    A(2,1)=A(1,2);
    A(2,2)=sum(sum(window_y.^2));
    tmp1=zeros(window_size,window_size);
    tmp2=zeros(window_size,window_size);
    for m=x(i)-win_half:x(i)+win_half
        for n=y(i)-win_half:y(i)+win_half
            tmp1(m-x(i)+win_half+1,n-y(i)+win_half+1)=m*(Ix(m,n)^2)+n*Ix(m,n)*Iy(m,n);
            tmp2(m-x(i)+win_half+1,n-y(i)+win_half+1)=n*(Iy(m,n)^2)+m*Ix(m,n)*Iy(m,n);
        end
    end
    B(1,1)=sum(sum(tmp1));
    B(2,1)=sum(sum(tmp2));
    t=A\B;
    x_corner(i)=t(1);
    y_corner(i)=t(2);
end

end

