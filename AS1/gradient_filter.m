function [ I ] = gradient_filter( im , win_size )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[row col]=size(im);
win_half=fix(win_size/2);
kernel=[-1,8,0,-8,1]/12;
Ix=imfilter(im,kernel,'symmetric');
Iy=imfilter(im,kernel','symmetric');
I=zeros(row,col);
sqnum=win_size^2;

for x = win_half + 1 : row - win_half 
    for y = win_half + 1 : col - win_half
        N=zeros(2,2);
        window_x=double(Ix(x-win_half:x+win_half,y-win_half:y+win_half));
        window_y=double(Iy(x-win_half:x+win_half,y-win_half:y+win_half));
        N(1,1)=sum(sum(window_x.^2))/sqnum;
        N(1,2)=sum(sum(window_x.*window_y))/sqnum;
        N(2,1)=N(1,2);
        N(2,2)=sum(sum(window_y.^2))/sqnum;
        I(x,y)=(trace(N)-sqrt(trace(N)^2-4*det(N)))/2;
    end
end

end

