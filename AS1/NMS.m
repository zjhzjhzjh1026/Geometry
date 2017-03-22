function [ J ] = NMS( im ,win_size )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[row col]=size(im);
win_half=fix(win_size/2);
J=zeros(row,col);

for x=win_half+1:row-win_half
    for y=win_half+1:col-win_half
        window=im(x-win_half:x+win_half,y-win_half:y+win_half);
        maxV=max(max(window));
        J(x,y)=im(x,y)*(im(x,y)==maxV);
    end
end

end

