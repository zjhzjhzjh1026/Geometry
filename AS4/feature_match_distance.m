function [ matches ] = feature_match_distance( im1 ,im2 ,x1 , y1 , x2 ,...
    y2 , win_size , cor_thre , dis_thre)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[r c]=size(im1);
row=size(x1,1);
col=size(x2,1);
matches=zeros(1,4);  %% store output
cor=zeros(row,col);
win_half=fix(win_size/2);
for i = 1:row
    if(x1(i)<win_half+2||x1(i)>r-win_half-2||y1(i)<win_half+2||y1(i)>c-win_half-2)
        continue;
    end;
    x_ceil=ceil(x1(i));
    y_ceil=ceil(y1(i));
    im1_compute_corr=double(im1(x_ceil-win_half-1:x_ceil+win_half+1,...
        y_ceil-win_half-1:y_ceil+win_half+1));
    X1=meshgrid((x1(i)-x_ceil+2):((x1(i)-x_ceil+win_size+1)));
    Y1=meshgrid((y1(i)-y_ceil+2):((y1(i)-y_ceil+win_size+1)))';
    window_x = interp2(im1_compute_corr,X1,Y1,'linear');
    for j = 1:col
        if(x2(j)<win_half+2||x2(j)>r-win_half-2||y2(j)<win_half+2||...
                y2(j)>c-win_half-2)
            continue;
        end;  
        x2_ceil=ceil(x2(j));
        y2_ceil=ceil(y2(j));
        im2_compute_corr=double(im2(x2_ceil-win_half-1:x2_ceil+...
            win_half+1,y2_ceil-win_half-1:y2_ceil+win_half+1));
        X2=meshgrid((x2(j)-x2_ceil+2):((x2(j)-x2_ceil+win_size+1)));
        Y2=meshgrid((y2(j)-y2_ceil+2):((y2(j)-y2_ceil+win_size+1)))';
        window_y = interp2(im2_compute_corr,X2,Y2,'linear');
        cor(i,j)=corr2(window_x,window_y);
    end
end

max_cor=max(max(cor));
while(max_cor>cor_thre)
    [x y]=find(cor==max_cor);
    cor(x(1),y(1))=-1;
    dist=sqrt((x1(x(1))-x2(y(1)))^2+(y1(x(1))-y2(y(1)))^2);
    if(dist>90) 
        max_cor=max(max(cor));
        continue;
    end;
    next_row=max(cor(x(1),:));
    next_col=max(cor(:,y(1)));
    next_best=max(next_col,next_row);

    if(1-max_cor<(1-next_best)*dis_thre)
       matches=[matches;x1(x(1)),y1(x(1)),x2(y(1)),y2(y(1))];
    end
    cor(:,y(1))=-1;
    cor(x(1),:)=-1;
    max_cor=max(max(cor));
end
[r c]=size(matches);
matches=matches(2:r,:);
end

