clear;clc;
im1=imread('price_center20.JPG');     %%read in data
im2=imread('price_center21.JPG');
figure(1);
subplot(1,2,1);imshow(im1);title('price\_center20.JPG');
subplot(1,2,2);imshow(im2);title('price\_center21.JPG');

im1_grey=rgb2gray(im1);
im2_grey=rgb2gray(im2);

window_size=9;
I1=gradient_filter(im1_grey,window_size);
I2=gradient_filter(im2_grey,window_size);
thre=6.5;
I1=I1.*(I1>thre);
I2=I2.*(I2>thre);
J1=NMS(I1,window_size);
J2=NMS(I2,window_size);
%%sum(sum(J1>0))
[x1,y1]=find(J1>0);
[x2,y2]=find(J2>0);

[x1_corner y1_corner]=find_corner( im1_grey , x1 , y1 , window_size);
[x2_corner y2_corner]=find_corner( im2_grey , x2 , y2 , window_size);
len1=size(x1_corner,1);
len2=size(x2_corner,1);
r1=ones(len1,1)*window_size;
r2=ones(len2,1)*window_size;
im_square1 = insertShape(im1,'rectangle',[y1-fix(window_size/2) x1-...
    fix(window_size/2) r1 r1],'LineWidth',1);
im_square2 = insertShape(im2,'rectangle',[y2-fix(window_size/2) x2-...
    fix(window_size/2) r2 r2],'LineWidth',1);
figure(2);
imshow(im_square1);
title('price\_center20.JPG Feature Detection Window Size=9');
figure(3);
imshow(im_square2);
title('price\_center21.JPG Feature Detection Window Size=9');
%%%%%%%%%%%%%%%%%%%%%%%%%%Problem 1
cor_thre=0.75;
dis_thre=0.92;
window_size2=15;
matches = feature_match_distance( im1_grey ,im2_grey ,x1_corner ,...
    y1_corner , x2_corner , y2_corner , window_size2, cor_thre ,dis_thre);

len3=size(matches,1);
rr=ones(len3,1)*window_size2;

figure(4);
im_square11 = insertShape(im1,'rectangle',[matches(:,2)-fix(...
    window_size2/2) matches(:,1)-fix(window_size2/2) rr rr],'LineWidth',1);
im_square11 = insertShape(im_square11,'line',[matches(:,2) matches(:,1)...
    matches(:,4) matches(:,3)],'LineWidth',1);
im_square22 = insertShape(im2,'rectangle',[matches(:,4)-fix(...
    window_size2/2) matches(:,3)-fix(window_size2/2) rr rr],'LineWidth',1);
im_square22 = insertShape(im_square22,'line',[matches(:,4) matches(:,3) matches(:,2) matches(:,1)],'LineWidth',1);
subplot(1,2,1);
imshow(im_square11);
title('Feature matched in price\_center20.JPG with Proximity Window ');
subplot(1,2,2);
imshow(im_square22);
title('Feature matched in price\_center21.JPG with Proximity Window ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%Problem 2
p1=[matches(:,2) matches(:,1)];
save('2Dpoints1.txt','p1','-ascii');
p2=[matches(:,4) matches(:,3)];
save('2Dpoints2.txt','p2','-ascii');