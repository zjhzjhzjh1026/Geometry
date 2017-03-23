clc;clear;close all;
x_1=importdata('2Dpoints1.txt');
x_2=importdata('2Dpoints2.txt');
n=size(x_1,1);   %%row
x1_ho=[x_1 ones(n,1)];  %%inhomogenous -> homogeneous
x2_ho=[x_2 ones(n,1)];
% rng('default');
rng(10);
%% MSAC
    consensus_min_cost = inf ;
    max_trail = inf ;
    trails=0;
    thre=chi2inv(0.95,2)*1;
while(trails < max_trail)
%     rand_num=[1,2,3,4];
    rand_num=randperm(n);
    H1_inv=four_points_algo(rand_num(1:4),x1_ho);
    H2_inv=four_points_algo(rand_num(1:4),x2_ho);
    H=H2_inv*inv(H1_inv);
    x2_ho_hat=H*x1_ho';
    %%sampson error
    error = sampson_error( x_1 , x_2 , H );
    [cost , inlier_num]= compute_cost( error , thre );
    if(cost < consensus_min_cost)
        consensus_min_cost = cost ;
        consensus_min_cost_model = H;
        inlier = inlier_num;
        w = inlier_num / n;
        max_trail = log(1 - 0.99) / log(1 - w^3);
    end;
    trails = trails + 1;
end;
%% Best Model
error = sampson_error( x_1 , x_2 , consensus_min_cost_model );
inlier_point = get_inlier( error , thre);
x1_inlier = x_1(inlier_point,:);
x2_inlier = x_2(inlier_point,:);
% x1_inlier = feat1_matched_inliers;
% x2_inlier = feat2_matched_inliers;
count=size(x2_inlier,1);
inlier=size(x2_inlier,1);
%% plot
im1=imread('price_center20.JPG');     %%read in data
im2=imread('price_center21.JPG');
window_size=15;
rr=ones(size(x1_inlier,1),1)*window_size;
figure(5);
im_square11 = insertShape(im1,'rectangle',[x1_inlier(:,1)-fix(...
    window_size/2) x1_inlier(:,2)-fix(window_size/2) rr rr],'LineWidth',1);
im_square11 = insertShape(im_square11,'line',[x1_inlier(:,1) x1_inlier(:,2)...
    x2_inlier(:,1) x2_inlier(:,2)],'LineWidth',1);
im_square22 = insertShape(im2,'rectangle',[x2_inlier(:,1)-fix(...
    window_size/2) x2_inlier(:,2)-fix(window_size/2) rr rr],'LineWidth',1);
im_square22 = insertShape(im_square22,'line',[x1_inlier(:,1) x1_inlier(:,2) x2_inlier(:,1) x2_inlier(:,2)],'LineWidth',1);
subplot(1,2,1);
imshow(im_square11);
title('Feature matched in price\_center20.JPG with inliers ');
subplot(1,2,2);
imshow(im_square22);
title('Feature matched in price\_center21.JPG with inliers ');

%% DLT
%%normlize data
x1_inlier_ho=[x1_inlier ones(count,1)];  %%inhomogenous -> homogeneous
x2_inlier_ho=[x2_inlier ones(count,1)];

miux1=mean(x1_inlier_ho,1);
varx1=var(x1_inlier_ho,1);
sx1=sqrt(2/sum(varx1));

miux2=mean(x2_inlier_ho,1);
varx2=var(x2_inlier_ho,1);
sx2=sqrt(2/sum(varx2));

T1=[sx1,0,-miux1(1)*sx1;0,sx1,-miux1(2)*sx1;0,0,1];
T2=[sx2,0,-miux2(1)*sx2;0,sx2,-miux2(2)*sx2;0,0,1];

x1_inlier_nor=x1_inlier_ho*T1';      %%left multiply ->right multiply
x2_inlier_nor=x2_inlier_ho*T2';
%% Compute Hv
e=[1,0,0];
v=x2_inlier_nor+sign(x2_inlier_nor(:,1)).*sqrt((sum(x2_inlier_nor'.^2)))'*e;
I=eye(3);
A=[];
for i=1:count
    Hv=I-2/(v(i,:)*v(i,:)')*(v(i,:)'*v(i,:));
    m=Hv(2:3,:);
    a=kron(m,x1_inlier_nor(i,:));
    A=[A ; a];
end
%% SVD  
[~,D,V]=svd(A);
c=size(V,2);
H_nor=reshape(V(:,c),3,3)';
H=inv(T2)*H_nor*T1;
H_mat= H/norm(H,'fro');      %%normalize to 1
format longg;
H_mat
% x2_inlier_ho_hat=H_mat*x1_inlier_ho';
% x2_inlier_hat=x2_inlier_ho_hat(1:2,:)./x2_inlier_ho_hat(3,:);
% error = sum(sampson_error( x1_inlier , x2_inlier , H_mat ))
%% Levenberg-Marquard Nonlinear Estimation
%% Sampson Correction
x_scene = sampson_correction( x1_inlier , x2_inlier , H_mat );
x_scene = [x_scene ones(size(x_scene,1),1)];
%%vectorize
H_vec=H_mat';
H_vec=H_vec(:)';
%%parameterization&normalization
H_vec_para=normalize(parameterization(H_vec));
x_scene_para=[];
for i =1:size(x_scene,1)
    tmp=normalize(parameterization(x_scene(i,:)));
    x_scene_para=[x_scene_para;tmp];
end;
%% LM Iteration
lambda=0.001;
k=1;   %%loop control
cost_diff=1000;
count=1;
abc=[];
Ai=matlabFunction(cal_A);
B1i=matlabFunction(cal_B1);
B2i=matlabFunction(cal_B2);
while(1)  %%¡¡Loop Control   
    %%use x_inlier_nor
    count =count+1;
    if(k==1) 
        error1=sum(sum((x_scene(:,1:2)-x1_inlier).^2));
        xe=H_mat*x_scene';
        xe=xe(1:2,:)./xe(3,:);
        error2=sum(sum((xe-x2_inlier').^2));
        initial_cost=error1+error2
        error=error1+error2;
        eps1=x1_inlier-x_scene(:,1:2);
        eps2=x2_inlier-xe';
        k=2;
        abc= [abc initial_cost];
    end;
    
    if(k==2)            %%Compute Jacobian
        A = zeros(2,8,inlier);
        B1 = zeros(2,2,inlier);
        B2 = zeros(2,2,inlier);
        for i = 1:inlier         
            A(:,:,i)=Ai(H_vec_para(1),H_vec_para(2),H_vec_para(3),H_vec_para(4)...
                ,H_vec_para(5),H_vec_para(6),H_vec_para(7),H_vec_para(8),x_scene_para(i,1),x_scene_para(i,2));
            B1(:,:,i)=B1i(x_scene_para(i,1),x_scene_para(i,2));
            B2(:,:,i)=B2i(H_vec_para(1),H_vec_para(2),H_vec_para(3),H_vec_para(4)...
                ,H_vec_para(5),H_vec_para(6),H_vec_para(7),H_vec_para(8),x_scene_para(i,1),x_scene_para(i,2));
        end;
        U=zeros(8,8);
        V=zeros(2,2,inlier);
        W=zeros(8,2,inlier);
        for i = 1 :inlier
            U = U + A(:,:,i)'*eye(2)*A(:,:,i);
            V(:,:,i) = B1(:,:,i)' * eye(2) * B1(:,:,i)+ B2(:,:,i)' * eye(2) * B2(:,:,i);
            W(:,:,i) = A(:,:,i)'*eye(2)*B2(:,:,i);
        end;
        k = 4;
    end;
    
    if(k==4)
        eps_a = zeros(8,1);
        eps_b = zeros(2,1,inlier);
        for i = 1 : inlier
            eps_a = eps_a + A(:,:,i)'*eye(2)*eps2(i,:)';
            eps_b(:,:,i) = B1(:,:,i)'*eye(2)*eps1(i,:)'+ B2(:,:,i)'*eye(2)*eps2(i,:)';
        end;
        S = U + lambda*eye(8);
        e = eps_a;
        for i = 1 : inlier
            S = S - W(:,:,i) * inv(V(:,:,i) +lambda*eye(2)) * W(:,:,i)';
            e = e - W(:,:,i) * inv(V(:,:,i) +lambda*eye(2)) * eps_b(:,:,i);
        end;
        delta_h = inv(S)*e;
        delta_x = [];
        for i=1 : inlier
            tmp = inv(V(:,:,i) +lambda*eye(2))*(eps_b(:,:,i)-W(:,:,i)'*delta_h);
            delta_x = [delta_x tmp];
        end;
    end;
    
    H_vec_para_tmp=H_vec_para+delta_h';  %% step 5
    x_scene_para_tmp = x_scene_para + delta_x';
    
    H_mat_tmp=deparameterization(H_vec_para_tmp);
    H_mat_tmp=reshape(H_mat_tmp,3,3)';
    X_scene_tmp=[];
    for i = 1 : inlier
        tmp=deparameterization(x_scene_para_tmp(i,:));
        X_scene_tmp=[X_scene_tmp;tmp];
    end;
    X_scene_tmp=X_scene_tmp';
    
    X_scene_tmp_imhomo = X_scene_tmp(1:2,:)./X_scene_tmp(3,:);
    x2_hat_homo=H_mat_tmp * X_scene_tmp;
    x2_hat_imhomo=x2_hat_homo(1:2,:)./x2_hat_homo(3,:);
    
    error1_tmp=sum(sum((X_scene_tmp_imhomo'-x1_inlier).^2));
    error2_tmp=sum(sum((x2_hat_imhomo-x2_inlier').^2));
    error_tmp = error1_tmp + error2_tmp
    X_scene_tmp=X_scene_tmp';  %%change back
    
    abc= [abc error_tmp];
    %%step 7
    cost_diff= error - error_tmp; 
    if(cost_diff>0)
        H_vec_para = H_vec_para_tmp ;
        x_scene_para = x_scene_para_tmp;
        eps1=x1_inlier-X_scene_tmp_imhomo';
        eps2=x2_inlier-x2_hat_imhomo';
        error = error_tmp;
        lambda=0.1*lambda;
        k=2;
        if(cost_diff<0.0000001)
            break;
        end;
    else
        lambda=10*lambda;
        k=4;
    end;   
end
format longg;
H_LM_mat=deparameterization(H_vec_para);
H_LM_mat=reshape(H_LM_mat,3,3)';
H_LM_mat= H_LM_mat/norm(H_LM_mat,'fro')