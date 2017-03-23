clc;clear;close all;
x_1=importdata('2Dpoints1.txt');
x_2=importdata('2Dpoints2.txt');
n=size(x_1,1);   %%row
x1_ho=[x_1 ones(n,1)];  %%inhomogenous -> homogeneous
x2_ho=[x_2 ones(n,1)];
%% MSAC
consensus_min_cost = inf ;
max_trail = inf ;
rng(8);
trails=0;
thre=chi2inv(0.95,1)*1;
poly_func = matlabFunction(cal_root);
while(trails < max_trail)
%     rand_num=[1,2,3,4,5,6,7];
    rand_num=randperm(n);
    F_all=seven_points_algo(rand_num(1:7),x1_ho,x2_ho,poly_func);
    %%sampson error
    for i = 1 : size(F_all,3)
        F = F_all(:,:,i);
        error = sampson_error_F( x1_ho , x2_ho , F );
        [cost , inlier_num]= compute_cost( error , thre );      
        if(cost < consensus_min_cost)
            consensus_min_cost = cost ;
            consensus_min_cost_model = F;
            inlier = inlier_num;
            w = inlier_num / n;
            max_trail = log(1 - 0.99) / log(1 - w^7);
        end;
    end;
    trails = trails + 1;
end;
% Best Model
error = sampson_error_F( x1_ho , x2_ho , consensus_min_cost_model );
[inlier_point,outlier_point] = get_inlier( error , thre);
x1_inlier = x_1(inlier_point,:);
x2_inlier = x_2(inlier_point,:);
x1_outlier = x_1(outlier_point,:);
x2_outlier = x_2(outlier_point,:);

count=size(x2_inlier,1);
inlier=size(x2_inlier,1);
%% plot
im1=imread('IMG_5030.JPG');     %%read in data
im2=imread('IMG_5031.JPG');
window_size=11;
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
title('Feature matched in IMG\_5030.JPG with inliers ');
subplot(1,2,2);
imshow(im_square22);
title('Feature matched in IMG\_5031.JPG with inliers ');

%% DLT
%%normlize data
x1_inlier_ho=[x1_inlier ones(count,1)];  %%inhomogenous -> homogeneous
x2_inlier_ho=[x2_inlier ones(count,1)];

miux1=mean(x1_inlier_ho,1);
varx1=var(x1_inlier_ho);
sx1=sqrt(2/sum(varx1));

miux2=mean(x2_inlier_ho,1);
varx2=var(x2_inlier_ho);
sx2=sqrt(2/sum(varx2));

T1=[sx1,0,-miux1(1)*sx1;0,sx1,-miux1(2)*sx1;0,0,1];
T2=[sx2,0,-miux2(1)*sx2;0,sx2,-miux2(2)*sx2;0,0,1];

x1_inlier_nor=x1_inlier_ho*T1';      %%left multiply ->right multiply
x2_inlier_nor=x2_inlier_ho*T2';

kron_bin = [];
for i=1:inlier
    tmp = kron(x2_inlier_nor(i,:),x1_inlier_nor(i,:));
    kron_bin = [kron_bin ; tmp];
end;
[~,~,V] = svd(kron_bin);
F_DLT_nor = reshape(V(:,end),3,3);
F_DLT_nor = F_DLT_nor';

[U,D,V] = svd(F_DLT_nor);   %%Rank Constraint
D(3,3) = 0;
F_DLT_nor_new = U*D*V';
F_DLT = T2' * F_DLT_nor_new * T1;
F_DLT = F_DLT /norm(F_DLT ,'fro');      %%normalize to 1
format longg;
F_DLT

%% Levenberg-Marquard Nonlinear Estimation
%% Triangularization&parameterization
[ w_u , w_v , s ] = fundamental_parameterization (F_DLT);
% [ F_DLT ] = fundamental_deparameterization( w_u , w_v , s );
X_tri = triangularization_v2(x1_inlier_ho , x2_inlier_ho , w_u , w_v , s );
[ P1 , P2 ] = F2P_new( w_u , w_v , s );
X_para=[];
for i =1:size(X_tri,1)
    tmp=normalize(parameterization(X_tri(i,:)));
    X_para=[X_para;tmp];
end;

%% LM Iteration
lambda=0.001;
k=1;   %%loop control
cost_diff=1000;
count=1;
abc=[];

Bi=matlabFunction(cal_B_2);
while(1)  %%¡¡Loop Control   
    count =count+1;
    if(k==1) 
        X_depara=[];
        for i =1:size(X_tri,1)
            tmp=deparameterization(X_para(i,:));
            X_depara=[X_depara;tmp];
        end;
        x1_hat = P1 * X_depara'; 
        x2_hat = P2 * X_depara';
        x1_hat = x1_hat(1:2,:)./x1_hat(3,:);
        eps1 = x1_inlier - x1_hat';
        x2_hat = x2_hat(1:2,:)./x2_hat(3,:);
        eps2 = x2_inlier - x2_hat';
        initial_cost=sum(sum(eps1.^2))+sum(sum(eps2.^2))
        
        error=initial_cost;
        k=2;
        abc= [abc initial_cost];
    end;
    
    if(k==2)            %%Compute Jacobian
        A = zeros(2,7,inlier);
        B1 = zeros(2,3,inlier);
        B2 = zeros(2,3,inlier);
        for i = 1:inlier          
            A(:,:,i)=cal_A_1( x2_hat(:,i)' , P2 , X_para(i,:))*func_jacob_direct ( w_u , w_v , s );
            B1(:,:,i)=cal_B_1( x1_hat(:,i)' , P1 , X_para(i,:))*Bi(X_para(i,1),X_para(i,2),X_para(i,3));
            B2(:,:,i)=cal_B_1( x2_hat(:,i)' , P2 , X_para(i,:))*Bi(X_para(i,1),X_para(i,2),X_para(i,3));
        end;
        U=zeros(7,7);
        V=zeros(3,3,inlier);
        W=zeros(7,3,inlier);
        for i = 1 :inlier
            U = U + A(:,:,i)'*eye(2)*A(:,:,i);
            V(:,:,i) = B1(:,:,i)' * eye(2) * B1(:,:,i)+ B2(:,:,i)' * eye(2) * B2(:,:,i);
            W(:,:,i) = A(:,:,i)'*eye(2)*B2(:,:,i);
        end;
         k = 4;
    end;
    
    if(k==4)
        eps_a = zeros(7,1);
        eps_b = zeros(3,1,inlier);
        for i = 1 : inlier
            eps_a = eps_a + A(:,:,i)'*eye(2)*eps2(i,:)';
            eps_b(:,:,i) = B1(:,:,i)'*eye(2)*eps1(i,:)'+ B2(:,:,i)'*eye(2)*eps2(i,:)';
        end;
        S = U + lambda*eye(7);
        e = eps_a;
        for i = 1 : inlier
            S = S - W(:,:,i) * inv(V(:,:,i) +lambda*eye(3)) * W(:,:,i)';
            e = e - W(:,:,i) * inv(V(:,:,i) +lambda*eye(3)) * eps_b(:,:,i);

        end;
        delta_a = S\e;
        delta_x = [];
        for i=1 : inlier
            tmp = inv(V(:,:,i) +lambda*eye(3))*(eps_b(:,:,i)-W(:,:,i)'*delta_a);
            delta_x = [delta_x tmp];
        end;
    end;
    
    w_u_tmp = w_u + delta_a(1:3);  %% step 5
    w_v_tmp = w_v + delta_a(4:6);
    s_tmp = s + delta_a(7);
    X_para_tmp = X_para + delta_x';
    
    X_depara_tmp=[];
    for i = 1 : inlier
        tmp=deparameterization(X_para_tmp(i,:));
        X_depara_tmp=[X_depara_tmp;tmp];
    end;
    
    [ P1_tmp , P2_tmp ] = F2P_new( w_u_tmp , w_v_tmp , s_tmp );

    x1_hat_tmp = P1_tmp * X_depara_tmp'; 
    x2_hat_tmp = P2_tmp * X_depara_tmp';
    x1_hat_tmp = x1_hat_tmp(1:2,:)./x1_hat_tmp(3,:);
    eps1_tmp = x1_inlier - x1_hat_tmp';
    x2_hat_tmp = x2_hat_tmp(1:2,:)./x2_hat_tmp(3,:);
    eps2_tmp =  x2_inlier - x2_hat_tmp';
    error_tmp=sum(sum(eps1_tmp.^2))+sum(sum(eps2_tmp.^2));
    %%step 7
    cost_diff= error - error_tmp; 
    if(cost_diff>0)
        error_tmp
        abc= [abc error_tmp];
        w_u = w_u_tmp ;
        w_v = w_v_tmp ;
        s = s_tmp ;
        P2 = P2_tmp;
        X_para = X_para_tmp;
        eps1=eps1_tmp;
        eps2=eps2_tmp;
        error = error_tmp;
        lambda=0.1*lambda;
        k=2;
        if(cost_diff<1e-5)
            break;
        end;
    else
        lambda=10*lambda;
        k=4;
    end;   
end
format longg;
rng(30);
[ F_LM ] = fundamental_deparameterization( w_u , w_v , s );
F_LM= F_LM/norm(F_LM,'fro')

%% Point to line mapping 
cnt_outlier = size(x1_outlier,1);  %% the number of outliers  
num = 3;
rand_num=randperm(cnt_outlier);
lines=point2line( rand_num , num , x1_outlier , F_LM);

%% Plot
figure(6);
rr=ones(1,1)*window_size;
im_square11 = insertShape(rgb2gray(im1),'rectangle',[x1_outlier(rand_num(1),1)-fix(...
    window_size/2) x1_outlier(rand_num(1),2)-fix(window_size/2) rr rr],'LineWidth',1,'Color','red', 'LineWidth', 3);
im_square11 = insertShape(im_square11,'rectangle',[x1_outlier(rand_num(2),1)-fix(...
    window_size/2) x1_outlier(rand_num(2),2)-fix(window_size/2) rr rr],'LineWidth',1,'Color','green', 'LineWidth', 3);
im_square11 = insertShape(im_square11,'rectangle',[x1_outlier(rand_num(3),1)-fix(...
    window_size/2) x1_outlier(rand_num(3),2)-fix(window_size/2) rr rr],'LineWidth',1,'Color','blue', 'LineWidth', 3);

im_square22 = insertShape(rgb2gray(im2),'rectangle',[x2_outlier(rand_num(1),1)-fix(...
    window_size/2) x2_outlier(rand_num(1),2)-fix(window_size/2) rr rr],'LineWidth',1,'Color','red', 'LineWidth', 3);
im_square22 = insertShape(im_square22,'rectangle',[x2_outlier(rand_num(2),1)-fix(...
    window_size/2) x2_outlier(rand_num(2),2)-fix(window_size/2) rr rr],'LineWidth',1,'Color','green', 'LineWidth', 3);
im_square22 = insertShape(im_square22,'rectangle',[x2_outlier(rand_num(3),1)-fix(...
    window_size/2) x2_outlier(rand_num(3),2)-fix(window_size/2) rr rr],'LineWidth',1,'Color','blue', 'LineWidth', 3);

subplot(1,2,1);
imshow(im_square11);
title('3 outliers in IMG\_5030.JPG');
subplot(1,2,2);
imshow(im_square22);
hold on;
width=1:size(im2,2);
i=1;
y = (-lines(3,i)-lines(1,i)*width)/lines(2,i);
line(width, y, 'Color','red', 'LineWidth', 1);
i=2;
y = (-lines(3,i)-lines(1,i)*width)/lines(2,i);
line(width, y, 'Color','green', 'LineWidth', 1);
i=3;
y = (-lines(3,i)-lines(1,i)*width)/lines(2,i);
line(width, y, 'Color','blue', 'LineWidth', 1);

title('Point to line mapping in IMG\_5031.JPG');