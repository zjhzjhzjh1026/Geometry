close all;clc;clear;
%% import data
x_in=importdata('hw3_points2D.txt');
X_in=importdata('hw3_points3D.txt');
[row col]=size(x_in);
x_all=[x_in ones(row,1)];  %%inhomogenous -> homogeneous
X_all=[X_in ones(row,1)];
%% Prob a. Outlier Rejection
%% calibrated normalization
K=[1545.0966799187809,0,639.5;0,1545.0966799187809,359.5;0,0,1];    
x_normalized=inv(K)*x_all';
x_normalized=x_normalized';
rng('default');
rng(30);
%% MSAC
    consensus_min_cost = inf ;
    max_trail = inf ;
    trails=0;
    thre=chi2inv(0.95,2)*1;
while(trails < max_trail)
%     rand_num=[1,2,3];
    rand_num=randperm(60);
    p_cam=finsterwald(rand_num(1:3),X_in , x_normalized);
    A=[X_in(rand_num(1),:);X_in(rand_num(2),:);X_in(rand_num(3),:)]';
    if(size(p_cam,3)<2)
        continue;
    end;
    for i=1:size(p_cam,3)
        B = p_cam(:,:,i);
        [Rt,flag] = svd_computeRt( A , B );
        if(flag)
            continue;
        end;
        P = K * Rt;
        error = compute_error( P , X_all , x_in);
        [cost , inlier_num]= compute_cost( error , thre );
        if(cost < consensus_min_cost)
            consensus_min_cost = cost ;
            consensus_min_cost_model = P;
            inlier = inlier_num;
            w = inlier_num / row;
            max_trail = log(1 - 0.99) / log(1 - w^3);
        end;
    end;
    trails = trails + 1;
end;
%% Best Model
error = compute_error (consensus_min_cost_model , X_all , x_in);
inlier_point = get_inlier( error , thre);
x_normalized_inlier= x_normalized(inlier_point,:);
X_inlier = X_in(inlier_point,:)';
X_all_inlier=X_all(inlier_point,:);
count=size(X_inlier,2);
%% Prob b. EPnP Estimation
miu_x = mean(X_inlier,2);
Sigma_x = compute_Sigma( X_inlier , miu_x);
[U , D , V] = svd(Sigma_x);
sigma_x=trace(Sigma_x);
%% control points
s=sqrt(sigma_x/3);
c1=miu_x;
c2=s*V(:,1)+miu_x;
c3=s*V(:,2)+miu_x;
c4=s*V(:,3)+miu_x;
%% parameterization
A_inv=V'/s;
b=X_inlier-c1;
alpha=A_inv*b;
alpha1=ones(1,count)-sum(alpha,1);
alpha=[alpha1;alpha]';
X_para=[];
for i=1:count
    X_temp=alpha(i,1)*c1+alpha(i,2)*c2+alpha(i,3)*c3+alpha(i,4)*c4;
    X_para=[X_para;X_temp'];
end;
%% control point in camera frame
M=[];
for i=1:count
    M_temp=[alpha(i,1),0,-alpha(i,1)*x_normalized_inlier(i,1),...
        alpha(i,2),0,-alpha(i,2)*x_normalized_inlier(i,1),...
        alpha(i,3),0,-alpha(i,3)*x_normalized_inlier(i,1),...
        alpha(i,4),0,-alpha(i,4)*x_normalized_inlier(i,1);
        0,alpha(i,1),-alpha(i,1)*x_normalized_inlier(i,2),...
        0,alpha(i,2),-alpha(i,2)*x_normalized_inlier(i,2),...
        0,alpha(i,3),-alpha(i,3)*x_normalized_inlier(i,2),...
        0,alpha(i,4),-alpha(i,4)*x_normalized_inlier(i,2)];
    M=[M;M_temp];
end;
[U,D,V]=svd(M);
C_cam=reshape(V(:,12),3,4);
%% Depara
X_cam=[];
for i=1:count
    X_temp=alpha(i,1)*C_cam(:,1)+alpha(i,2)*C_cam(:,2)...
        +alpha(i,3)*C_cam(:,3)+alpha(i,4)*C_cam(:,4);
    X_cam=[X_cam;X_temp'];
end;
X_cam=X_cam';
%% Scale
miu_x_cam = mean(X_cam,2);
Sigma_x_cam = compute_Sigma( X_cam , miu_x_cam);
[U , D , V] = svd(Sigma_x_cam);
sigma_x_cam=trace(Sigma_x_cam);
if(miu_x_cam(3,1)<0)
    beta=-sqrt(sigma_x/sigma_x_cam);
else
    beta=sqrt(sigma_x/sigma_x_cam);
end;
X_cam=beta*X_cam;
%% Calculate Rt
Rt = EPnP_Rt( X_inlier , X_cam );
P_DLT = K * Rt;
format longg;
R_DLT=Rt(:,1:3)         %%display
t_DLT=Rt(:,4)
% mse = sum(compute_error( P_DLT ,X_all , x_in))/60
%% Prob c. Nonlinear Estimation
w_hat = angle_axis_para( R_DLT );
w_hat = angle_axis_norm(w_hat);
t=t_DLT; 
%% LM Iteration
lambda=0.001;
k=1;   %%loop control
cost_diff=1000;
count=1;
abc=[];
while(cost_diff>0.0000001)  %%¡¡Loop Control   
    count =count+1;
    if(k==1)
        x_hat=(Rt*X_all_inlier');    %%step 1
        x_hat=x_hat./x_hat(3,:);
        x_hat=x_hat(1:2,:)';
        eps=x_normalized_inlier(:,1:2)-x_hat;
        eps=reshape(eps',2*size(eps,1),1);
        k=2;
        sigma_x_DN = zeros(2*inlier);
        te=inv(K);
        tmp=te(1:2,1:2)*eye(2)*te(1:2,1:2)';
        for i=1:inlier
            sigma_x_DN((2*i-1):2*i,(2*i-1):2*i)=tmp;
        end;
        eps'*inv(sigma_x_DN)*eps
        abc=[abc eps'*inv(sigma_x_DN)*eps];
    end;
    
    if(k==2)
        J = [];                           %%step2
        J  = compute_Jacobian( (Rt*X_all_inlier')' , X_para , w_hat );
        k = 4;
    end;
    
    if(k==4)
        sigma_x_DN = zeros(2*inlier);
        te=inv(K);
        tmp=te(1:2,1:2)*eye(2)*te(1:2,1:2)';
        for i=1:inlier
            sigma_x_DN((2*i-1):2*i,(2*i-1):2*i)=tmp;
        end;
        left= J'*inv(sigma_x_DN)*J+lambda*eye(6);  %%step4
        right=J'*inv(sigma_x_DN)*eps;
        delta=inv(left)*right;    %%??
    end;
    
    w_hat_tmp=w_hat+delta(1:3);  %%step 5
    t_tmp=t+delta(4:6);
    
    R_test=angle_axis_depara( w_hat_tmp );
    P_test=[R_test,t_tmp];
    x0=(P_test*X_all_inlier');    
    x0=x0./x0(3,:);
    x0=x0(1:2,:)';
    eps0=x_normalized_inlier(:,1:2)-x0;
    eps0=reshape(eps0',2*inlier,1);
    eps0'*inv(sigma_x_DN)*eps0
    abc= [abc eps0'*inv(sigma_x_DN)*eps0];
    %%step 7
    cost_diff=eps'*inv(sigma_x_DN)*eps-eps0'*inv(sigma_x_DN)*eps0;
    if(cost_diff>0)
        w_hat=w_hat_tmp;
        t=t_tmp;
        Rt=[R_test,t];
        eps=eps0;
        lambda=0.1*lambda;
        k=2;
    else
        lambda=10*lambda;
        k=4;
    end;   
end
% plot(abc);
w_hat
t
R=angle_axis_depara(w_hat)