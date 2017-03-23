close all;clc;clear;
%% import data
x_in=importdata('hw2_points2D.txt');
X_in=importdata('hw2_points3D.txt');
[row col]=size(x_in);
x_all=[x_in ones(row,1)];
X_all=[X_in ones(row,1)];

%% normlize data
miux=mean(x_all,1);
varx=var(x_all,1);
sx=sqrt(2/sum(varx));

miuX=mean(X_all,1);
varX=var(X_all,1)*row/(row-1);
sX=sqrt(3/sum(varX));

T=[sx,0,-miux(1)*sx;0,sx,-miux(2)*sx;0,0,1];
U=[sX,0,0,-miuX(1)*sX;0,sX,0,-miuX(2)*sX;0,0,sX,-miuX(3)*sX;0,0,0,1];

x_nor=x_all*T';      
X_nor=X_all*U';

%% Compute Hv
e=[1,0,0];
v=x_nor+sign(x_nor(:,1)).*sqrt((sum(x_nor'.^2)))'*e;
I=eye(3);
A=[];
for i=1:row
    Hv=I-2/(v(i,:)*v(i,:)')*(v(i,:)'*v(i,:));
    m=Hv(2:3,:);
    a=kron(m,X_nor(i,:));
    A=[A ; a];
end

%% SVD  
[~,D,V]=svd(A);
c=size(V,2);
p_nor=reshape(V(:,c),4,3)';
p=inv(T)*p_nor*U;
p_mat= p/norm(p,'fro');      %%normalize
format shortg; 

%%
pT=p_mat';
p=pT(:)';

%% data normalize, cov propagation
sigma_x_DN=T(1,1)^2*eye(2);

pT=p_nor';
p_hat=pT(:)';
p_hat=parameterization(p_hat);  %%1*11
p_hat_norm=normalize(p_hat);

%% LM Iteration
lambda=0.001;
k=1;   %%loop control
cost_diff=1000;
count=1;
while(cost_diff>0.00001)  %%¡¡Terminate    
    count =count+1;
    if(k==1)
        x_hat=(p_nor*X_nor');    %%step 1
        x_hat=x_hat./x_hat(3,:);
        x_hat=x_hat(1:2,:)';
        eps=x_nor(:,1:2)-x_hat;
        eps=reshape(eps',100,1);
        k=2;

    end;
    
    if(k==2)
        J=[];                           %%step2
        for i = 1:row
            J=[J; compute_matrix_deriv(x_hat(i,:) ,...
                p_hat_norm , X_nor(i,:))];
        end;
        k=4;
    end;
    
    if(k==4)
        sigma_x_DN=sx^2*eye(100);
        left= J'*inv(sigma_x_DN)*J+lambda*eye(11);  %%step4
        right=J'*inv(sigma_x_DN)*eps;
        delta=inv(left)*right;   
    end;
    
    p0=p_hat_norm+delta';  %%step 5
    
    p0_de=deparameterization( p0 );%%step 6
    p0_de_nor=reshape(p0_de,4,3)';
    x0=(p0_de_nor*X_nor');    
    x0=x0./x0(3,:);
    x0=x0(1:2,:)';
    eps0=x_nor(:,1:2)-x0;
    eps0=reshape(eps0',100,1);
    %%eps'*inv(sigma_x_DN)*eps
    eps0'*inv(sigma_x_DN)*eps0                            %%step 7
    cost_diff=eps'*inv(sigma_x_DN)*eps-eps0'*inv(sigma_x_DN)*eps0;
    if(cost_diff>0)
        p_hat_norm=p0;
        eps=eps0;
        lambda=0.1*lambda;
        k=2;
    else
        lambda=10*lambda;
        k=4;
    end;
        
end

ptemp=inv(T)*p0_de_nor*U;      
ptemp2= ptemp/norm(ptemp,'fro');
ptemp2