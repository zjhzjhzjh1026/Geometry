function [ X_tri ] = triangularization_v2(x1 , x2 , wu, wv , s)

F  = fundamental_deparameterization( wu , wv , s );
X_tri=[];
poly_func = matlabFunction(cal_t);
for i = 1 : size(x1,1)
    T_1 = [x1(i,3),0,-x1(i,1);0,x1(i,3),-x1(i,2);0,0,x1(i,3)];
    T_2 = [x2(i,3),0,-x2(i,1);0,x2(i,3),-x2(i,2);0,0,x2(i,3)]; 
    Fs = inv(T_2)' * F * inv(T_1);
    [~ , ~ ,V] = svd(Fs);
    e1 = V(:,end);
    [~ , ~ ,V] = svd(Fs');
    e2 = V(:,end);
    e1_nor = sqrt(1/((e1(1))^2+(e1(2))^2))*e1;
    e2_nor = sqrt(1/((e2(1))^2+(e2(2))^2))*e2;
    R_1 = [e1_nor(1),e1_nor(2),0;-e1_nor(2),e1_nor(1),0;0,0,1];
    R_2 = [e2_nor(1),e2_nor(2),0;-e2_nor(2),e2_nor(1),0;0,0,1];
    Fs = R_2 * Fs * R_1';
    root_poly = poly_func(e1_nor(3),e2_nor(3),Fs(2,2),Fs(2,3),...
        Fs(3,2),Fs(3,3));
    root_poly_reverse = root_poly(end:-1:1);
    t = roots(root_poly_reverse);
    min_cost = inf;
    for j = 1:size(t,1)
        t_r = real(t(j));
        s_t = cal_cost (e1_nor(3),e2_nor(3),Fs(2,2),Fs(2,3),...
        Fs(3,2),Fs(3,3),t_r);
        if(min_cost > s_t)
            min_cost = s_t;
            num = t_r;
        end;
    end;
    a=Fs(2,2);
    b=Fs(2,3);
    c=Fs(3,2);
    d=Fs(3,3);
    t_final = num;
    l_1=[t_final*e1_nor(3),1,-t_final]';
    l_2=[-e2_nor(3)*(c*t_final+d),a*t_final+b,c*t_final+d];
    x1_hat = [-l_1(1)*l_1(3);-l_1(2)*l_1(3);l_1(1)^2+l_1(2)^2];
    x2_hat = [-l_2(1)*l_2(3);-l_2(2)*l_2(3);l_2(1)^2+l_2(2)^2];
    x1_hat = inv(T_1) * R_1' *x1_hat;
    x2_hat = inv(T_2) * R_2' *x2_hat;
    l_epi = F * x1(i,:)';
    l_epr_ver = [-l_epi(2)*x2_hat(3),l_epi(1)*x2_hat(3),l_epi(2)*x2_hat(1)-l_epi(1)*x2_hat(2)];
    
    [P1 , P2 ] = F2P_new( wu , wv , s );
    plane = P2' * l_epr_ver';
    a=plane(1);
    b=plane(2);
    c=plane(3);
    d=plane(4);    
    X_temp = [d*x1_hat(1),d*x1_hat(2),d*x1_hat(3),-(a*x1_hat(1)+b*x1_hat(2)+c*x1_hat(3))];
    X_tri = [X_tri;X_temp];
end;

end

