function [ pp ] = compute_pi( u,m,n,a2,cos_a,j1,j2,j3 )

v=u*m+n;

s1=sqrt(a2/(u^2+v^2-2*u*v*cos_a)); 
s2=u*s1;
s3=v*s1;

p1=s1*j1;
p2=s2*j2;
p3=s3*j3;
pp=[p1,p2,p3];

end

