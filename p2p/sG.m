function r=sG(p,u) % compute pde-part of residual
f=bwh_rhs(p,u); sc=u(p.nu+21); 
Mf=p.mat.M*f;
Ku=sc*p.mat.K*u(1:p.nu);
r=Ku-Mf; 


