function r=sG_h_sq(p,u) % compute pde-part of residual
ov=ones(1,p.np);
dw=u(p.nu+15);dh=u(p.nu+16);
h=u(2*p.np+1:3*p.np);
sK = p.mat.sK;
%dh=sparse(2*dh*h)
K=[[sK     N      N];
   [N    dw*sK    N];
   [N      N      N]];
Dx=p.mat.Dx;
part2 = [ov;ov;2*dh*(Dx*u(1+2*p.np:p.nu))];
part1 = [ov;ov;2*dh*h.*(sK*h)];
f=bwh_rhs(p,u); sc=u(p.nu+21); 
Mf=p.mat.M*f;
Ku=sc*(K*u(1:p.nu)+part1+part2);
r=Ku-Mf; 

% (h^2)_xx = 2*h*(h_xx) + 2(h_x)^2

