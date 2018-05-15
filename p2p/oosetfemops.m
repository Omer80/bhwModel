function p=oosetfemops(p) % in problem-dir, since highly problem dependent
grid=p.pdeo.grid; par=p.u(p.nu+1:end); dw=par(15); dh=par(16);

[K,M,~]=p.pdeo.fem.assema(grid,1,1,1); 
%[Q,G,H,R]=p.pdeo.fem.assemb(grid); 
N=sparse(grid.nPoints, grid.nPoints);
D=diag([1 1 1]); p.mat.M=kron(D,M);  po=getpte(p);
p.mat.K=[[K     N      N]; ...
         [N    dw*K    N]; ...
         [N     N   dh*K]]; 
Dx=makeDx(p); p.mat.Dx=Dx;
p.mat.sK=K;

%if size(po,1)==1
%    Kx=convection(p.pdeo.fem,grid,1); p.mat.Kx=kron(D,Kx); p.mat.Ky=0*p.mat.Kx;
%else
%    Kx=convection(p.pdeo.fem,grid,[1;0]); Ky=convection(p.pdeo.fem,grid,[0;1]);
%    p.mat.Kx=kron(D,Kx); p.mat.Ky=kron(D,Ky);
end
