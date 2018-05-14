function p=init(p,lx,nx,par,ndim,ic) 
p=bwh_stanparam(p); screenlayout(p); % set standard parameters and screenlayout 
p.nc.neq=3; p.nc.ilam=1; % number eq, cont-param
p.fuha.sG=@sG;p.fuha.sGjac=@sGjac; 
p.sw.jac=0; % Numerical Jacobian 0, Analytical Jacobian 1
switch ndim % domain and BC, depending on spatial dim 
    case 1;pde=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; p.DimVol = 2*lx*par(20); % Calculating the dimensional domain 1D
        bc=pde.grid.neumannBC('0'); p.x0i=10; % index for ts-plot
        p.plot.cl={'black','blue'}; p.plot.pcmp=[1 2];
    case 2; pde=stanpdeo2D(lx,lx,2*lx/nx); p.vol=2*lx^2; p.DimVol = (2*lx*par(20))^2; % Calculating the dimensional domain 2D
        bc=pde.grid.neumannBC('0'); p.x0i=30; p.plot.pstyle=2; p.plot.cm=hot; 
        p.plot.axis='image'; 
    case 3 
        ny=round(nx/sqrt(3));sw=[];sw.sym=2;
        pde=stanpdeo2D(lx,lx/sqrt(3),nx,ny,sw);
        p.vol=2*lx^2; p.DimVol = (2*lx*par(20)/sqrt(3))^2; % Calculating the dimensional domain 2D
        bc=pde.grid.neumannBC('0'); p.x0i=30; p.plot.pstyle=2; p.plot.cm=hot; 
        p.plot.axis='image'; 
    case 4; pde=stanpdeo3D(lx,lx/2,lx/4,2*lx/nx); p.vol=0.5*lx^3; 
        p.plot.ng=20; p.plot.lev=[-0.1 0.1]; p.x0i=200; 
        p.plot.levc={'blue','red'}; p.plot.alpha=0.5; 
        bc=pde.grid.dirichletBC('1','0'); 
end 
pde.grid.makeBoundaryMatrix(bc); p.nc.sf=1e3; % OOPDE setting of BC 
p.pdeo=pde; p.sw.sfem=-1; p.np=pde.grid.nPoints; p.nu=p.np*p.nc.neq; 
p.sol.xi=1/p.nu;
switch ic
    case 1 % bare soil        
        b1=0*ones(p.np,1); b2=b1; s1=0.124215276835*ones(p.np,1); s2=s1;
    case 2 % random        
        b1=0.1*rand(p.np,1)+0.4; b2=0.1*rand(p.np,1)+0.4; s1=0.1*rand(p.np,1)+0.1; s2=0.1*rand(p.np,1)+0.1; 
    case 3 % random b1, b2 zero        
        b1=rand(p.np,1); b2=0*ones(p.np,1); s1=0.2*rand(p.np,1); s2=s1; 
    case 4 % random b2, b1 zero        
        b1=0*ones(p.np,1); b2=rand(p.np,1); s1=0.2*rand(p.np,1); s2=s1;
    case 5
        b1=0.1*ones(p.np,1); b2=b1; s1=0.124215276835*ones(p.np,1); s2=s1;
    case 6
        b1=0*ones(p.np,1); b2=0.1*ones(p.np,1); s1=0.124215276835*ones(p.np,1); s2=s1;        
end

p.u=[b2;s1;s2; par']; % initial guess (here trivial) and pars
p=setfemops(p);% setfemops calls oosetfemops in problem dir 





