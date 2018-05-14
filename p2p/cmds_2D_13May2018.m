close all; format compact; keep pphome;clc;
%% init - setting of 14.86x(14.86/sqrt(3)) domain size from AUTO prediction
ic=1 ; % choose initial condition - 1 bare soil, 2 random, 3 b1 random, 4 b2 random
ndim=1; % For 2D domain size of lx * lx/sqrt(3)
dir='1Dbs'; p=[]; lx=14.862921667; % From auto continuation over af_set12
nx=50; % domain size and spat.resolution 
% loading parameters from file
par=loadparms('bwh_set2.mat'); par(1)=0.1; 
p=init(p,lx,nx,par,ndim,ic); p=setfn(p,dir); 
p.nc.dsmax=0.05; r=resi(p,p.u); norm(r,inf); p.plot.pmod=20; p.nc.neig=50; 
[p.u,res]=nloop(p,p.u);fprintf('first res=%g\n',res); 
%% some more settings - setup of illupack
p.file.smod=1; p.sw.bifcheck=0; p.sw.spcalc=1;
p.sw.verb=2; p=setbelilup(p,0,1e-3,5,1e-4,500); p.sol.ds = 0.01;
%% Continuatation
p.nc.lammax = 0.7;
p.sw.bifcheck=2; p.sw.spcalc=1;
tic; p=cont(p,100); toc;
%% Find bifurcations
p.sw.bifcheck=2; p.sw.spcalc=1;
p.nc.dsmin=0.0;
p.nc.nsteps=5000;p.nc.lammax = 1.0;
tic; p=findbif(p,1); toc
%%
p.sw.bifcheck=0; p.sw.spcalc=1;
p.nc.lammax = 1.0; 
tic; p=cont(p,100); toc;
%% test bpt
p=swibra('2D_b2only_bs','bpt1','tmp',-0.01); p.sw.bifcheck=0; p.sw.spcalc=1; 
p.plot.pmod = 1;
tic; p=cont(p,10); toc
%% bifurcation to subcritical hom b2 branch
p=swibra('2D_b2only_bs','bpt1','2D_b2only_hom',0.01); p.sw.bifcheck=0; p.sw.spcalc=1; 
p.plot.pmod = 1;
tic; p=cont(p,10); toc
%% follow back hom b2 branch to findbif last bifurcation point 
p = loadp('2D_b2only_hom','pt30','2D_b2only_hom-b');
%% continuation to backward on the mix hom branch
p=resetc(p); p.fuha.savefu(p); 
p.nc.lammax=1.1; p.sol.ds=-p.sol.ds;
p.sw.bifcheck=2; p.sw.spcalc=1;
p.nc.dsmin=0.0;
p.nc.nsteps=5000;
tic; p=findbif(p,10); toc
%% Test bifurcations from mix hom branch
p=mswibra('2D_b2only_hom-b','bpt1','tmp');
p=taugen(p,[-1 1]);
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=cont(p,22); toc
%% Gap patterns
p=mswibra('2D_b2only_hom-b','bpt1','2D_b2only_hex1');
p=taugen(p,[-1 -1]);
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=cont(p,17); toc
%% Stripes
p=mswibra('2D_b2only_hom-b','bpt1','2D_b2only_str1');
p=taugen(p,[0 1]);
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=cont(p,25); toc
%% Dots patterns
p=mswibra('2D_b2only_hom-b','bpt1','2D_b2only_hex2');
p=taugen(p,[-1 1]);
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=cont(p,22); toc
%% Bifurcation to gap b2 and dot b1
p=mswibra('2Dmix-b','bpt1','2Dmix-hex1');
p=taugen(p,[-1 -1]);
plotaf2(p,7);pause;
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=contaf2D(p,50); toc
%% Following stripes
p=mswibra('2Dmix-b','bpt1','2Dmix-str1');
p=taugen(p,[1 0]);
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=contaf2D(p,50); toc
%% Following negative hexagons
p=mswibra('2Dmix-b','bpt1','2Dmix-hex2');
p=taugen(p,[1 -1]);
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=contaf2D(p,50); toc
%% Creating pure b2 state from the mixed branch
p = loadp('2Dmix-b','pt0','2Db2');
n=p.np; % load some SS as IC 
po=getpte(p); x=po(1,:)'; x1=min(x); x2=max(x);  
p.u(1:n)=0*p.u(1:n);
[p.u,res]=nloop(p,p.u);fprintf('first res=%g\n',res); 
p=resetc(p); p.fuha.savefu(p); 
p.sw.bifcheck=0; p.sw.spcalc=1;
p.nc.lammax=1.1;p.sol.ds=-p.sol.ds;
tic; p=contaf2D(p,50); toc
%%
p = loadp('2Dmix-hex2','pt22','2Dmix-hex2-1');
p.fuha.lss=@lss; p.fuha.blss=@lss;
%%
p.sw.bifcheck=0; p.sw.spcalc=1;
tic; p=pmcont(p,10); toc
%% non-hom bifurcation: 1
p=swibra('2Dmix-b','bpt1','2Dmix-pat1',-0.01); p.sw.bifcheck=0; p.sw.spcalc=2; 
p.plot.pmod = 1;
tic; p=contaf2D(p,20); toc
%% non-hom bifurcation: 5
p=swibra('2Dmix-b','bpt5','2Dmix-pat2',-0.01); p.sw.bifcheck=0; p.sw.spcalc=2; 
p.plot.pmod = 1;
tic; p=contaf2D(p,20); toc
%% Trying to converge back to stable patterned state
p=loadp('1Dmix-pat2','pt36','tint');
plotsol('1Dmix-pat2','pt36',7,[1 2]);
%%
n=p.np; % load some SS as IC 
po=getpte(p); x=po(1,:)'; x1=min(x); x2=max(x); 
per=0.1*rand(p.np,1);
p.u(1:n)=p.u(1:n)+0.1*(1+rand(p.np,1));
p.u(n+1:2*n)=p.u(n+1:2*n)+0.1*(1+rand(p.np,1)); plotsol(p); % perturb IC 
%% 
t0=0; ts=[]; dt=0.1; nt=10000; nc=0; pmod=1000; smod=1000; nffu=@afrhs; p.mat.Kadv=0; 
%% time integration, repeat this cell until convergence 
[p,t0,ts,nc]=tintxs(p,t0,ts,dt,10*nt,nc,pmod,smod,nffu);
plotaf2(p,7);
%% 
clf(2); p=setfn(p,'2Dmix-pat3'); p=resetc(p); p.fuha.savefu(p); 
p.sw.bifcheck=0; p.sol.restart=1; %  p.sol.ds=-p.sol.ds;
tic; p=contaf2D(p,50); toc
%%
p = loadp('2Dmix-pat3','pt0','2Dmix-pat3-f');p=resetc(p); p.fuha.savefu(p); 
p.sw.bifcheck=0; p.sol.restart=1; p.sol.ds=-p.sol.ds;
tic; p=contaf2D(p,50); toc
%% Generating the b2 only hex1 pattern
p=loadp('2Dmix-str1','pt25','tint');
n=p.np; % load some SS as IC 
po=getpte(p); x=po(1,:)'; x1=min(x); x2=max(x); 
%per=0.1*rand(p.np,1); 
p.u(1:n)=0*p.u(1:n); %0.1*(1+rand(p.np,1)); plotsol(p); % perturb IC 
plotaf2(p,7); % perturb IC 
p.plot.pcmp=2;
%% 
t0=0; ts=[]; dt=0.1; nt=10000; nc=0; pmod=1000; smod=1000; nffu=@afrhs; p.mat.Kadv=0; 
%% time integration, repeat this cell until convergence 
[p,t0,ts,nc]=tintxs(p,t0,ts,0.1*dt,10*nt,nc,pmod,smod,nffu);
plotaf2(p,7);
%% Continuation of patterned state after convergence till bpt5
clf(2); p=setfn(p,'2Db2-hex1'); p=resetc(p); p.fuha.savefu(p); 
p.sw.bifcheck=0; p.sol.restart=1; %  p.sol.ds=-p.sol.ds;
tic; p=contaf2D(p,50); toc
%% backward continuation with findbif to find last bifurcation point before it's stable
p=loadp('2Dhom-mix','pt210','2Dhom-mix-1');    
p=setbelilup(p,0,1e-3,5,1e-4,500);             
p.file.smod=10;p.sw.bifcheck=2; p.sw.spcalc=1; 
p.sol.ds=0.001;                                
p.nc.dsmin=0.0;                                
tic; p=findbif(p,100); toc;   
%% bifurcation to mix hom state
p=swibra('2Dhom-b1','bpt5','2Dhom-mix', -0.01); 
p=setbelilup(p,0,1e-3,5,1e-4,500); 
p.file.smod=10;p.sw.bifcheck=2; p.sw.spcalc=1; 
tic; p=cont(p,100); toc;
%% test a bifurcation point
p=swibra('2Dbs','bpt10','tmp', -0.01); 
p.file.smod=1;p.plot.pmod=1;p.sw.bifcheck=0; p.sw.spcalc=1; 
tic; p=contaf2D(p,5); toc
plotaf2(p,7);
%% Mesh adaptation
p.fuha.e2rs=@e2rs; p=meshada(p,'ngen',5,'sig',0.3); % Repeat two times
%% Trying to converge to new patterned state by perturbing unstable patterned one
p=loadp('2Dpat-hex4','pt20','tint');
n=p.np; % load some SS as IC 
po=getpte(p); x=po(1,:)'; x1=min(x); x2=max(x); 
%per=0.1*rand(p.np,1); 
p.u(n+1:2*n)=p.u(n+1:2*n)+per;p.u(1:n)=p.u(1:n)+per; %0.1*(1+rand(p.np,1)); plotsol(p); % perturb IC 
plotaf(p,7,[1 2]); % perturb IC 
%% 
t0=0; ts=[]; dt=0.01; nt=10000; nc=0; pmod=1000; smod=1000; nffu=@afrhs; p.mat.Kadv=0; 
%% time integration, repeat this cell until convergence 
[p,t0,ts,nc]=tintxs(p,t0,ts,dt,10*nt,nc,pmod,smod,nffu);
%% Continuation of patterned state after convergence till bpt5
clf(2); p=setfn(p,'2Db2-hex1'); p=resetc(p); p.fuha.savefu(p); 
p.sw.bifcheck=0; p.sol.restart=1;   p.sol.ds=-0.01;
%%
p=cont(p,50);
%% plot
pt='pt100';
figure(7); clf;figure(8); clf;
plotsol('2Dpat-mix1',pt,7,1,2);plotsol('2Dpat-mix1',pt,8,2,2);
%% follow hom.branch b1 and b2
p=swibra('hom-b1','bpt8','hom-mix', -0.01); p.sw.bifcheck=2; p.sw.spcalc=1;
p.nc.neig=40;
tic; p=cont(p,200); toc
%% follow hom.branch b1 and b2
p=swibra('hom-mix','bpt22','tmp', -0.01); p.sw.bifcheck=2; p.sw.spcalc=1;
p.nc.neig=40;
tic; p=pmcont(p,20); toc