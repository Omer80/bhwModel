%%
% STANPARAM: settings of p-structure to 'standard' values. 
% We set some fields to [], because this file is also intended 
% as a quick reference for the meaning of fields. 
%
%  p=stanparam(p)
function p=bwh_stanparam(p)
%%%%%%%% function handles for which there are 'standard choices' 
p.fuha.headfu=@stanheadfu; % headline 
p.fuha.ufu=@stanufu;    % printout, stop if lam<>p.sw.lammin,p.sw.lammax
p.fuha.outfu=@bwhbra;  % branch output 
p.fuha.lss=@lss;        % linear systems solver, default use \ 
p.fuha.blss=@lss;      % bordered linear system solver, default \
p.fuha.savefu=@stansavefu; % save function
p.fuha.postmmod=@stanpostmeshmod; % post-mesh-modification 
p.fuha.e2rs=@stane2rs;  % elements2refine selection function 
%%%%%%%% function handles which must/should be user provided     
% sfem= 0  setting: [c,f,a,b]=p.fuha.G(p,u)  (rhs)  
%                   [c,a,b]=p.fuha.Gjac(p,u) (for analytical Jacobians) 
% sfem=+-1 setting: r=p.fuha.sG(p,u), 
%                   Gu=p.fuha.sGjac(p,u,r)
%                   Guuphi=p.fuha.spjac(p,u) (spectral cont) 
%%%%%%%% numerical param   
p.nc.neq=4;             % #of PDEs (components in PDE system)
p.nc.tol=1e-8;          % tol for corrector 
p.nc.imax=10;           % #of newton-iterations in corrector 
p.nc.dsmin=0.0001; p.nc.dsmax=5; % min/max stepsize 
p.nc.lammin=-1e6; p.nc.lammax=1.5; % Bif.diagram bounds 
p.nc.nsteps=10;         % number of continuation steps (multi-steps for pmcont)
p.nc.ntot=1e4;          % max total number of continuation steps 
p.nc.del=1e-8;          % perturbation size for finite differences
p.nc.lamdtol=0.5;       % parametrization switch tolerance (for p.sw.para=1)
p.nc.dsinciter=p.nc.imax/2;   % increase ds if iter < dsinciter 
p.nc.dsincfac=2;        % by this factor 
p.nc.errbound=0;        % if >0 and errchecksw>0 this bound is determines call to mesh refinement
p.nc.dlammax=1;         % max difference in lambda
p.nc.neig=10;           % # eigenvalues for eigs for G_u, must be extended 
                        % to vector if  p.nc.eigref is vector 
p.nc.neigdet=0;         % # eigenvalues for eigs for A, 0 for LU-method
p.nc.eigref=0;          % 0 (initially), but extended to vector by, e.g., 
                        % initeig if looking for Hopf bifurcations                           
p.nc.eigint=[-1 0.1];   % 'sarn' eigenvalue region real part interval
p.nc.bisecmax=10;       % max # of bisections during special point locations 
p.nc.dsminbis=1e-9;     % dsmin for bif-and foldloc via bisection
p.nc.mu1=0.1;           % start bisec if ineg changed, and |re mu|<mu1
p.nc.mu2=1e-3;          % assume re(mu)=0 if |re(mu)|<mu2 at end of bisec 
p.nc.tol_ev=1e-4;       % threshold condition for instability flag 
%%%%%%%% LSS
p.bel.tol=1e-4;         % tol for lssbel (bordered elim)
p.bel.imax=10;          % max # of iterations in bel 
p.bel.bw=0;             % border width
p.ilup.droptol=1e-4;    % droptol for lssAMG        (may be adapted during lssAMG) 
p.ilup.droptolmin=1e-8; % min droptol for lssAMG
p.ilup.maxit=100;        % max # of GMRES iterations (may be adapted during lssAMG) 
p.ilup.maxitmax=1000;   % upper bound for max # of GMRES iterations
p.ilup.noconvhandling=0;% switch how to deal with no conv of AMG: 0: stop;  
                        % 1: change droptol, 2:  change maxit 
%%%%%%%% mesh-adaption 
p.nc.amod=0;            % adapt mesh each amod-th step (0 for never) 
p.nc.bddistx=0.1;       % boundary distance x for mesh-ref. for pBC
p.nc.bddisty=0.1; 
p.nc.ngen=5;            % max number of refinements during adaption 
p.nc.sig=0.5;           % sig-value for mesh-refinement
p.nc.posdef_tol = -1e-4; % threshold for negative solution continuation break
%%%%%%%% switches
p.sw.bifcheck=1;        % 0: off. 1: LU, 2: via spcalc
p.sw.cdbb=0; % continue directly behind bif: 0: use u1 computed in cont;                        
                        % 1: use ubb computed in bifdetec 
p.sw.bifloc=1;          % 0 for tangent, 1 for secant, 2 for quadratic in bif.localization
p.sw.foldcheck=0;       % 0: no fold detection, 1: fold detection on.
p.sw.spcont=0;          % 0=normal continuation, 2/1=fold/branch point (zero eigenvalue) continuation
p.sw.spcalc=1;          % 1/0 to calculate/not calculate stability EVals 
p.sw.spjac=1;           % 1/0 to use fuha.spjac or not during fold/branch point cont
p.sw.sfem=0;            % switch for simple Jacobian and residual, e.g. semilinear problem:
                        % 0: off, 1: use p.fuha.sG and p.fuha.sGjac, requires call to 
                        %            setfemops and setting p.eqn.c, p.eqn.b, p.eqn.a 
p.sw.para=1;            % parametrization switch 
p.sw.secpred=0;         % for BP continuation use secant predictor vi p.sw.secpred=1
p.sw.jac=1;             % Jacobian switch, % 1: Gu analytically, else FD
p.sw.qjac=1;            % 1: qu analytically, else FD 
p.sw.newt=0;            % 0=newton, 1=chord, 2=...  
p.sw.norm='inf';        % norm and tol for corrector 
p.sw.errcheck=0;        % 0: off, 1: put err-est to p.sol.err, but no further action 
                        % 2: meshadac if p.sol.err>p.nc.errbound, >2: as for 2 but refine mesh 
p.sw.eigmeth='eigs';    % eigendata computation method 'eigs' or 'sarn'
p.sw.evopts.disp=0;     % don't display anything during EVal calculations
p.sw.eigsstart=1;       % 0 to use random start for eigs, 1 for [1;...;1]  
p.sw.inter=1;p.sw.verb=1;   % interaction/verbosity switch: 0=none, 1=some, 2=much
p.sw.bprint=[];         % #branchcompos for printout by outfu
p.sw.bcper=0;           % switch for periodic BC 
%%%%%%%% plotting 
p.plot.pstyle=1;        % solution plot style: 1 for mesh-plot of u, 2 for contour-plot, 
                        % 3 for surface, with lightning, 4 for only the mesh
p.plot.pfig=1;          % screen layout of profile figure
p.plot.brfig=2;         % same for branch figure
p.plot.ifig=6;          % same for info figure
p.plot.pmod=1;          % plot every pmod-th step, 
p.plot.pcmp=1; p.plot.bpcmp=5; % component of sol. resp. branch for plotting 
p.plot.cm='cool';       % colormap 
p.plot.lpos=[0 0 10];   % light-pos for pstyle=3
p.plot.axis='tight';    % choose, e.g., 'tight', 'equal', 'image' 
p.plot.fs=16;           % fontsize for sol-plots, 
p.plot.labelsw=0;       % 1/0 for labels/no labels in solplot
p.plot.spfig=4;         % figure number for spectral output with spcalc (or specGu)
p.plot.brafig=3;        % figure number for standard post-computation plotting of branch
p.plot.spfig2=6;        % figure number for spectral output of bordered matrix Gua used in jaccheck
p.plot.udict={'b','w','h'};
p.plot.auxdict={'p','coverage','max(b)','min(b)','mean(b)'};
p.plot.alpha=1;       % transparency for 3D plots (small alpha=transparent) 
p.plot.ng=20;           % # points for interpol. for isosurf plot
p.plot.levc={'blue','red'}; % colors for isoplots
p.plot.fancybd=1;       % 0: old plotbra, 1 labels to point via line, 2 via annotate
p.plot.lsw=1;           % Switch for default reg/FP/HP/BP/usrlam labels, 0: all off!
% binary setting, i.e. usr=2^0, BP=2^1, ..., reg=2^4=16. eg: lsw=3=1+2 means usr and BP
% lsw=16 means: reg, lsw=18 means reg+BP, lsw=31 means all
%%%%%%%% file handling 
p.file.pnamesw=0;       % 1 to automatically set prefix of file names to variable name
p.file.mdir='meshes';   % dir so save meshes if ms=0; 
p.file.dirchecksw=0;    % 1 for user check if directories need to be created, 
                        % or written into if exist;  0 no check
p.file.smod=10;          % save ev. smod-th step, 0 for none 
p.file.count=0; p.file.bcount=1; p.file.hcount=1; p.file.fcount=1; % counter for step/branching/fold points 
p.file.msave=1;         % mesh saving: 1 save mesh in every point, 0 only save mesh if it is changed 
p.file.single=0;        % 0 save p.u, p.tau, p.sol.muv, and p.branch as double, 1 save as single
%%%%%%%% timing: here only for reference, values set in cont/pmcont 
p.time.timesw=1;        % 1 for output at end of cont
p.time.tot=0; p.time.totst=0; p.time.st=0; % times: total, total for steps, current step
p.time.bif=0; p.time.spec=0; p.time.newton=0; % times for: bifcheck, spcalc, nloop 
%%%%%%%% pmcont 
p.pm.resfac=0.2;        % resi-improvement for pmcont 
p.pm.mst=10;            % #predictors for pmcont 
p.pm.imax=1;            % base-max-iterations for pmcont 
%%%%%%%% fsolve 
p.fsol.fsol=0;          % use fsolve? 0:no, 1:pde, 2:ext, 3: both
p.fsol.meth=1;          % 1=trust-reg., else Levenberg-Marquardt
p.fsol.disp=3;          % 0=off, 1=final, 2=notify, 3=iter
p.fsol.tol=p.nc.tol^2; p.fsol.imax=5; 
p.fsol.opt=optimset('Jacobian','on','PrecondBandWidth',0);
%%%%%%%% runtime data, here initialized for startup 
p.sol.deta=0;           % determinant of linearization -- change is used to detect bif
p.sol.err=0;            % a posteriori error estimate
p.sol.meth=' ';         % string code for continuation step method 'arc', or 'nat'. See p.sw.para
p.sol.res=0;            % residual of solution in norm p.sw.norm
p.sol.iter=0;           % number of iterations used for this solution
p.sol.ineg=-1;          % number of negative (i.e. unstable) eigenvalues
p.sol.muv=[];           % eigenvalues (vector, or matrix if eigref is vector)
p.sol.lamd=0;           % step length in primary parameter (from tangent vector)
p.sol.restart=1;        % 1: make initial steps, 0: use tangent p.tau 
p.tau=1;                % tangent vector, here trivially initialized
p.branch=[];            % branch, here trivially initialized
p.u=[];                 % solution, here trivially initialized
p.sol.xiq=0;            % weight of auxiliary eqns for arclength scalar product (needs to exist)
p.mat.fill=1;           % fill operator (nontrivial in case of periodic b.c.)
p.mat.drop=1;           % drop operator (nontrivial in case of periodic b.c.)
p.mat.M=[]; p.mat.K=[]; % empty mass and stiffness matrices 
p.mesh=[];              % the FEM mesh
p.mesh.sympoi=0;        % if 1 then make poimesh symmetric
p.sol.ptype=-99;        % initial point type; types are: -1 (initial point), 
                        % -2 (swibra point), 0 (normal), 1 (branch point), 
                        % 2 (fold point), 3 (Hopf point), 4 (normal point on Hopf branch)
p.usrlam=[];            % desired lam-values
p.nc.nq=0;              % number of auxiliary equations
%%%% my extra parameter %%%%
p.sw.spbif=0;           % using spectral cont. for bifdetection if spbif==1
p.sw.mnorm=0;           % 0: use standart xi-norm, 1: use M-norm 
p.sol.mxiq=0;           % weight of auxiliary eqns for arclength scalar product when M-norm is used (needs to exist)
p.sol.mxi=0.5;          % weight for the norm of tau when M-norm is used