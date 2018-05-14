%%
% STANUFU: standard "user function" called after each cont.step
%
%  [p,cstop]=stanufu(p,brout,ds)
% brout=branch data from internal calculations
% ds   =current stepsize
% Returns cstop flag for error post-processing
%
% See also cont, stanparam
function [p,cstop]=spufu(p,brout,ds)
p=ulamcheck(p); % test if a desired lambda value has been passed
brplot=brout(length(bradat(p))+p.plot.bpcmp); %y-axis value in bif-figure
fprintf('%4i %s %s %5.2e %4i  %s %s ', ...
    p.file.count,printcon(getlam(p)),printcon(brplot),p.sol.res,p.sol.iter, ...
    p.sol.meth,printcon(ds));
if(p.sw.errcheck>0) 
    fprintf('%5.2e ',p.sol.err);
end
if(p.sw.spcalc==1) 
    fprintf(' %2i ', p.sol.ineg); 
end
npr=length(p.sw.bprint);
for i=1:npr; fprintf('%s ',printcon(brout(length(bradat(p))+p.sw.bprint(i)))); end;
% put anything else here
fprintf('\n');
cstop=0;
if(getlam(p)<p.nc.lammin)
    fprintf('  lam=%g < lammin=%g, stopping\n',getlam(p),p.nc.lammin); cstop=1;
end
if(getlam(p)>p.nc.lammax)
    fprintf('  lam=%g > lammax=%g, stopping\n',getlam(p),p.nc.lammax); cstop=1;
end
n=p.np; par=p.u(p.nu+1:end); u=p.u; u=[u(1); u(n+1); u(2*n+1)];
uv=[u;par];% u
J=bwh_sp_jac(p,uv);
dw=par(15); dh=par(16); kv=0:0.001:1; kl=length(kv); lamv=zeros(3,kl); 
for i=1:kl 
    k=kv(i); 
    K=[[k^2 0 0];[0 dw*k^2 0]; [0 0 dh*k^2]]; 
    A=J-K; 
    lam=eig(A); 
    [~, ix]=sort(real(lam)); 
    for j=1:3; lamv(j,i)=lam(ix(j)); end 
end 
figure(10); clf; plot(kv, real(lamv(3,:)), kv, imag(lamv(3,:)));  
%hold on; plot(kv, real(lamv(2,:)), kv, imag(lamv(2,:))); 
pause 
end
