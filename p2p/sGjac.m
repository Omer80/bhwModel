function Gu=sGjac(p,u)
%% w
[f1b2,f1w1,f1w2,f2b2,f2w1,f2w2,f3b2,f3w1,f3w2]=bwh_jac(p,u);
%% 
n=p.np; sc=u(p.nu+21);
Fu=[[spdiags(f1b2,0,n,n),spdiags(f1w1,0,n,n),spdiags(f1w2,0,n,n)];
    [spdiags(f2b2,0,n,n),spdiags(f2w1,0,n,n),spdiags(f2w2,0,n,n)];
    [spdiags(f3b2,0,n,n),spdiags(f3w1,0,n,n),spdiags(f3w2,0,n,n)]];
Gu=sc*p.mat.K-p.mat.M*Fu; 
end 