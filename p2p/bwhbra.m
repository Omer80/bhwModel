%%
% STANBRA: standard function for output beyond bradat(p) to p.branch.
%
%  out=stanbra(p,u)
%
% Here as template forusers: out=[ aux vars, max|u_1|, min|u_1| ]
% Hence, for k<= # aux var: bpcmp=k means aux var k (usually parameter k).
%
% See also stanparam.
function out=bwhbra(p,u)
%u=p.mat.fill*u(1:p.nu); 
n=p.np;
B=u(1:p.np); W=u(p.np+1:2*p.np);H=u(2*p.np+1:3*p.np);
par=u(p.nu+1:end);
% Loading to parameters names
prec=par(1);

Bnorm=p.mat.M(1:n, 1:n)*B;
Bcoverage = numel(find(Bnorm>max(Bnorm(:))/2))/numel(Bnorm);

out=[prec;% Need to add the continuation parameter to use auxdict
     Bcoverage;
     max(B); min(B);mean(B)];
end
