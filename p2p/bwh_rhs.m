function f=bwh_rhs(p,u) 
% rhs of nondimensional BWH model with trade-off parameter
B=u(1:p.np); W=u(p.np+1:2*p.np);H=u(2*p.np+1:3*p.np);
par=u(p.nu+1:end);
% Loading to parameters names
prec=par(1);eta=par(2);nuw=par(3);rhow=par(4);nuh=par(5);
rhoh=par(6);q=par(7);ff=par(8);chi=par(9);gam=par(10);alpha=par(11);

del_to = 0.3;

q_min   = q*(1.0-del_to);
q_max   = q*(1.0+del_to);
eta_min = eta*(1.0-del_to);
eta_max = eta*(1.0+del_to);
K_min   = (1.0-del_to);
K_max   = (1.0+del_to);

K_to   = K_max    + chi*(K_min-K_max);
q_to   = (q_max   + chi*(q_min-q_max))/K_to;
eta_to = (eta_max + (1-chi)*(eta_min-eta_max))*K_to;
gam    = gam*K_to;

G = W.*(1 + eta_to*B).^2;
I = alpha*((B + q_to*ff)./(B + q_to));
evapw = ((nuw)./(1 + rhow*B)).*W;
evaph = ((nuh)./(1 + rhoh*B)).*H;
tras  = gam*B.*G;

f1=G.*B.*(1-B)-B;
f2=I.*H-evapw-tras;
f3=prec-I.*H-evaph;
f=[f1;f2;f3]; 
end