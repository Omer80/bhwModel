function [f1B,f1W,f1H,f2B,f2W,f2H,f3B,f3W,f3H]=bwh_jac(p,u)
% jacobian of nondimensional AgroForestry model with NO light competition
B=u(1:p.np); W=u(p.np+1:2*p.np);H=u(2*p.np+1:3*p.np);
par=u(p.nu+1:end);
% Loading to parameters names
eta=par(2);nuw=par(3);rhow=par(4);nuh=par(5);
rhoh=par(6);q=par(7);ff=par(8);chi=par(9);gam=par(10);alpha=par(11);


f1B= - B.*W.*(B.*((3*chi)/5 - 13/10).*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1).^2 - W.*(B - 1).*(B.*((3*chi)/5 - 13/10).*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1).^2 - 2*B.*W.*((3*chi)/5 - 13/10).*(B - 1)*((13*eta)/10 + (3*eta*(chi - 1))/5).*(B.*((3*chi)/5 - 13/10)*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1) - 1;
f1W= -B.*(B - 1).*(B.*((3*chi)/5 - 13/10)*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1).^2;
f1H= zv;
f2B= (H.*alpha)./(B - ((13*q)/10 - (3*chi*q)/5)/((3*chi)/5 - 13/10)) + (W.*nuw*rhow)./(B.*rhow + 1).^2 + W.*gam.*((3*chi)/5 - 13/10).*(B.*((3*chi)/5 - 13/10).*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1).^2 - (H.*alpha.*(B - (ff*((13*q)/10 - (3*chi*q)/5))/((3*chi)/5 - 13/10)))./(B - ((13*q)/10 - (3*chi*q)/5)/((3*chi)/5 - 13/10))^2 + 2*B.*W.*gam.*(((3*chi)/5 - 13/10).^2).*((13*eta)/10 + (3*eta*(chi - 1))/5).*(B.*((3*chi)/5 - 13/10).*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1);
f2W= B.*gam.*((3*chi)/5 - 13/10).*(B.*((3*chi)/5 - 13/10)*((13*eta)/10 + (3*eta*(chi - 1))/5) - 1).^2 - nuw./(B.*rhow + 1);
f2H= (alpha*(B - (ff*((13*q)/10 - (3*chi*q)/5))./((3*chi)/5 - 13/10)))./(B - ((13*q)/10 - (3*chi*q)/5)/((3*chi)/5 - 13/10));
f3B= (H.*nuh.*rhoh)./(B.*rhoh + 1)^2 - (H.*alpha)./(B - ((13*q)/10 - (3*chi*q)/5)/((3*chi)/5 - 13/10)) + (H.*alpha.*(B - (ff*((13*q)/10 - (3*chi*q)/5))/((3*chi)/5 - 13/10)))./(B - ((13*q)/10 - (3*chi*q)/5)/((3*chi)/5 - 13/10)).^2;
f3W= zv;
f3H= - nuh./(B.*rhoh + 1) - (alpha.*(B - (ff*((13*q)/10 - (3*chi*q)/5))./((3*chi)/5 - 13/10)))./(B - ((13*q)/10 - (3*chi*q)/5)/((3*chi)/5 - 13/10));
end