function par = loadparms( parameters_fname )
% setting parameters

p = load(parameters_fname);
par(1)=p.p;
par(2)=p.eta;
par(3)=p.nuw;
par(4)=p.rhow;
par(5)=p.nuh;
par(6)=p.rhoh;
par(7)=p.q;
par(8)=p.f;
par(9)=p.chi;
par(10)=p.gamma;
par(11)=p.alpha;
par(15)=p.dw;
par(16)=p.dh;
par(21)=1.0; % Scale parameter that multiplies all of the diffusion coeffs
end

