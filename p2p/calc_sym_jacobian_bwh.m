close all; format compact; keep pphome;clc;
%% 
syms prec lamb eta gam rhow rhoh nuw nuh q ff alpha chi
syms B W H
%% Copy the system equations to this block
del_to = 0.3;

q_min   = q*(1.0-del_to);
q_max   = q*(1.0+del_to);
eta_min = eta*(1.0-del_to);
eta_max = eta*(1.0+del_to);
K_min   = (1.0-del_to);
K_max   = (1.0+del_to);

K_to   = K_max    + chi*(K_min-K_max);
q_to   = (q_max    + chi*(q_min-q_max))/K_to;
eta_to = (eta_max  + (1-chi)*(eta_min-eta_max))*K_to;
gam_to    = gam*K_to;

G = W*(1 + eta_to*B)*(1 + eta_to*B);
I = (alpha*((B + q_to*ff)/(B + q_to)));
evapw = ((nuw)/(1 + rhow*B))*W;
evaph = ((nuh)/(1 + rhoh*B))*H;
tras  = gam_to*B*G;

f1=G*B*(1-B)-B;
f2=I*H-evapw-tras;
f3=prec-I*H-evaph; 
%% Varifying the equations
disp(f1);
disp(f2);
disp(f3);
%% Symbolic calculation of the Jacobian
jac = jacobian([f1 f2 f3],[B W H]);
%% Display the components of the Jacobian
disp(jac(3,3));