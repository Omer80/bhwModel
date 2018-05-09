#include "auto_f2c.h"
#include <math.h>
#include <stdio.h>
#define F2C -1

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*   The OZ model                      	                                  */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* System generated locals */
  //integer dfdu_dim1, dfdp_dim1;
  
  /* Local variables */
  doublereal B,W,H,Bx,Wx,Hx, X;
  doublereal dummy_b,dummy_w, dummy_j;
  doublereal eta, nuw, nuh, rhow, rhoh, gam, alpha, q,ff,p, delw, delh;
  doublereal pi, num_periods,L,powH;
  doublereal e1,e2,e3;
  pi = atan(1.) * 4.0;  
  /* Defining the parameters */
 
  eta   = par[1+F2C];
  nuw   = par[2+F2C];
  nuh   = par[3+F2C];
  rhow  = par[4+F2C];
  rhoh  = par[5+F2C];
  gam = par[6+F2C];
  alpha = par[7+F2C];
  q     = par[8+F2C];
  ff    = par[9+F2C];
  dummy_b = par[13+F2C]; // Advection terms for Hopf fixing
  dummy_w = par[14+F2C]; // Advection terms for Hopf fixing
  dummy_j = par[15+F2C]; // Advection terms for Hopf fixing
  p     = par[16+F2C];
  delw  = par[17+F2C];
  delh  = par[18+F2C];
  powH  = par[19+F2C];

  //num_periods = par[17+F2C];
  L     = 1; //    = par[18+F2C];

  // 2*pi/L = kf/2 for locking. We multiply by num_periods to make the domain larger
  
  /* Function Body */
  B      = u[0];
  W      = u[1];
  H      = u[2];  // J = H^2
  Bx     = u[3];  // x stands for space derivative
  Wx     = u[4];
  Hx     = u[5];
 

  // f[i] = u[i]_x
  // We multiply by L because the derivatives are relative to "AUTO"s space.
  // x_real = [0,L], and x_auto = [0,1]. Therefore x_real = L * x_auto
  // d/dx_auto = L * d/dx_real
  f[0] = L * Bx;
  f[1] = L * Wx;
  f[2] = L * Hx;
  powH=2;
  f[3] = -L * (W*pow((1 + eta*B),2)*B*(1-B) - B + dummy_b * Bx) ;//+ dummy_b * Bx;
  f[4] = -L * (1/delw)*((alpha*((B + q*ff)/(B + q)))*H - ((nuw)/(1 + rhow*B))*W - (gam*B*pow((1 + eta*B),2))*W + (dummy_b + dummy_w)  * Wx) ;//+ (dummy_b + dummy_w)  * Wx;
  f[5] = -L * (1/(delh * powH * pow(H, (powH - 1))))*(p - (alpha*((B + q*ff)/(B + q))) * H - ((nuh)/(1 + rhoh*B)) * H + (delh * powH * (powH - 1) * pow(H, (powH - 2))* Hx*Hx)  + (dummy_b + dummy_j) * Hx);// + (dummy_b + dummy_j) * Hx;
  //f[4] = -L * (1/delw)*((alpha*((B + q*ff)/(B + q)))*H - ((nuw)*(1 - rhow*B))*W - (gam*B*pow((1 + eta*B),2))*W + (dummy_b + dummy_w)  * Wx) ;//+ (dummy_b + dummy_w)  * Wx;
  //f[5] = -L * (1/(delh * powH * pow(H, (powH - 1))))*(p - (alpha*((B + q*ff)/(B + q))) * H - ((nuh)*(1 - rhoh*B)) * H + (delh * powH * (powH - 1) * pow(H, (powH - 2))* Hx*Hx)  + (dummy_b + dummy_j) * Hx);// + (dummy_b + dummy_j) * Hx;
  //f[5] = -L * (1/(delh))*(p - (alpha*((B + q*ff)/(B + q))) * H - ((nuh)/(1 + rhoh*B)) * H  ) ;//+ (dummy_b + dummy_j) * Hx;
 
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal x,
           doublereal *u, doublereal *par)
{
  /* local variables */
  doublereal B,W,H,Bx,Wx,Hx,X;
  doublereal dummy_b,dummy_w, dummy_j;
  doublereal E, K, M, NW, NH, Lambda, Gamma, RW, RH, DeltaB, DeltaW, N, A, Q, DeltaH; // The fixed dimensional parameters
  doublereal eta, nuw, nuh, rhow, rhoh, gam, alpha, q,ff,p, delw, delh, Lb,FF;
  doublereal L, amp, pi, num_periods, dm, dp,offset,tet,batch,tat,powH;
  doublereal e1,e2,e3;
  
  /* defining the numerical values of the dimensional parameters */
  // Gilad's model parameters
  //E      = 3.5;  // m^2 / kg Root’s augmentation per unit biomass
  //K      = 1.0;  // kg / m^2 Maximum standing biomass
  //M      = 1.2; // 1 / yr Rate of biomass loss due to mortality
  //NW     = 4.0; // 1 / yr  Soil water evaporation rate
  //NH     = 6.0; // 1 / yr  Surface water evaporation rate
  //Lambda = 0.032;  // (m^2 / kg) / yr Biomass growth rate
  //Gamma  = 20;   // (m^2 / kg) / yr Soil water consumption rate
  //RW      = 0.95;  //  Soil water evaporation reduction due to shading
  //RH      = 1.05;  //  Soil water evaporation reduction due to shading
  //DeltaB = 0.000625;  // m^2 / yr Seed dispersal coefficient
  //DeltaW = 0.0625;  // m^2 / yr Transport coefficient for soil water
  //DeltaH = 0.05; // m^2 / yr (kg / m^2)^{-1} Bottom friction coefficient between surface water and ground surface
  //Q      = 0.05; // kg / m^2 Biomass reference value beyond which infiltration rate under a patch approaches its maximum
  //A      = 40.0; // 1/yr  Infiltration rate in fully vegetated soil
  

  // Yuval's model parameters
  //E      = 7.0;  // m^2 / kg Root’s augmentation per unit biomass
  //K      = 0.4;  // kg / m^2 Maximum standing biomass
  //M      = 10.5; // 1 / yr Rate of biomass loss due to mortality
  //NW     = 15.0; // 1 / yr  Soil water evaporation rate
  //NH     = 16.0; // 1 / yr  Surface water evaporation rate
  //Lambda = 0.9;  // (m^2 / kg) / yr Biomass growth rate
  //Gamma  = 12.0;   // (m^2 / kg) / yr Soil water consumption rate
  //RW      = 0.7;  //  Soil water evaporation reduction due to shading
  //RH      = 1.05;  //  Soil water evaporation reduction due to shading
  //DeltaB = 1.2;  // m^2 / yr Seed dispersal coefficient
  //DeltaW = 150.0;  // m^2 / yr Transport coefficient for soil water
  //DeltaH = 0.05; // m^2 / yr (kg / m^2)^{-1} Bottom friction coefficient between surface water and ground surface
  //Q      = 0.05; // kg / m^2 Biomass reference value beyond which infiltration rate under a patch approaches its maximum
  //A      = 40.0; // 1/yr  Infiltration rate in fully vegetated soil

  // Modified Yuval's model parameters
  //E      = 3.0;  // m^2 / kg Root’s augmentation per unit biomass
  //K      = 0.9;  // kg / m^2 Maximum standing biomass
  //M      = 5.5; // 1 / yr Rate of biomass loss due to mortality
  //NW     = 15.0; // 1 / yr  Soil water evaporation rate
  //NH     = 16.0; // 1 / yr  Surface water evaporation rate
  //Lambda = 0.9;  // (m^2 / kg) / yr Biomass growth rate
  //Gamma  = 12.0;   // (m^2 / kg) / yr Soil water consumption rate
  //RW      = 0.7;  //  Soil water evaporation reduction due to shading
  //RH      = 1.05;  //  Soil water evaporation reduction due to shading
  //DeltaB = 1.2;  // m^2 / yr Seed dispersal coefficient
  //DeltaW = 150.0;  // m^2 / yr Transport coefficient for soil water
  //DeltaH = 0.05; // m^2 / yr (kg / m^2)^{-1} Bottom friction coefficient between surface water and ground surface
  //Q      = 0.05; // kg / m^2 Biomass reference value beyond which infiltration rate under a patch approaches its maximum
  //A      = 40.0; // 1/yr  Infiltration rate in fully vegetated soil
 
  // Dimensional parameters based on Hezi's simulation on Apr 12
  //E      = 3.0;  // m^2 / kg Root’s augmentation per unit biomass
  //K      = 0.33333;  // kg / m^2 Maximum standing biomass
  //M      = 3.0; // 1 / yr Rate of biomass loss due to mortality
  //NW     = 6.0; // 1 / yr  Soil water evaporation rate
  //NH     = 4.2; // 1 / yr  Surface water evaporation rate
  //Lambda = 0.3;  // (m^2 / kg) / yr Biomass growth rate
  //Gamma  = 14.0;   // (m^2 / kg) / yr Soil water consumption rate
  //RW      = 1.0;  //  Soil water evaporation reduction due to shading
  //RH      = 1.0;  //  Soil water evaporation reduction due to shading
  //DeltaB = 0.00025;  // m^2 / yr Seed dispersal coefficient
  //DeltaW = 0.00025;  // m^2 / yr Transport coefficient for soil water
  //DeltaH = 0.05; // m^2 / yr (kg / m^2)^{-1} Bottom friction coefficient between surface water and ground surface
  //Q      = 0.433; // kg / m^2 Biomass reference value beyond which infiltration rate under a patch approaches its maximum
  //A      = 240.0; // 1/yr  Infiltration rate in fully vegetated soil
  
  // Modified Dimensional parameters based on Hezi's simulation on Apr 12
  // REAR_PARAMETERS_FROM_HERE
  E      = 1.5;  // m^2 / kg Root’s augmentation per unit biomass
  K      = 2.0/3.0;  // kg / m^2 Maximum standing biomass
  M      = 2.0; // 1 / yr Rate of biomass loss due to mortality
  NW     = 1.5; // 1 / yr  Soil water evaporation rate
  NH     = 4.5; // 1 / yr  Surface water evaporation rate
  Lambda = 0.03;  // (m^2 / kg) / yr Biomass growth rate
  Gamma  = 14.0;   // (m^2 / kg) / yr Soil water consumption rate
  RW      = 0.3;  //  Soil water evaporation reduction due to shading
  RH      = 0.8;  //  Soil water evaporation reduction due to shading
  DeltaB = 0.1;  // m^2 / yr Seed dispersal coefficient
  DeltaW = 2.5;  // m^2 / yr Transport coefficient for soil water
  DeltaH = 4.0; // m^2 / yr (kg / m^2)^{-1} Bottom friction coefficient between surface water and ground surface
  Q      = 1.2; // kg / m^2 Biomass reference value beyond which infiltration rate under a patch approaches its maximum
  A      = 120.0; // 1/yr  Infiltration rate in fully vegetated soil
  FF     = 0.01;
  
  printf ( "Critical percipitation: P_c = %4.2f ",  (((M*NW)/Lambda)*(1 + (NH/(A*FF)))));


  /*  defining the numerical value of the non-dimensional parameters */

  // Calculating from the dimensional parameters
  eta   = E*K;
  nuw   = NW / M;
  nuh   = NH / M;
  rhow  = RW;
  rhoh  = RH;
  gam   = (Gamma * K) / M;
  alpha = A / M;
  q     = Q / K;
  ff    = FF;
  //p   = 1.9;
  delw  = DeltaW / DeltaB;
  delh  = (DeltaH * M) / (DeltaB * Lambda);
  powH  = 2;
  p = nuw * ( (alpha * ff) + nuh)/ (alpha * ff) - 0.1; // Critical value of p for stable bare-soil, minus 1
  
  printf ( "Critical percipitation: p_c = %4.2f \n",  (p + 0.1) );
  printf ( "eta=%4.2f,gam=%4.2f,alpha=%4.2f,ff=%4.2f,q=%4.2f, \n",  eta,gam,alpha,ff,q);
  
  // Direct input of non-dimensional parameters
  
  //eta   = 1.0;
  //nuw   = 2.0 ;
  //nuh   = 4.4 ;
  //rhow  = 1.0;
  //rhoh  = 4.0;
  //gam   = 14.0 ;
  //alpha = 18.0 ;
  //q     = 0.05;
  //ff    = 0.01;
  ////p   = 1.9;
  //delw  = 10.0  ;
  //delh  = 20000.0 ;
  //powH  = 2;
  //p = nuw * ( (alpha * ff) + nuh)/ (alpha * ff) - 0.01;
  
  
  pi = atan(1.) * 4.0;  
  num_periods = 26;
  L       = 320;
  
  dummy_b = (doublereal).00000000000000002;
  dummy_w = (doublereal).00000000000000003;
  dummy_j = (doublereal).00000000000000005;

  /* load into internal parameters */

 
  par[1+F2C]  = eta;
  par[2+F2C]  = nuw;
  par[3+F2C]  = nuh;
  par[4+F2C]  = rhow;
  par[5+F2C]  = rhoh;
  par[6+F2C]  = gam;
  par[7+F2C]  = alpha;
  par[8+F2C]  = q;
  par[9+F2C]  = ff;
  par[13+F2C] = dummy_b;
  par[14+F2C] = dummy_w;
  par[15+F2C] = dummy_j;
  par[16+F2C] = p;
  par[17+F2C] = delw;
  par[18+F2C] = delh;
  par[19+F2C] = powH;


  //par[17+F2C] = num_periods;
  par[11+F2C] = L;
  X     =  x*L;

  //par[12] = num_periods;
  //par[10] = L;
 
  // The exact Solution
  // the derivatives are relative to the "real" space: L * d/dx_real = d/dx_auto. That's why we divide by L.
  // The L multiplying and the L dividing cancel out, but we'd rather write it this way so we can always remember...
  //dp = 1/((kf+k0)*(kf+k0)-k0*k0)/((kf+k0)*(kf+k0)-k0*k0);
  //dm = 1/((kf-k0)*(kf-k0)-k0*k0)/((kf-k0)*(kf-k0)-k0*k0);
  //amp   =  2.0*sqrt(epsilon + gamma*gamma*epsilon*(dp+dm)/4.0)/sqrt(3.0)/1.0;
  //amp = 1.0*sqrt(epsilon/2.0)/sqrt(3.0);
  //amp = 0.0289334;
  //  amp   =  2.0*sqrt(-epsilon )/sqrt(3.0);
  
  tet=2*pi/L*num_periods;
  amp=1.02;
  //offset=0.021;
  batch=8.5*L/num_periods;
  tat=0.5;

  tet=2*pi/L*num_periods;
  amp=1.02;
  
  // Bare-soil stable state initial conditions
  B     =  0;
  W     =  ((alpha * ff)/nuw)*(p / (alpha * ff + nuh));
  H     =  (p / (alpha * ff + nuh));
  Bx    =  0;
  Wx    =  0;
  Hx    =  0;
 

  
  u[0] = B;  
  u[1] = W;  
  u[2] = H;  
  u[3] = Bx;
  u[4] = Wx;
  u[5] = Hx;
 
 

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp, 
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac, 
          doublereal *fb, doublereal *dbc) 
{ 
  fb[0] = u1[3];   // Bx_right   = 0 
  fb[1] = u0[3];   // Bx_left    = 0 
  fb[2] = u1[4];   // Wx_right   = 0 
  fb[3] = u0[4];   // Wx_left    = 0 
  fb[4] = u1[5];   // Jx_right   = 0 
  fb[5] = u0[5];   // Jx_left    = 0 
 
  return 0; 
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
    return 0;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

