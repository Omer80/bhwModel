#include "auto_f2c.h"
#define F2C -1
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*   The BWH model with tradeoff parameter chi                            */
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
  doublereal p,lamb_max,lamb_min,eta,nuw,nuh,rhow,rhoh;
  doublereal gamma,alpha,ff,q,w_wp,w_fos,w_fc;
  doublereal mu_s_max,omegaf,chi,beta,a,dw,dh;
  doublereal pi, num_periods,L,powH;
  doublereal sigma,G,gamma_to,lamb,midlamb,mu_s,mu,tras;
  pi = atan(1.) * 4.0;  
  
  pi = atan(1.) * 4.0;  
  /* Defining the parameters */
  p=par[1+F2C];
  lamb_max=par[2+F2C];
  lamb_min=par[3+F2C];
  eta=par[4+F2C];
  nuw=par[5+F2C];
  nuh=par[6+F2C];
  rhow=par[7+F2C];
  rhoh=par[8+F2C];
  gamma=par[9+F2C];
  alpha=par[13+F2C];
  ff=par[14+F2C];
  q=par[15+F2C];
  w_wp=par[16+F2C];
  w_fos=par[17+F2C];
  w_fc=par[18+F2C];
  mu_s_max=par[19+F2C];
  omegaf=par[20+F2C];
  chi=par[21+F2C];
  beta=par[22+F2C];
  a=par[23+F2C];
  dw=par[24+F2C];
  dh=par[25+F2C];
  dummy_b=par[26+F2C];
  dummy_w=par[27+F2C];
  dummy_j=par[28+F2C];
  
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
  
  //L  = num_periods* (2.0*pi)/(1.0/resonance);
  /* Function Body */
  sigma= 100.0;
  G    = (1.0+eta*B)*(1.0+eta*B);
  midlamb = (lamb_max+lamb_min)/2.0;
  lamb = lamb_max + pow(chi,beta) * (lamb_min - lamb_max);
  gamma_to = gamma*(lamb/midlamb);
  mu_s = mu_s_max + pow((1.0-chi),beta) * (0.0 - mu_s_max);
  mu   = 1.0-mu_s*(1.0/(1.0 + exp(sigma*(W-(w_wp+w_fos)/2.0))));
  tras = gamma_to*G*W*B;
  
  f[0] = L * Bx;
  f[1] = L * Wx;
  f[2] = L * Hx;
  powH=2;
  f[3] = -L * (lamb*W*G*B*(1-B) - mu*B + dummy_b * Bx) ;//+ dummy_b * Bx;
  f[4] = -L * (1/dw)*((alpha*((B + q*ff)/(B + q)))*H - ((nuw)/(1 + rhow*B))*W - tras + (dummy_b + dummy_w)  * Wx) ;//+ (dummy_b + dummy_w)  * Wx;
  f[5] = -L * (1/(dh * powH * pow(H, (powH - 1))))*(p - (alpha*((B + q*ff)/(B + q))) * H - ((nuh)/(1 + rhoh*B)) * H + (dh * powH * (powH - 1) * pow(H, (powH - 2))* Hx*Hx)  + (dummy_b + dummy_j) * Hx);// + (dummy_b + dummy_j) * Hx;
 
  
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal x,
           doublereal *u, doublereal *par)
{
  /* local variables */
  doublereal B,W,H,Bx,Wx,Hx, X;
  doublereal dummy_b,dummy_w, dummy_j;
  doublereal p,lamb_max,lamb_min,eta,nuw,nuh,rhow,rhoh;
  doublereal gamma,alpha,ff,q,w_wp,w_fos,w_fc;
  doublereal mu_s_max,omegaf,chi,beta,a,dw,dh;
  doublereal pi, num_periods,L,powH;

  /*  defining the numerical value of the parameters */
  pi = atan(1.) * 4.0;
  /* Loading from file */
  FILE *myFile;
  myFile = fopen("tlm_parameters.txt", "r");
  //read file into array
  double numberArray[22];
  int i = 0;
  while (fscanf(myFile, "%lf", &numberArray[i]) != EOF) // && i!=20)
  {
	  //printf("%.15f ",numberArray[i]);
	  //printf("\n");
	  i++;
  }
  
  fclose(myFile);
  
  p=numberArray[0];
  lamb_max=numberArray[1];
  lamb_min=numberArray[2];
  eta=numberArray[3];
  nuw=numberArray[4];
  nuh=numberArray[5];
  rhow=numberArray[6];
  rhoh=numberArray[7];
  gamma=numberArray[8];
  alpha=numberArray[9];
  ff=numberArray[10];
  q=numberArray[11];
  w_wp=numberArray[12];
  w_fos=numberArray[13];
  w_fc=numberArray[14];
  mu_s_max=numberArray[15];
  omegaf=numberArray[16];
  chi=numberArray[17];
  beta=numberArray[18];
  a=numberArray[19];
  dw=numberArray[20];
  dh=numberArray[21];

  dummy_b = (doublereal).00000000000000002;
  dummy_w = (doublereal).00000000000000003;
  dummy_j = (doublereal).00000000000000005;
  
  printf("Starting continuation for chi= %.4f \n",chi);

  /* load into internal parameters */
  par[1+F2C]=p;
  par[2+F2C]=lamb_max;
  par[3+F2C]=lamb_min;
  par[4+F2C]=eta;
  par[5+F2C]=nuw;
  par[6+F2C]=nuh;
  par[7+F2C]=rhow;
  par[8+F2C]=rhoh;
  par[9+F2C]=gamma;
  par[11+F2C]=20.0;
  par[13+F2C]=alpha;
  par[14+F2C]=ff;
  par[15+F2C]=q;
  par[16+F2C]=w_wp;
  par[17+F2C]=w_fos;
  par[18+F2C]=w_fc;
  par[19+F2C]=mu_s_max;
  par[20+F2C]=omegaf;
  par[21+F2C]=chi;
  par[22+F2C]=beta;
  par[23+F2C]=a;
  par[24+F2C]=dw;
  par[25+F2C]=dh;
  par[26+F2C] = dummy_b;
  par[27+F2C] = dummy_w;
  par[28+F2C] = dummy_j;

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

  //fb[0] = u0[3];   // X_left     = 0
  
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

