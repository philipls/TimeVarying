/*  Building:  
mex wiensustn400x.c -lgsl -lgslcblas  -lm

*/

#include "mex.h"
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#define ssz 400
#define bsz (2 * ssz + 1)
#define NP 7  /* Size of parameter vector */
#define A1 3   /* Locations of criteria in parameter vector */
#define A2  4
#define SIGMA1 5
#define SIGMA2 6
/* #define trace */

typedef double big[bsz];
typedef double small[ssz];

double max(double a, double b) {
    double mx;
    if (a >= b) {
        mx = a;
    } else {
        mx = b;
    }
    return mx;
}              

void wiensustn400a(double *T, double *G1, double *G2, double *Ig1, double *Ig2, double *Stat,
         double *P, double h) {

 /* 
    ===============================================================================
    wiensustn400 (400 step version).
    ----------
   
    Time-changed diffusion, adds time offset. 
    P = [mu, beta, a1, a2, sigma1, sigma2, t0].
    beta is a range parameter for a 3-stage gamma drift rate function.
    sigma1 and sigma2 are variable and fixed components of diffusion. 
    ================================================================================
 */
    double eps = 0.0005;
    double pi = 3.141592654;
    double m, mu, alpha, beta, a1, a2, sigma1, sigma2, sigma1_sq, sigma2_sq, t0, dPhi, w2;
    int i, j, k;
    char ch;

    small Th, IntTimeBase, IntMu, IntSigma2, Mu, Sigma2, Phi, DPhi, DPhiOnDPsi,
          Psib_1, Psib_2, DPsib_1, DPsib_2,
          D01_1, D01_2, D02_1, D02_2, 
          Psi0_1, Psi0_2,  
          D1_11, D1_12, D1_21, D1_22,
          D2_11, D2_12, D2_21, D2_22,
          V1_11, V1_12, V1_21, V1_22,
          Psi_11, Psi_12, Psi_21, Psi_22;
 

    big TimeBase, Tbase;

   /*
    -------------------------------------------------------------------------------
    Check input parameters and other housekeeping.
    -------------------------------------------------------------------------------
   */
    
   mu = P[0];
   beta = P[1];
   m = P[2];
   a1 = P[3];
   a2 = P[4];
   sigma1 = P[5];
   sigma2 = P[6];

   sigma1_sq = sigma1 * sigma1;
   sigma2_sq = sigma2 * sigma2;
   
#ifdef trace
   mexPrintf("mu =    %6.2f \n", mu);
   mexPrintf("beta =  %6.2f \n", beta);
   mexPrintf("m =    %6.2f \n", m);
   mexPrintf("a1 =    %6.2f \n", a1);
   mexPrintf("a2 =    %6.2f \n", a2);
   mexPrintf("sigma1 = %6.2f \n", sigma1);
   mexPrintf("sigma2 = %6.2f \n", sigma2);
   mexPrintf("h    =  %6.3f \n", h);
#endif

   
   /* 
      ---------------------------------------------------------------------
      Compute various time vectors. T is standard time vector that is returned
      by the function.  T1 is h/2 scale vector used in kernel calculations.
      Th is h/2	shifted vector used in calculating RTs.
      --------------------------------------------------------------------- 
   */   
    for (i = 0; i < bsz; i++) {
        Tbase[i] = max((i * h / 2 - t0), 0);
        /* Tbase[i] = i * h / 2;
        mexPrintf("variables: %6.3f \n", Tbase[i]); */
    }    

    for (i = 0; i < ssz; i++) {
		T[i] = (i+1) * h;
		Th[i] = T[i] - h / 2;
    }
    
    /* -------------------------------------------------------------------
      Compute the common time base (m-stage gamma), then integrate it.
       ----------------------------------------------------------------- */
    
    for (j = 0; j < bsz; j++) {
        TimeBase[j] = gsl_sf_gamma_inc_P(m, beta * Tbase[j]);
         /* mexPrintf("j = %6d, TimeBase[j] =  %6.3f \n", j, TimeBase[j]); */  
    }
    for (i = 0; i < ssz; i++) {
       j = 2 * (i + 1);
       if (i > 0) {
	    IntTimeBase[i] = IntTimeBase[i-1]
			+ (TimeBase[j-2] + 4 * TimeBase[j-1] + TimeBase[j]) * h / 6;
            Mu[i] = mu * TimeBase[j];
            Sigma2[i] = sigma1_sq * TimeBase[j] + sigma2_sq;
            IntMu[i] = mu * IntTimeBase[i];
            IntSigma2[i] = sigma1_sq * IntTimeBase[i] + sigma2_sq * T[i]; 
       } else {
	        IntTimeBase[i] = (TimeBase[j-2] + 4 * TimeBase[j-1] + TimeBase[j]) * h / 6;
            Mu[i] = 0;
            Sigma2[i] = 0;
            IntMu[i] = mu * IntTimeBase[i];
            IntSigma2[i] = sigma1_sq * IntTimeBase[i] + sigma2_sq * T[i];
       }
        /* mexPrintf("variables: %6.3f %6.3f %6.3f %6.3f, %6.3f  \n",
                      T[i], Mu[i], Sigma2[i], IntMu[i], IntSigma2[i]); */
    }
    
     /*
       -------------------------------------------------------------------
       Compute Phi(t), Psib1(a, t), Psib2(a, t) and their derivatives.
       -------------------------------------------------------------------
     */  
     for (i = 0; i < ssz; i++) {
           DPhi[i] = Sigma2[i];
           Phi[i] = IntSigma2[i];
           Psib_1[i] = a1 - IntMu[i];
           Psib_2[i] = a2 - IntMu[i];
           DPsib_1[i] = -Mu[i];
	       DPsib_2[i] = -Mu[i];
           DPhiOnDPsi[i] = Sigma2[i];
           /* mexPrintf("variables: %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f \n",
                      DPhi[i], Phi[i], Psib_1[i], Psib_2[i], DPsib_1[i], DPsib_2[i], 
                      DPhiOnDPsi[i]); */
     }
 
     /*
       --------------------------------------------------------------------
       Compute constant kernels 2*Psi[a1(kh),kh|0,0], 2*Psi[a2(kh),kh|0,0]
       --------------------------------------------------------------------
     */
     /* ignore the Phi_prime_x multiple on density because 1.0 */
     for (i = 0; i < ssz; i++) {
         D01_1[i] = exp(-0.5 * Psib_1[i] * Psib_1[i] / Phi[i]) /
                    sqrt(2.0 * pi * Phi[i]);

         D01_2[i] = exp(-0.5 * Psib_2[i] * Psib_2[i] / Phi[i]) /
                    sqrt(2.0 * pi * Phi[i]);
         
         D02_1[i] = DPsib_1[i] - Psib_1[i] * DPhiOnDPsi[i] / Phi[i];      
         D02_2[i] = DPsib_2[i] - Psib_2[i] * DPhiOnDPsi[i] / Phi[i];
         
         Psi0_1[i] = D01_1[i] * D02_1[i];
         Psi0_2[i] = D01_2[i] * D02_2[i];
         /*  mexPrintf("variables: %6.3f %6.3f %6.3f, %6.3f, %6.3f, %6.3f  \n",
                      T[i], Phi[i], D01_1[i], D02_1[i], D01_2[i], D02_2[i]); */
     }             
     
     G1[0] = -Psi0_1[0];
     G2[0] =  Psi0_2[0];
     Ig1[0] = G1[0] * h;
     Ig2[0] = G2[0] * h;
     /* mexPrintf("variables: %6.3f, %6.3f, %6.3f, %6.3f  \n",
                      G1[0], G2[0], Ig1[0], Ig2[0]); */
     k = 0;
     while (k < ssz && (Ig1[k] + Ig2[k]) < 1.0 - eps) {                       
          /* mexPrintf("k: %6d \n", k); */
        /*   
          -----------------------------------------------------------------------
          Compute variable kernels 2*Psi[a(kh),kh|a(jh),jh] for fixed k;  j=0:k-1
          -----------------------------------------------------------------------
        */
        
        G1[k] = -Psi0_1[k];
        G2[k] =  Psi0_2[k];
        
        for (j = 0; j < k; j++) {
             V1_11[j] = Psib_1[k] - Psib_1[j];
             V1_12[j] = Psib_1[k] - Psib_2[j];                
             V1_21[j] = Psib_2[k] - Psib_1[j];     
             V1_22[j] = Psib_2[k] - Psib_2[j];
             
             if (Phi[k] - Phi[j] < 0) {
                 mexErrMsgTxt("wiensustn400: Neg sqrt");
             }              
	     dPhi = Phi[k] - Phi[j];           
             w2 = sqrt(2 * pi * dPhi); 
             
             D1_11[j] = exp( -0.5 * V1_11[j] * V1_11[j] / dPhi) / w2;
             D1_12[j] = exp( -0.5 * V1_12[j] * V1_12[j] / dPhi) / w2;
             D1_21[j] = exp( -0.5 * V1_21[j] * V1_21[j] / dPhi) / w2;
             D1_22[j] = exp( -0.5 * V1_22[j] * V1_22[j] / dPhi) / w2;
             
             D2_11[j] = DPsib_1[k] - V1_11[j] * DPhiOnDPsi[k] / dPhi;
             D2_12[j] = DPsib_1[k] - V1_12[j] * DPhiOnDPsi[k] / dPhi;
             D2_21[j] = DPsib_2[k] - V1_21[j] * DPhiOnDPsi[k] / dPhi;           
             D2_22[j] = DPsib_2[k] - V1_22[j] * DPhiOnDPsi[k] / dPhi; 
             
             Psi_11[j] = D1_11[j] * D2_11[j];
             Psi_12[j] = D1_12[j] * D2_12[j];            
             Psi_21[j] = D1_21[j] * D2_21[j];            
             Psi_22[j] = D1_22[j] * D2_22[j];

             G1[k] = G1[k] + h * (Psi_11[j] * G1[j] + Psi_12[j] * G2[j]); 
             G2[k] = G2[k] - h * (Psi_21[j] * G1[j] + Psi_22[j] * G2[j]);
             if (k > 0) {
                  Ig1[k] = Ig1[k-1] + (G1[k-1] + G1[k]) * h / 2;             
                  Ig2[k] = Ig2[k-1] + (G2[k-1] + G2[k]) * h / 2;
             } else {
                  Ig1[k] = G1[k] * h / 2;
                  Ig2[k] = G2[k] * h / 2;
             }   
         }
         k++;
     }
        
  /* 
     -----------------------------------------------------------------------------
      If all mass accounted for before ssz, fill arrays with appropriate values.
     -----------------------------------------------------------------------------
  */
     if (k < ssz - 1) {
          for (j = k; j < ssz; j++) {
               G1[j] = 0;
               G2[j] = 0;
               Ig1[j] = Ig1[k];
               Ig2[j] = Ig2[k];
          }
     }              
      
  /*
     -----------------------------------------------------------------------------
     Calculate means and response probabilities.
     First step of Simpson's rule integration is done at initialisation because 
     T = 0 is not represented in any of the arrays.  
     -----------------------------------------------------------------------------
  */
  
  Stat[0] = 4 * G1[0] + G1[1]; 
  Stat[1] = 4 * G2[0] + G2[1];
  Stat[2] = 4 * G1[0] * Th[0] + G1[1] * Th[1]; 
  Stat[3] = 4 * G2[0] * Th[0] + G2[1] * Th[1];
  for (j = 2; j <= ssz-1; j += 2) {
      Stat[0] = Stat[0] + G1[j-1] + 4 * G1[j] + G1[j+1];
      Stat[1] = Stat[1] + G2[j-1] + 4 * G2[j] + G2[j+1];
      Stat[2] = Stat[2] + G1[j-1] * Th[j-1]  + 4 * G1[j] * Th[j] + G1[j+1] * Th[j+1];
      Stat[3] = Stat[3] + G2[j-1] * Th[j-1]  + 4 * G2[j] * Th[j] + G2[j+1] * Th[j+1];
  }
  /* Simpson's rule constants of integration applied here */
  Stat[0] = Stat[0] * h / 3;
  Stat[1] = Stat[1] * h / 3;
 
  if (Stat[0] > 0) {
      Stat[2] = Stat[2] * h / (3 * Stat[0]);
  } else {
      Stat[2] = 0;
  }
  if (Stat[1] > 0) {
      Stat[3] = Stat[3] * h / (3 * Stat[1]);
  } else {
      Stat[3] = 0;
  }
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 /*
     -----------------------------------------------------------------------
     Matlab gateway routine to call wiensustn400.
     -----------------------------------------------------------------------
 */
 
double *T, *G1, *G2, *Ig1, *Ig2, *Stat, *P;
double h;
 
unsigned n, m;
     
    if (nrhs != 2) {
         mexErrMsgTxt("wiensustn400: Requires 2 input args.");
    } else if (nlhs != 6) {
        mexErrMsgTxt("wiensustn400: Requires 6 ouput args."); 
    }

    /*
      -----------------------------------------------------------------------
      Check all input argument dimensions.
      -----------------------------------------------------------------------   
    */
	
    /* P */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);       
    if (!mxIsDouble(prhs[0]) || !(m * n == NP)) {
        mexErrMsgTxt("wiensustn400: Wrong size P");
    } else {
        P = mxGetPr(prhs[0]);
    }
   if (P[A1] <= 0 || P[A2] >= 0) {
        mexPrintf("a1, a2: %6.2f, %6.2f \n", P[A1], P[A2]); 
        mexErrMsgTxt("wiensustn400: Bad criterion values");
    }
    if (P[SIGMA1] <= 0 || P[SIGMA2] <= 0) {
        mexPrintf("sigma1, sigma2: %6.2f, %6.2f \n", P[SIGMA1], P[SIGMA2]); 
        mexErrMsgTxt("wiensustn400: Bad standard deviations");
    } 
    
    /* h */
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);       
    if (!mxIsDouble(prhs[1]) || !(m * n == 1)) {
        mexErrMsgTxt("wiensustn400: h must be scalar");
    } else {
        h = mxGetScalar(prhs[1]);
    } 
    
    /*
      -----------------------------------------------------------------------
      Create output arrays (All 1 x ssz).
      -----------------------------------------------------------------------   
    */
    
    /* T */
    plhs[0] = mxCreateDoubleMatrix(1, ssz, mxREAL);
    T = mxGetPr(plhs[0]);
      
    /* G1 */
    plhs[1] = mxCreateDoubleMatrix(1, ssz, mxREAL);
    G1 = mxGetPr(plhs[1]);         
    
    /* G2 */
    plhs[2] = mxCreateDoubleMatrix(1, ssz, mxREAL);
    G2 = mxGetPr(plhs[2]);         
    
    /* Ig1 */
    plhs[3] = mxCreateDoubleMatrix(1, ssz, mxREAL);
    Ig1 = mxGetPr(plhs[3]);    
    
    /* Ig2 */
    plhs[4] = mxCreateDoubleMatrix(1, ssz, mxREAL);
    Ig2 = mxGetPr(plhs[4]);    
    
    /* Stat */
    plhs[5] = mxCreateDoubleMatrix(1, 4, mxREAL);
    Stat = mxGetPr(plhs[5]);    
    /*
      -----------------------------------------------------------------------
      Run the C-function wiensustn400.
      -----------------------------------------------------------------------   
    */    
    wiensustn400a(T, G1, G2, Ig1, Ig2, Stat, P, h);          
}           
           
    
