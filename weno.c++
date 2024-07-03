#if 0
g++ -Wall -o weno -std=gnu++11 weno.cc
exit
#endif

#include <cmath>
#include <cstdio>
#include <iostream>
using namespace std;

#define Ni 10
#define Nj 10
#define Nk 10

#define dxyz 0.1
// avoid any extrema in our test monomials
#define x0 1.
#define y0 1.
#define z0 1.

#define CCTK_GFINDEX3D(cctkGH, i, j, k) ((i) + (j) * Ni + (k) * (Ni*Nj))
#define DECLARE_CCTK_PARAMETERS do {} while(0)
#define CCTK_REAL double
#define cGH void
#define restrict /**/

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MAX3(a,b,c) MAX(MAX((a), (b)), (c))

const int GRHydro_stencil = 3;
const double weno_eps = 1e-10;

template <typename T> static inline T SQR (T const & x) { return x*x; }
template <typename T> static inline T CUBE (T const & x) { return x*x*x; }

/**
   WENO5 reconstruction operator.
   Supports standard WENO5 (with and without adaptive epsilon), and WENO-z.
*/
template <bool do_wenoz, bool do_wham, bool do_adaptive_epsilon, int dir>
void apply(const int nx, const CCTK_REAL* const restrict a,
      CCTK_REAL* const restrict aminus, CCTK_REAL* const restrict aplus,
      const cGH* const cctkGH, const int j, const int k)
{
   DECLARE_CCTK_PARAMETERS;
   
#define A(i_) (a[ijk[i_]])
#define Aplus(i_) (aplus[ijk[i_]])
#define Aminus(i_) (aminus[ijk[i_]])
#define Acenter(i_) (aplus[ijk[i_]])
   
   for (int i=GRHydro_stencil-1; i < nx-GRHydro_stencil+1; ++i)
   {
      const int ijk[5] = {
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i-2, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i-2, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i-2), 
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i-1, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i-1, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i-1),
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i  , j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i  , k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i  ),
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i+1, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i+1, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i+1),
                            dir ==0 ? (int)CCTK_GFINDEX3D(cctkGH, i+2, j, k) : dir ==1 ? (int)CCTK_GFINDEX3D(cctkGH, j, i+2, k) : (int)CCTK_GFINDEX3D(cctkGH, j, k, i+2)
                         };
                       
   
         
   
      static_assert (! (do_wenoz && do_adaptive_epsilon), "Adaptive_epsilon not supported for WENO-Z");

      if (do_wenoz)
      {
         static const CCTK_REAL 
                         weno_coeffs[3][5] = { { 2.0/6.0, -7.0/6.0, 11.0/6.0, 0,        0 }, 
                                               { 0,       -1.0/6.0, 5.0/6.0,  2.0/6.0,  0 },
                                               { 0,        0,       2.0/6.0,  5.0/6.0, -1.0/6.0 } };
      
         const CCTK_REAL beta1 = 13.0/12.0*SQR(A(0)-2.0*A(1)+A(2)) + 1.0/4.0*SQR(A(0)-4.0*A(1)+3.0*A(2));
         const CCTK_REAL beta2 = 13.0/12.0*SQR(A(1)-2.0*A(2)+A(3)) + 1.0/4.0*SQR(A(1)-A(3));
         const CCTK_REAL beta3 = 13.0/12.0*SQR(A(2)-2.0*A(3)+A(4)) + 1.0/4.0*SQR(3.0*A(2)-4.0*A(3)+A(4));
            
            
         //    Compute weights according to weno-z alorithm
         const CCTK_REAL wbarplus1 = 1.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta1));
         const CCTK_REAL wbarplus2 = 3.0/5.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta2));
         const CCTK_REAL wbarplus3 = 3.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta3));

         const CCTK_REAL wbarminus1 = 3.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta1));
         const CCTK_REAL wbarminus2 = 3.0/5.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta2));
         const CCTK_REAL wbarminus3 = 1.0/10.0 * (1.0 + abs(beta1-beta3) / (weno_eps + beta3));
         
         const CCTK_REAL iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3);
         
         const CCTK_REAL wplus1 = wbarplus1 * iwbarplussum;
         const CCTK_REAL wplus2 = wbarplus2 * iwbarplussum;
         const CCTK_REAL wplus3 = wbarplus3 * iwbarplussum;
         
         const CCTK_REAL iwbarminussum = 1.0 / (wbarminus1 + wbarminus2 + wbarminus3);
         
         const CCTK_REAL wminus1 = wbarminus1 * iwbarminussum;
         const CCTK_REAL wminus2 = wbarminus2 * iwbarminussum;
         const CCTK_REAL wminus3 = wbarminus3 * iwbarminussum;
         
         //    Calculate the reconstruction
         Aplus(2) = 0;
         Aminus(2) = 0;
         for (int n=0; n < 5; ++n) {
               Aplus(2) += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);
               Aminus(2) += (wminus1 * weno_coeffs[2][4-n]
                           + wminus2 * weno_coeffs[1][4-n]
                           + wminus3 * weno_coeffs[0][4-n]) * A(n);
         }
      } else if(do_wham) {
         
         // equ. of https://arxiv.org/abs/0704.2608 (WHAM paper) beta same for
         // center to interface, average to center and center to average
         // reconstruction (appendix A3 before equ A18)
         static const CCTK_REAL beta_shu[3][6] = { { 4.0/3.0,  -19.0/3.0, 25.0/3.0, 11.0/3.0, -31.0/3.0, 10.0/3.0 },
                                      { 4.0/3.0,  -13.0/3.0, 13.0/3.0, 5.0/3.0,  -13.0/3.0, 4.0/3.0 },
                                      { 10.0/3.0, -31.0/3.0, 25.0/3.0, 11.0/3.0, -19.0/3.0, 4.0/3.0 } };

         // these are from equ. (18) and (19) when substituting into each other
         // and sorting by the index i-r+j         
         /**static const CCTK_REAL weno_coeffs[3][5] = { { 1.0/3.0, -7.0/6.0, 11.0/6.0, 0,        0 },
                                                    {   0,       -1.0/6.0,  5.0/6.0, 1.0/3.0 , 0 },
                                                    {   0,       0,         1.0/3.0, 5.0/6.0, -1.0/6.0 } };**/
          //Only de-averaging for our particular test case
          static const CCTK_REAL weno_coeffs[3][5] = { { -1.0/24.0, 1.0/12.0, 23.0/24.0, 0,        0 },
                                                    {   0,       -1.0/24.0,  13.0/12.0, -1.0/24.0 , 0 },
                                                    {   0,       0,         23.0/24.0, 1.0/12.0, -1.0/24.0 } };
                                      
         //    Compute smoothness indicators
         //    This is from Tchekhovskoy et al 2007 (WHAM code paper).
         CCTK_REAL beta1  = beta_shu[0][0]*SQR(A(0))
                  + beta_shu[0][1]*A(0)*A(1)
                  + beta_shu[0][2]*SQR(A(1))
                  + beta_shu[0][3]*A(0)*A(2)
                  + beta_shu[0][4]*A(1)*A(2)
                  + beta_shu[0][5]*SQR(A(2));
         
         CCTK_REAL beta2  = beta_shu[1][0]*SQR(A(1))
                  + beta_shu[1][1]*A(1)*A(2)
                  + beta_shu[1][2]*SQR(A(2))
                  + beta_shu[1][3]*A(1)*A(3)
                  + beta_shu[1][4]*A(2)*A(3)
                  + beta_shu[1][5]*SQR(A(3));
         
         CCTK_REAL beta3  = beta_shu[2][0]*SQR(A(2))
                  + beta_shu[2][1]*A(2)*A(3)
                  + beta_shu[2][2]*SQR(A(3))
                  + beta_shu[2][3]*A(2)*A(4)
                  + beta_shu[2][4]*A(3)*A(4)
                  + beta_shu[2][5]*SQR(A(4));
         
         
         if (do_adaptive_epsilon) {
            const CCTK_REAL vnorm = (SQR(A(0)) + SQR(A(1)) + SQR(A(2)) + SQR(A(3)) + SQR(A(4)));
               
            beta1 += 100.0*weno_eps*(vnorm + 1.0);
            beta2 += 100.0*weno_eps*(vnorm + 1.0);
            beta3 += 100.0*weno_eps*(vnorm + 1.0);
               
            const CCTK_REAL ibetanorm = 1.0 / (beta1 + beta2 + beta3);
               
            beta1 *= ibetanorm;
            beta2 *= ibetanorm;
            beta3 *= ibetanorm;
         }
         /**
         const CCTK_REAL wbarplus1 = 1.0/16.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarplus2 = 5.0/8.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarplus3 = 5.0/16.0 / SQR(weno_eps + beta3);**/
         
         //Updated weights for Average to Center
         const CCTK_REAL wbarplus1 = -9.0/80.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarplus2 = 49.0/40.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarplus3 = -9.0/80.0 / SQR(weno_eps + beta3);
         
         const CCTK_REAL iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3);
         
         const CCTK_REAL wplus1 = wbarplus1 * iwbarplussum;
         const CCTK_REAL wplus2 = wbarplus2 * iwbarplussum;
         const CCTK_REAL wplus3 = wbarplus3 * iwbarplussum;

         const CCTK_REAL wbarminus1 = 5.0/16.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarminus2 = 5.0/8.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarminus3 = 1.0/16.0 / SQR(weno_eps + beta3);
         
         const CCTK_REAL iwbarminussum = 1.0 / (wbarminus1 + wbarminus2 + wbarminus3);
         
         const CCTK_REAL wminus1 = wbarminus1 * iwbarminussum;
         const CCTK_REAL wminus2 = wbarminus2 * iwbarminussum;
         const CCTK_REAL wminus3 = wbarminus3 * iwbarminussum;
                                         
         //    Calculate the reconstruction
         Aplus(2) = 0;
         Aminus(2) = 0;
         Acenter(2) = 0;
         for (int n=0; n < 5; ++n) {
               /**Aplus(2) += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);**/
               /**Aminus(2) += (wminus1 * weno_coeffs[2][4-n]
                           + wminus2 * weno_coeffs[1][4-n]
                           + wminus3 * weno_coeffs[0][4-n]) * A(n);**/
                  Aminus(2) += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);
              //Interpolating Cell Center (2nd index not verified? Should it be 2-n?)
              //Aplus(2)       += (wplus1 * weno_coeffs[0][n]
                //          + wplus2 * weno_coeffs[1][n]
                  //        + wplus3 * weno_coeffs[2][n]) * A(n);
             /* Acenter(2)       += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);*/
                          
         }

         /**
         *
         *Reconstruction of values at the interface
         *
         */
         /**
         *WENO coeffs for center->left
         */
         static const CCTK_REAL weno_coeffs_recon[3][5] = { { 3.0/8.0, -5.0/4.0, 15.0/8.0, 0,        0 },
                                                    {   0,       -1.0/8.0,  3.0/4.0, 3.0/8.0 , 0 },
                                                    {   0,       0,         3.0/8.0, 3.0/4.0, -1.0/8.0 } };;
        //    Compute smoothness indicators
         //    This is from Tchekhovskoy et al 2007 (WHAM code paper).
         CCTK_REAL beta4  = beta_shu[0][0]*SQR(Aminus(0))
                  + beta_shu[0][1]*Aminus(0)*Aminus(1)
                  + beta_shu[0][2]*SQR(Aminus(1))
                  + beta_shu[0][3]*Aminus(0)*Aminus(2)
                  + beta_shu[0][4]*Aminus(1)*Aminus(2)
                  + beta_shu[0][5]*SQR(Aminus(2));

         CCTK_REAL beta5  = beta_shu[1][0]*SQR(Aminus(1))
                  + beta_shu[1][1]*Aminus(1)*Aminus(2)
                  + beta_shu[1][2]*SQR(Aminus(2))
                  + beta_shu[1][3]*Aminus(1)*Aminus(3)
                  + beta_shu[1][4]*Aminus(2)*Aminus(3)
                  + beta_shu[1][5]*SQR(Aminus(3));
         
         CCTK_REAL beta6  = beta_shu[2][0]*SQR(Aminus(2))
                  + beta_shu[2][1]*Aminus(2)*Aminus(3)
                  + beta_shu[2][2]*SQR(Aminus(3))
                  + beta_shu[2][3]*Aminus(2)*Aminus(4)
                  + beta_shu[2][4]*Aminus(3)*Aminus(4)
                  + beta_shu[2][5]*SQR(Aminus(4));

          if (do_adaptive_epsilon) {
            const CCTK_REAL vnorm2 = (SQR(Aminus(0)) + SQR(Aminus(1)) + SQR(Aminus(2)) + SQR(Aminus(3)) + SQR(Aminus(4)));
               
            beta4 += 100.0*weno_eps*(vnorm2 + 1.0);
            beta5 += 100.0*weno_eps*(vnorm2 + 1.0);
            beta6 += 100.0*weno_eps*(vnorm2 + 1.0);
               
            const CCTK_REAL ibetanorm2 = 1.0 / (beta4 + beta5 + beta6);
               
            beta4 *= ibetanorm2;
            beta5 *= ibetanorm2;
            beta6 *= ibetanorm2;
         }

         //Updated weights for Center to Left
         const CCTK_REAL wbarplus1recon = 1.0/16.0 / SQR(weno_eps + beta4);
         const CCTK_REAL wbarplus2recon = 5.0/8.0 / SQR(weno_eps + beta5);
         const CCTK_REAL wbarplus3recon = 5.0/16.0 / SQR(weno_eps + beta6);
         
         const CCTK_REAL iwbarplussumrecon = 1.0 / (wbarplus1recon + wbarplus2recon + wbarplus3recon);
         
         const CCTK_REAL wplus1recon = wbarplus1recon * iwbarplussumrecon;
         const CCTK_REAL wplus2recon = wbarplus2recon * iwbarplussumrecon;
         const CCTK_REAL wplus3recon = wbarplus3recon * iwbarplussumrecon;
                                         
         //    Calculate the reconstruction
         Aplus(2) = 0;
         for (int n=0; n < 5; ++n) {
              //Interpolating Cell Center to Left Interface
              Aplus(2)       += (wplus1 * weno_coeffs_recon[0][n]
                          + wplus2 * weno_coeffs_recon[1][n]
                          + wplus3 * weno_coeffs_recon[2][n]) * A(n);
                          
         }

          
        

         
      } else {
         
         static const CCTK_REAL beta_shu[3][6] = { { 4.0/3.0,  -19.0/3.0, 25.0/3.0, 11.0/3.0, -31.0/3.0, 10.0/3.0 },
                                      { 4.0/3.0,  -13.0/3.0, 13.0/3.0, 5.0/3.0,  -13.0/3.0, 4.0/3.0 },
                                      { 10.0/3.0, -31.0/3.0, 25.0/3.0, 11.0/3.0, -19.0/3.0, 4.0/3.0 } };
         static const CCTK_REAL weno_coeffs[3][5] = { { 3.0/8.0, -5.0/4.0, 15.0/8.0, 0,      0 },
                                         { 0,       -1.0/8.0, 3.0/4.0,  3.0/8.0, 0 },
                                         { 0,       0,        3.0/8.0,  3.0/4.0, -1.0/8.0 } };
                                      
         //    Compute smoothness indicators
         //    This is from Tchekhovskoy et al 2007 (WHAM code paper).
         CCTK_REAL beta1  = beta_shu[0][0]*SQR(A(0))
                  + beta_shu[0][1]*A(0)*A(1)
                  + beta_shu[0][2]*SQR(A(1))
                  + beta_shu[0][3]*A(0)*A(2)
                  + beta_shu[0][4]*A(1)*A(2)
                  + beta_shu[0][5]*SQR(A(2));
         
         CCTK_REAL beta2  = beta_shu[1][0]*SQR(A(1))
                  + beta_shu[1][1]*A(1)*A(2)
                  + beta_shu[1][2]*SQR(A(2))
                  + beta_shu[1][3]*A(1)*A(3)
                  + beta_shu[1][4]*A(2)*A(3)
                  + beta_shu[1][5]*SQR(A(3));
         
         CCTK_REAL beta3  = beta_shu[2][0]*SQR(A(2))
                  + beta_shu[2][1]*A(2)*A(3)
                  + beta_shu[2][2]*SQR(A(3))
                  + beta_shu[2][3]*A(2)*A(4)
                  + beta_shu[2][4]*A(3)*A(4)
                  + beta_shu[2][5]*SQR(A(4));
         
         
         if (do_adaptive_epsilon) {
            const CCTK_REAL vnorm = (SQR(A(0)) + SQR(A(1)) + SQR(A(2)) + SQR(A(3)) + SQR(A(4)));
               
            beta1 += 100.0*weno_eps*(vnorm + 1.0);
            beta2 += 100.0*weno_eps*(vnorm + 1.0);
            beta3 += 100.0*weno_eps*(vnorm + 1.0);
               
            const CCTK_REAL ibetanorm = 1.0 / (beta1 + beta2 + beta3);
               
            beta1 *= ibetanorm;
            beta2 *= ibetanorm;
            beta3 *= ibetanorm;
         }
         
         const CCTK_REAL wbarplus1 = 1.0/16.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarplus2 = 5.0/8.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarplus3 = 5.0/16.0 / SQR(weno_eps + beta3);
         
         const CCTK_REAL iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3);
         
         const CCTK_REAL wplus1 = wbarplus1 * iwbarplussum;
         const CCTK_REAL wplus2 = wbarplus2 * iwbarplussum;
         const CCTK_REAL wplus3 = wbarplus3 * iwbarplussum;

         const CCTK_REAL wbarminus1 = 5.0/16.0 / SQR(weno_eps + beta1);
         const CCTK_REAL wbarminus2 = 5.0/8.0 / SQR(weno_eps + beta2);
         const CCTK_REAL wbarminus3 = 1.0/16.0 / SQR(weno_eps + beta3);
         
         const CCTK_REAL iwbarminussum = 1.0 / (wbarminus1 + wbarminus2 + wbarminus3);
         
         const CCTK_REAL wminus1 = wbarminus1 * iwbarminussum;
         const CCTK_REAL wminus2 = wbarminus2 * iwbarminussum;
         const CCTK_REAL wminus3 = wbarminus3 * iwbarminussum;
                                         
         //    Calculate the reconstruction
         Aplus(2) = 0;
         Aminus(2) = 0;
         for (int n=0; n < 5; ++n) {
               Aplus(2) += (wplus1 * weno_coeffs[0][n]
                          + wplus2 * weno_coeffs[1][n]
                          + wplus3 * weno_coeffs[2][n]) * A(n);
               Aminus(2) += (wminus1 * weno_coeffs[2][4-n]
                           + wminus2 * weno_coeffs[1][4-n]
                           + wminus3 * weno_coeffs[0][4-n]) * A(n);
              
         }
      }
   }
}

CCTK_REAL data[Nk][Nj][Ni];
CCTK_REAL data_plus_weno[MAX3(Nk,Nj,Ni)];
CCTK_REAL data_minus_weno[MAX3(Nk,Nj,Ni)];
CCTK_REAL data_plus_wenoz[MAX3(Nk,Nj,Ni)];
CCTK_REAL data_minus_wenoz[MAX3(Nk,Nj,Ni)];
CCTK_REAL data_plus_wham[MAX3(Nk,Nj,Ni)];
CCTK_REAL data_minus_wham[MAX3(Nk,Nj,Ni)];
//Data for cell center interpolation
CCTK_REAL data_center_wham[MAX3(Nk,Nj,Ni)];

#define ORDER 2

#if (ORDER == 1)
CCTK_REAL fun(double x, double y, double z) {
  return 2.*x;
}
CCTK_REAL Fun(double x, double y, double z) {
  return (SQR(x+0.5*dxyz) - SQR(x-0.5*dxyz)) / dxyz;
}
#elif (ORDER == 2)
CCTK_REAL fun(double x, double y, double z) {
  return 3.*SQR(x);
}
CCTK_REAL Fun(double x, double y, double z) {
  return (CUBE(x+0.5*dxyz) - CUBE(x-0.5*dxyz)) / dxyz;
}
#endif

void test(CCTK_REAL (*f)(double x, double y, double z), const char *label) {
  printf("Testing %d order %s\n", ORDER, label);

  // make up some dummy data
  for(int k = 0 ; k < Nk ; k++) {
    for(int j = 0 ; j < Nj ; j++) {
      for(int i = 0 ; i < Ni ; i++) {
        // coords of centers of cells
        double x = x0 + i * dxyz;
        double y = y0 + j * dxyz;
        double z = z0 + k * dxyz;
        // should be Fun for reconstruction and fun for interpolation
        data[k][j][i] = f(x, y, z);
      }
    }
  }

  // reconstruct along x
  apply<false, false, false, 0>(Ni, &data[0][0][0], data_minus_weno, data_plus_weno, NULL, 0, 0); // WENO as in GRHydro
  apply<true,  false, false, 0>(Ni, &data[0][0][0], data_minus_wenoz, data_plus_wenoz, NULL, 0, 0); // WENOZ as in GRHydro
  //apply<false, true,  false, 0>(Ni, &data[0][0][0], data_minus_wham, data_plus_wham, NULL, 0, 0); // WENNO as in WHAM paper
  //Modified apply function for testing de-averaging
  apply<false, true,  false, 0>(Ni, &data[0][0][0], data_center_wham, data_plus_wham, NULL, 0, 0); // WENNO as in WHAM paper

  // some output
  /**
  printf("GRhydro-WENO:\n");
  for(int i = 0 ; i < Ni ; i++) {
    double x = x0 + i * dxyz;
    printf("%d %g %g : %g =?= %g  %g =?= %g\n", i, x, data[0][0][i], data_minus_weno[i], fun(x-0.5*dxyz, 0., 0.), data_plus_weno[i], fun(x+0.5*dxyz, 0., 0.));
  }
  printf("GRhydro-WENOZ:\n");
  for(int i = 0 ; i < Ni ; i++) {
    double x = x0 + i * dxyz;
    printf("%d %g %g : %g =?= %g  %g =?= %g\n", i, x, data[0][0][i], data_minus_wenoz[i], fun(x-0.5*dxyz, 0., 0.), data_plus_wenoz[i], fun(x+0.5*dxyz, 0., 0.));
  }
  printf("WHAM-WENO:\n");
  for(int i = 0 ; i < Ni ; i++) {
    double x = x0 + i * dxyz;
    printf("%d %g %g : %g =?= %g  %g =?= %g\n", i, x, data[0][0][i], data_minus_wham[i], fun(x-0.5*dxyz, 0., 0.), data_plus_wham[i], fun(x+0.5*dxyz, 0., 0.));
  }**/

  printf("WHAM-WENO:\n");
  for(int i = 0 ; i < Ni ; i++) {
    double x = x0 + i * dxyz;
    printf("%d %g %g : %g =?= %g  %g =?= %g\n", i, x, data[0][0][i], data_center_wham[i], fun(x, 0., 0.), data_plus_wham[i], fun(x+0.5*dxyz, 0., 0.));
  }

  puts("");
}

int main(void) {

  test(fun, "interpolation");
  test(Fun, "reconstruction");
  

  return 0;
}
