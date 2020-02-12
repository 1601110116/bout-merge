/**************************************************************************
 * Various differential operators defined on BOUT grid
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <globals.hxx>
#include <bout.hxx>
#include <boutexception.hxx>
#include <integrops.hxx>
#include <difops.hxx>
#include <vecops.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <invert_parderiv.hxx>
#include <fft.hxx>
#include <msg_stack.hxx>
#include <stdexcept>

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

#include <interpolation.hxx>

#include <math.h>
#include <stdlib.h>

/*******************************************************************************
 * iSign_kpar()
 * Integral operator calculates non-local parallel heat flux for given
 * parallel temperature profile
 *******************************************************************************/

const Field3D iSign_kpar(const Field3D &vin, BoutReal k0, const int bndry_c, const int N,
    const Field3D &SBC_value, const BoutReal alpha, const BoutReal beta)
{

  Field3D qout = 0.0;

  const BoutReal alp=alpha;     // default: 5.0
  const BoutReal bet=beta;      // default: 1.04

  // InvertPar *invpar;
  // invpar = InvertPar::Create();
  // much better if static:

  static InvertPar* invpar = InvertPar::Create();

  //-range of Lorentzian terms
  const int nMin=0;
  const int nMax=N;

  //-spectral shift parameter
  //replace the default by the right value here
  //-use empirical rule k0=min(kpar)/20., min(kpar)=2PI/Lpar, where Lpar is parallel domain size
  // Note: Lpar is in same units as your equations, if normalized then Lpar normalized same way

  //const BoutReal k0=1.0; //-this is the default
  //const BoutReal k0=(2*PI/Lpar)/20.;

  if ( nMax < 7)
    throw BoutException(
        "ERROR: At least 7 Lorentzians required for collisionless Landau closure fitting!\n");

  #ifdef DEBUG
    output.write("iSign_kpar_wcoll: bndry_c = %d, N = %d\n", bndry_c, N);
  #endif

  for (int n=nMin; n<nMax; n++)
  //for (int n=1; n<=1; n++) /*just for testing*/
    {

      BoutReal alpn=k0*pow(alp,n);
      BoutReal alpn2=alpn*alpn;

      //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(alpn2); //-free term
      invpar->setCoefB(-1.0);  //-d2/dz2 term

      /*
      switch(bndry_c) {
        case 0: {
          qout += invpar->solve(bet*alpn*Grad_par(vin), 0);
          break;
        }
        case 1: {
          qout += invpar->solve(bet*alpn*Grad_par(vin), 1);
          break;
        }
        case 2:
          break;
        case 3: {
          if (n == 0)
            qout += invpar->solve(bet*alpn*Grad_par(vin), 0);
          else
            qout += invpar->solve(bet*alpn*Grad_par(vin), 1);
          break;
        }
        defaut: {
          throw BoutException("iSign_kpar: unsupported boundary condition flag!\n");
          break;
        }
      }
      */
      if (bndry_c == 0) {
        // Dirichlet
        //output.write("iSign_kpar: [0]\n");
        qout += invpar->solve(bet*alpn*Grad_par(vin), 0);
      } else if (bndry_c == 1) {
        // Default: Zero Parallel Gradient
        //output.write("iSign_kpar: [1]\n");
        qout += invpar->solve(bet*alpn*Grad_par(vin), 1);
      } else if (bndry_c == 2) {
        // Sheath Boundary Condition
        // q^0 = SBC_value; q^n = 0
        //output.write("iSign_kpar: [2], SBC\n");
        if (n == 0)
          qout += invpar->solve(bet*alpn*Grad_par(vin), 2, SBC_value);
        else
          qout += invpar->solve(bet*alpn*Grad_par(vin), 0);
      } else if (bndry_c == 3) {
        // q^0: Neumann; q^n = 0
        //output.write("iSign_kpar: [3]\n");
        if (n == 0)
          qout += invpar->solve(bet*alpn*Grad_par(vin), 0);
        else
          qout += invpar->solve(bet*alpn*Grad_par(vin), 1);
      } else if (bndry_c == 4) {
        // q^n: Neumann; q^0 = 0
        //output.write("iSign_kpar: [4]\n");
        if (n == 0)
          qout += invpar->solve(bet*alpn*Grad_par(vin), 1);
        else
          qout += invpar->solve(bet*alpn*Grad_par(vin), 0);
      } else {
        // unsupported flag
        throw BoutException("iSign_kpar: unsupported boundary condition flag!\n");
      }
    }

#if 0
  delete invpar; //-deallocate the memory associated with the pointer
  //invpar->~InvertPar(); //-this should be equivalent to delete
#endif

  return(qout);
}

/*******************************************************************************/

/*******************************************************************************
 * iSign_kpar_wcoll()
 * Integral operator calculates non-local parallel heat flux with collisions for given
 * parallel temperature profile
 *******************************************************************************/

const Field3D iSign_kpar_wcoll(const Field3D &vin, const Field3D &k0, const int bndry_c,
    const int N, const Field3D &SBC_value)
{

   Field3D qout = 0.0;
   Field3D ddyvin;

   // ddyvin = DDY(vin);
   ddyvin = Grad_par(vin);
   ddyvin.applyBoundary();
   mesh->communicate(ddyvin);

   // InvertPar *invpar;
   // invpar = InvertPar::Create();
   // much better if static:

   static InvertPar* invpar = InvertPar::Create();

   //-range of Lorentzian terms
   const int nMin=0;
   const int nMax= (N == 0 ? 3 : N);

   // set coefficient for Lorentzian fitting
   const BoutReal *alpha;
   const BoutReal *beta;

   // N = 3:
   // Old coef. from  M.V. Umansky, et. al.,  Journal of Nuclear Materials 463, 506 (2015)
   // WARNING: deprecated
   const BoutReal alpha0[] = {0.0018, 0.0769, 2.4498};
   const BoutReal beta0[] = {0.1192, 0.4913, 2.1495};
   // new coef
   // N = 3:
   const BoutReal alpha3[] = {0.01315, 0.924, 14.1365};
   const BoutReal beta3[] = {0.2044, 1.3587, 8.9643};
   // N = 7:
   const BoutReal alpha7[] = { \
            0.007438, 0.6161, 5.9804, 37.9822, 234.3654, \
            1466.4331, 14981.4634};
   const BoutReal beta7[] = {0.1678, 1.1106, 5.6457, 33.1536, \
            202.738, 1254.2144, 9275.3323};
   // N = 12:
   const BoutReal alpha12[] = {0.001424, 0.20736, 2.5653, 14.927, \
             79.3050, 419.2399, 2215.7233, 11709.7857, \
             61885.2763, 327392.6096, 1773350.1566, 16903628.3745};
   const BoutReal beta12[] = {0.09419, 0.6741, 2.9628, 14.43958, \
             75.1106, 395.8293, 2090.8877, 11049.1471, \
             58392.0969, 308695.7371, 1645460.1472, 10794779.4293};
   const BoutReal alpha20[] = {
        0.000550164043066914, 0.10275855014401332, 1.5574968086648138, \
        8.981597596447177, 44.21499619448371, 215.21860764463176, \
        1046.9746330430294, 5093.378441306047, 24777.417109831247, \
        120516.8713406688, 586165.152129492, 2851132.6549336165, \
        13867844.765921066, 67444549.79120941, 327987870.7244059, \
        1595138256.7227132, 7758815048.312665, 37807952481.99112, \
        190502487382.14917, 1778418212793.4631};
   const BoutReal beta20[] = {
        0.0680957687739,  0.506661589811,  2.11683114157,  \
        9.29868746793,  44.186253098,  213.959779668,  \
        1039.85190758,  5057.62188658,  24603.1598059,  \
        119676.949309,  582095.153209,  2831266.12002,  \
        13771369.0447,  66979787.6208,  325739140.621,  \
        1584151555.15,  7704746851.2,  37497152109.4,  \
        1.84685042437e+11,  1.14854367316e+12};
   switch (N) {
     case 0: {
       // WARNING: deprecated
       alpha = alpha0;
       beta = beta0;
       break;
     }
     case 3: {
       alpha = alpha3;
       beta = beta3;
       break;
     }
     case 7: {
       alpha = alpha7;
       beta = beta7;
       break;
     }
     case 12: {
       alpha = alpha12;
       beta = beta12;
       break;
     }
     case 20: {
       alpha = alpha20;
       beta = beta20;
       break;
     }
     default: {
       throw BoutException(
         "ERROR: Invalid choice of number of Lorentzian terms. Must be in [3, 7, 12]\n");
       break;
     }
   }

   #ifdef DEBUG
     output.write("iSign_kpar_wcoll[3d]: bndry_c = %d, N = %d\n", bndry_c, N);
   #endif

   for (int n=nMin; n<nMax; n++)
      //for (int n=1; n<=1; n++) /*just for testing*/
   {

      //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(beta[n] * beta[n]); //-free term
      invpar->setCoefB(-1.0 / (k0 * k0));  //-d2/dz2 term

      //output.write("iSign_kpar_wcoll[3d]: bndry_c = %d, N = %d, n= %d\n", bndry_c, N, n);
      /* try {
      switch(bndry_c) {
        case 0: {
          output.write("iSign_kpar_wcoll: [3d, 0]\n");
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
          break;
        }
        case 1: {
          output.write("iSign_kpar_wcoll: [3d, 1]\n");
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
          break;
        }
        case 2: {
          break;
        }
        case 3: {
          output.write("iSign_kpar_wcoll: [3d, 3]\n");
          if (n == 0)
            qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
          else
            qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
          break;
        }
        defaut: {
          output.write("Default used!\n");
          // not output "Default used!\n" due to the throw statement
          // throw BoutException("iSign_kpar_wcoll: unsupported boundary condition flag!\n");
          throw invalid_argument( "received unsupport option\n" );
          return -1;
          break;
        }
      }
     }catch(invalid_argument& e) {
        output << "Error encountered during iSign_kpar_wcoll [3d, df, err]\n";
        output << e.what() << endl;
        return 1;
     } */
      if (bndry_c == 0) {
        // Dirichlet
        //output.write("iSign_kpar_wcoll: [3d, 0]\n");
        qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
      } else if (bndry_c == 1) {
        // Default: Zero Parallel Gradient
        //output.write("iSign_kpar_wcoll: [3d, 1]\n");
        qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
      } else if (bndry_c == 2) {
        // Sheath Boundary Condition
        // q^0 = SBC_value; q^n = 0
        //output.write("iSign_kpar_wcoll: [3d, 2], SBC\n");
        if (n == 0)
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 2, SBC_value);
        else
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
      } else if (bndry_c == 3) {
        // q^0: Neumann; q^n = 0
        //output.write("iSign_kpar_wcoll: [3d, 3]\n");
        if (n == 0)
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
        else
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
      } else if (bndry_c == 4) {
        // q^n: Neumann; q^0 = 0
        //output.write("iSign_kpar_wcoll: [3d, 4]\n");
        if (n == 0)
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
        else
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
      } else {
        // unsupported flag
        //output.write("Default used!\n");
        throw BoutException("iSign_kpar_wcoll[3d]: unsupported boundary condition flag!\n");
      }


   }

#if 0
   delete invpar; //-deallocate the memory associated with the pointer
   //invpar->~InvertPar(); //-this should be equivalent to delete
#endif

   return(qout);
}

const Field2D iSign_kpar_wcoll(const Field2D &vin, const Field2D &k0, const int bndry_c, const int N, const Field2D &SBC_value)
{

   Field2D qout = 0.0;
   Field2D ddyvin;

   // ddyvin = DDY(vin);
   ddyvin = Grad_par(vin);
   mesh->communicate(ddyvin);

   // InvertPar *invpar;
   // invpar = InvertPar::Create();
   // much better if static:

   static InvertPar* invpar = InvertPar::Create();

   //-range of Lorentzian terms
   const int nMin=0;
   const int nMax= (N == 0 ? 3 : N);

   // set coefficient for Lorentzian fitting
   const BoutReal *alpha;
   const BoutReal *beta;

   // N = 3:
   // Old coef. from  M.V. Umansky, et. al.,  Journal of Nuclear Materials 463, 506 (2015)
   const BoutReal alpha0[] = {0.0018, 0.0769, 2.4498};
   const BoutReal beta0[] = {0.1192, 0.4913, 2.1495};
   // new coef
   // N = 3:
   const BoutReal alpha3[] = {0.01315, 0.924, 14.1365};
   const BoutReal beta3[] = {0.2044, 1.3587, 8.9643};
   // N = 7:
   const BoutReal alpha7[] = { \
            0.007438, 0.6161, 5.9804, 37.9822, 234.3654, \
            1466.4331, 14981.4634};
   const BoutReal beta7[] = {0.1678, 1.1106, 5.6457, 33.1536, \
            202.738, 1254.2144, 9275.3323};
   // N = 12:
   const BoutReal alpha12[] = {0.001424, 0.20736, 2.5653, 14.927, \
             79.3050, 419.2399, 2215.7233, 11709.7857, \
             61885.2763, 327392.6096, 1773350.1566, 16903628.3745};
   const BoutReal beta12[] = {0.09419, 0.6741, 2.9628, 14.43958, \
             75.1106, 395.8293, 2090.8877, 11049.1471, \
             58392.0969, 308695.7371, 1645460.1472, 10794779.4293};
   const BoutReal alpha20[] = {
        0.000550164043066914, 0.10275855014401332, 1.5574968086648138, \
        8.981597596447177, 44.21499619448371, 215.21860764463176, \
        1046.9746330430294, 5093.378441306047, 24777.417109831247, \
        120516.8713406688, 586165.152129492, 2851132.6549336165, \
        13867844.765921066, 67444549.79120941, 327987870.7244059, \
        1595138256.7227132, 7758815048.312665, 37807952481.99112, \
        190502487382.14917, 1778418212793.4631};
   const BoutReal beta20[] = {
        0.0680957687739,  0.506661589811,  2.11683114157,  \
        9.29868746793,  44.186253098,  213.959779668,  \
        1039.85190758,  5057.62188658,  24603.1598059,  \
        119676.949309,  582095.153209,  2831266.12002,  \
        13771369.0447,  66979787.6208,  325739140.621,  \
        1584151555.15,  7704746851.2,  37497152109.4,  \
        1.84685042437e+11,  1.14854367316e+12};
   switch (N) {
     case 0: {
       alpha = alpha0;
       beta = beta0;
       break;
     }
     case 3: {
       alpha = alpha3;
       beta = beta3;
       break;
     }
     case 7: {
       alpha = alpha7;
       beta = beta7;
       break;
     }
     case 12: {
       alpha = alpha12;
       beta = beta12;
       break;
     }
     case 20: {
       alpha = alpha20;
       beta = beta20;
       break;
     }
     default: {
       throw BoutException(
         "ERROR: Invalid choice of number of Lorentzian terms. Must be in [3, 7, 12]\n");
       break;
     }
   }

   #ifdef DEBUG
     output.write("iSign_kpar_wcoll[2d]: bndry_c = %d, N = %d\n", bndry_c, N);
   #endif

   for (int n=nMin; n<nMax; n++)
      //for (int n=1; n<=1; n++) /*just for testing*/
   {

      //-invert parallel Poisson operator: rho = A*phi + B * Grad2_par2(phi)
      invpar->setCoefA(beta[n] * beta[n]); //-free term
      invpar->setCoefB(-1.0 / (k0 * k0));  //-d2/dz2 term

      // TODO: different bndry flag like collisional closure
      // NOTE: ERROR: Field2d += Field3d

      /*
      switch(bndry_c) {
        case 0: {
        output.write("iSign_kpar_wcoll[2d, 0]: bndry_c = %d, N = %d\n", bndry_c, N);
          qout += invpar->solve(alpha[n] / k0 * ddyvin);
          break;
        }
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
          break;
        case 1:
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
          break;
        case 2:
          break;
        case 3:
          if (n == 0)
            qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
          else
            qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
          break;
        defaut: {
          output.write("iSign_kpar_wcoll[2d, df]: bndry_c = %d, N = %d\n", bndry_c, N);
          throw BoutException("iSign_kpar_wcoll: unsupported boundary condition flag!\n");
          break;
        }
      } */

      if (bndry_c == 1) {
        // Default: Zero Parallel Gradient
        qout += invpar->solve(alpha[n] / k0 * ddyvin);
        //output.write("iSign_kpar_wcoll[2d, 0]: bndry_c = %d, N = %d\n", bndry_c, N);
      /* } else if (bndry_c == 0) {
        // Dirichlet
        qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
      } else if (bndry_c == 2) {
        // Sheath Boundary Condition
        throw BoutException("NotImplementedError: iSign_kpar_wcoll[2d] bndry_c == 2!\n");
      } else if (bndry_c == 3) {
        // q^0: Neumann; q^n = 0
        if (n == 0)
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 0);
        else
          qout += invpar->solve(alpha[n] / k0 * ddyvin, 1);
      */
      } else {
        // unsupported flag
        //output.write("iSign_kpar_wcoll[2d, df]: bndry_c = %d, N = %d\n", bndry_c, N);
        throw BoutException("iSign_kpar_wcoll[2d]: unsupported boundary condition flag!\n");
      }
   }

#if 0
   delete invpar; //-deallocate the memory associated with the pointer
   //invpar->~InvertPar(); //-this should be equivalent to delete
#endif

   return(qout);
}

