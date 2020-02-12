
/*******************************************************************
 * qpar with collisions and different boundary condition
 *
 * J.G.Chen, 10/13/2018
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <interpolation.hxx>
#include <derivs.hxx>
#include <math.h>
#include <msg_stack.hxx>
#include <integrops.hxx>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <unistd.h>

// based on qparwcoll.cxx in bout_glf repo by C.H. Ma
const char CXXVERSION[] = "v0.1_03042018_J.G.Chen";

// Evolving variables
Field3D N, P; // Density, Pressure
Vector3D V;   // velocity

// parameters
BoutReal gamma_ratio;   // Ratio of specific heats
BoutReal nu;            // Viscosity
bool include_viscosity;
Field3D a;
Field2D te, te1, ddyterg;
Field3D ddyte;
Field2D vthe, lambdae, k0, qpar_Brag;
Field2D nue, lmfp;
BoutReal n0, logLam, R0, q95, alpha;
Field2D qpar_fs, qpar_eff, qpar_fs2, qpar_eff2; 
Field3D te3d, k03d, qpar_nonlocal_3d_wocoll, qpar_landau_3d;
Field2D qsbc, Cs;   // heat flux sheath boundary condition: qsbc = gamma * n0 * T * Cs
Field3D qsbc3d, Cs3d;
Field2D qpar_nonlocal_0, qpar_nonlocal_1, qpar_nonlocal_3, qpar_nonlocal_4;
Field3D qpar_nonlocal_3d_0, qpar_nonlocal_3d_1, qpar_nonlocal_3d_2, qpar_nonlocal_3d_3, qpar_nonlocal_3d_4;
Field3D qpar_nonlocal_3d;
Field3D qpar_landau_3d_0;
Field3D qpar_landau_3d_1;
Field3D qpar_landau_3d_2;
Field3D qpar_landau_3d_3;
Field3D qpar_landau_3d_4;

int physics_init(bool restarting)
{
   output.write("Parallel heat flux test\n");
   output.write("\tFile     : %s\n", __FILE__);
   output.write("\tCXXVERSION : %s\n", CXXVERSION);
   output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

   //////////////////////////////////////
   // 2D initial profiles

   logLam = 16.;
   std::clock_t tstart;

   // Read initial conditions

//   mesh->get(te, "te");
   mesh->get(te, "Teexp");
   mesh->get(te1,"Tiexp"); // use tiexp for perturbed temp input
   dump.add(te, "te", 0);
   dump.add(te1, "te1", 0);
   te=te+te1;
   output.write("Te = %e -> %e [eV]\n", min(te), max(te));
   te.applyBoundary("neumann");
   mesh->communicate(te);
//   dump.add(te, "te", 0);

   // read options

   Options *options = Options::getRoot();
   options = options->getSection("gas");
   options->get("density", n0, 1e13);

   // used in free streaming expression
   options->get("q95", q95, 4.0);
   options->get("R0", R0, 167.);        // cm
   options->get("alpha", alpha, 1.0);   // flux limiting coefficient

   // calculating intermidiate parameters

   ddyterg = 1.6e-12 * DDY(te);
   ddyterg.applyBoundary("neumann");
   mesh->communicate(ddyterg);
   SAVE_ONCE(ddyterg);
   nue = 2.91e-6 * n0 * logLam / (te * sqrt(te)); //-[1/s], electron coll. frequency
   SAVE_ONCE(nue);
   vthe = 4.19e7 * sqrt(te);          //-[cm/s]
   SAVE_ONCE(vthe);
   lmfp = vthe/nue; //-[cm], electron collisional mfp
   SAVE_ONCE(lmfp);
   k0 = 0.5 / lmfp;

   // calculating qpar
   qpar_Brag = - 3.2 * n0 * vthe * lmfp * ddyterg;
   qpar_Brag.applyBoundary();
   mesh->communicate(qpar_Brag);
   SAVE_ONCE(qpar_Brag);
   qpar_fs = - alpha * n0 * vthe * q95 * R0 * ddyterg;
   qpar_fs.applyBoundary();
   mesh->communicate(qpar_fs);
   SAVE_ONCE(qpar_fs);
   qpar_eff = (qpar_Brag * qpar_fs) / (qpar_Brag + qpar_fs + 1e-10);
   SAVE_ONCE(qpar_eff);

   qpar_fs2 = 1.6e-12 * n0 * vthe * te;
   mesh->communicate(qpar_fs2);
   SAVE_ONCE(qpar_fs2);
   qpar_eff2 = (qpar_Brag * qpar_fs2) / (qpar_Brag + qpar_fs2 + 1e-10);
   SAVE_ONCE(qpar_eff2);

   // heat flux sheath boundary condition
   // ion sound speed, [cm/s]
   // Cs = 9.79e5 * sqrt(gamma * Z * Te/ mu) [cm/s] ---- NRL formular
   Cs = 9.79e5 * sqrt(te);
   qsbc = 4.8 * n0 * (te * 1.6e-12) * Cs;
   SAVE_ONCE(qsbc);
   qsbc = qsbc / (- 1.6 * n0 * vthe);
//   SAVE_ONCE(qsbc);
   // SAVE_ONCE(qsbc3d);

   output.write("Flag1\n");
   tstart = std::clock();
   qpar_nonlocal_0 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te1, k0, 1, 20);
   output.write("2d, 0: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   /* tstart = std::clock();
   output.write("Flag2\n");
   qpar_nonlocal_1 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te1, k0, 1, 20);
   output.write("N = 12: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   tstart = std::clock();
   output.write("Flag3\n");
   qpar_nonlocal_3 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te1, k0, 3, 20);
   output.write("N = 20: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   qpar_nonlocal_4 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te1, k0, 4, 20);
   */
   SAVE_ONCE(qpar_nonlocal_0);
   /* SAVE_ONCE(qpar_nonlocal_1);
   SAVE_ONCE(qpar_nonlocal_3);
   SAVE_ONCE(qpar_nonlocal_4);
   */

   te3d = te1;
   te3d.applyBoundary("neumann");
   mesh->communicate(te3d);
   SAVE_ONCE(te3d);
   k03d = min(k0, true);
   k03d.applyBoundary("neumann");
   mesh->communicate(k03d);
   SAVE_ONCE(k03d);
   output.write("k0 = %e --> %e\n", min(k0), max(k0));
   output.write("k03d = %e --> %e\n", min(k03d), max(k03d));
   ddyte = Grad_par(te3d);
   mesh->communicate(ddyte);
   SAVE_ONCE(ddyte);

   qpar_landau_3d_0 = - 1.6 * n0 * vthe * iSign_kpar(1.6e-12 * te3d, 6.28 / 200000., 0, 20);
   qpar_landau_3d_1 = - 1.6 * n0 * vthe * iSign_kpar(1.6e-12 * te3d, 6.28 / 200000., 1, 20);
   qpar_landau_3d_2 = - 1.6 * n0 * vthe * iSign_kpar(1.6e-12 * te3d, 6.28 / 200000., 2, 20, qsbc);
   qpar_landau_3d_3 = - 1.6 * n0 * vthe * iSign_kpar(1.6e-12 * te3d, 6.28 / 200000., 3, 20);
   qpar_landau_3d_4 = - 1.6 * n0 * vthe * iSign_kpar(1.6e-12 * te3d, 6.28 / 200000., 4, 20);
   qpar_landau_3d_0.applyBoundary("neumann");
   qpar_landau_3d_1.applyBoundary("neumann");
   qpar_landau_3d_2.applyBoundary("neumann");
   qpar_landau_3d_3.applyBoundary("neumann");
   qpar_landau_3d_4.applyBoundary("neumann");
   mesh->communicate(qpar_landau_3d_0);
   mesh->communicate(qpar_landau_3d_1,qpar_landau_3d_2,qpar_landau_3d_3,qpar_landau_3d_4);

   output.write("Flag4\n");
   qpar_nonlocal_3d = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d);
   output.write("3d, default: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   tstart = std::clock();
   qpar_nonlocal_3d_0 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d, 0, 20);
   output.write("3d, 0: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   tstart = std::clock();
   qpar_nonlocal_3d_1 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d, 1, 20);
   qpar_nonlocal_3d_2 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d, 2, 20, qsbc);
   output.write("3d, 1: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   tstart = std::clock();
   qpar_nonlocal_3d_3 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d, 3, 20);
   output.write("3d, 3: %f\n", ( std::clock() - tstart ) / (double) CLOCKS_PER_SEC);
   tstart = std::clock();
   qpar_nonlocal_3d_4 = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d, 4, 20);

   qpar_nonlocal_3d.applyBoundary("neumann");
   qpar_nonlocal_3d_0.applyBoundary("neumann");
   qpar_nonlocal_3d_1.applyBoundary("neumann");
   qpar_nonlocal_3d_2.applyBoundary("neumann");
   qpar_nonlocal_3d_3.applyBoundary("neumann");
   qpar_nonlocal_3d_4.applyBoundary("neumann");
   mesh->communicate(qpar_nonlocal_3d,qpar_nonlocal_3d_0,qpar_nonlocal_3d_1);
   mesh->communicate(qpar_nonlocal_3d_2,qpar_nonlocal_3d_3,qpar_nonlocal_3d_4);

   SAVE_ONCE(qpar_nonlocal_3d);
   SAVE_ONCE(qpar_nonlocal_3d_0);
   SAVE_ONCE(qpar_nonlocal_3d_1);
   SAVE_ONCE(qpar_nonlocal_3d_2);
   SAVE_ONCE(qpar_nonlocal_3d_3);
   SAVE_ONCE(qpar_nonlocal_3d_4);

//   output.write("Flag5\n");
//   qpar_nonlocal_3d_wocoll = - 1.6 * n0 * vthe * iSign_kpar_wcoll(1.6e-12 * te3d, k03d * 0.05);
//   SAVE_ONCE(qpar_nonlocal_3d_wocoll);

   SAVE_ONCE(qpar_landau_3d_0);
   SAVE_ONCE(qpar_landau_3d_1);
   SAVE_ONCE(qpar_landau_3d_2);
   SAVE_ONCE(qpar_landau_3d_3);
   SAVE_ONCE(qpar_landau_3d_4);

//   return 1;

   usleep(1000000);

   output.write("before solve_for\n");
   SOLVE_FOR(a);

   return 0;
}

int physics_run(BoutReal t)
{
//   output.write("enter physics_run\n");
   ddt(a) = 1.;

   return 0;
}


