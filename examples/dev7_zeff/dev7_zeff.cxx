/*******************************************************************************
 * High-Beta Flute-Reduced MHD with 6-field of (N_i, T_e, T_i, U, Psi, Vipar)
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * diffusion_par can open the parallel thermal conductivity
 * T.Y. Xia
 *******************************************************************************/

/* updated by BZhu 02/12/2020
 * 1. perturbed magnetic field is calculated with Apar instead of Psi
 *    note that emass=true is not implemented yet
 *    02/12/2020 -- psi option recovered
 * 2. improved resistivity model
 * 3. corrected parallel and gyro-viscosity terms
 * 4. add hyper-diffusion terms for all six variables
 * 5. flux driven source
 * 	fixed amplitude one, PI controller TBD
 * 6. add zonal field solver
 * 	two needed, one for phi and one for Apar(if electron inertia is kept)
 * 7. add neutral model (experimental)
 * */

// based on elm_6f_landau.cxx v0.2.3_100417_J.G.Chen"
const char CXXVERSION[] = "v0.2.2_02282018_J.G.Chen";

// debug mode to output flags
#define DEBUG_6F 0

#include "bout.hxx"
#include "derivs.hxx"
#include "initialprofiles.hxx"
#include "integrops.hxx"
#include "interpolation.hxx"
#include "invert_laplace.hxx"
#include "invert_parderiv.hxx"
#include "sourcex.hxx"
//#include "neutral.hxx"
#include <boutmain.hxx>
#include <dcomplex.hxx>
#include <fft.hxx>
#include <math.h>
#include <msg_stack.hxx>
#include <invert_laplace.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/invert/laplacexy2.hxx>

/********** Physical constants ************************/
const BoutReal PI = 3.14159265;
const BoutReal MU0 = 4.0e-7 * PI;
const BoutReal Mp = 1.6726e-27;   // [kg] proton mass
const BoutReal ratio_pe = 1836.2; // proton/electron mass ratio
const BoutReal KB = 1.38065e-23;  // Boltamann constant
const BoutReal ee = 1.602e-19;    // ln(Lambda)
const BoutReal eV_K = 11605.0;    // 1eV = 11605K

/********** Magnetic configuration ************************/
int mag_config; // magnetic geometry: 1-circular, 2-circular with limiter, 3-single null, 4-double null
Field2D Rxy, Bpxy, Btxy, B0, hthe, I; // I: shear factor
BoutReal ixsep, ixsep2;
BoutReal jysep1, jysep2, jysep1_2, jysep2_1; // index for x-point on y direction
Vector2D B0vec; // B0 field vector

/********** Primary variables ************************/
BoutReal AA, Zi, Mi, Zeff;              // main ion info (atomic mass, charge and  mass)

Field2D J0, P0; // Current and pressure
BoutReal J0_factor, P0_factor;
Vector2D b0xcv; // Curvature term
Field2D phi0;   // When diamagnetic terms used

Field2D N0, Ti0, Te0, Ne0, N_imp0, T_imp0; // number density and temperature
Field2D Pi0, Pe0, P_imp0, density_tmp;
Field2D q95;
BoutReal q95_input;
bool local_q;
bool n0_fake_prof, n0_p0_0p3, T0_fake_prof, Nimp_lowlimit, quasi_neutral_Ni;
//BoutReal LnLambda; // ln(Lambda)
Field2D LnLambda;

// V0 field vectors
Vector2D Ve0;     // equilibrium ExB velocity: Ve0 = Ve0_net + Ve0_dia
Vector2D Ve0_net; // net flow
Vector2D Ve0_dia; // diamagnetic drift velocity.
Vector2D V0vec;   // V0 field vector in convection
Vector2D V0eff;   // effective V0 field vector in Ohm's law
Vector3D Vexb;    // total ExB velocity
Vector3D Vbtilde; // perturbed B vec: Vbtilde = -b0vec cross Grad Psi

// Er0 field vectors
Vector2D Er0;        // total equilibrium Er0
BoutReal Er0_factor; // Er0 *= Er0_factor
Vector2D Er0_dia;
Vector2D Er0_net; // Er0 = Er0_net + Er0_dia

// 3D evolving variables
Field3D Ni, Te, Ti, P, Pi, Pe;
Field3D U, Vipar, Vepar, Apar, Psi;
// derived variables
Field3D Jpar, phi; // Parallel current, electric potential
Field2D phiDC,VortDC,TeDC;

Field3D ubyn;

Field3D Jpar2;                         // Delp2 of Parallel current
Field3D tmpU2, tmpA2;                  // Grad2_par2new of Parallel vector potential
Field3D tmpN2, tmpTi2, tmpTe2, tmpVp2; // Grad2_par2new of Parallel density

Field3D nu_e, nu_i;         // Electron/ion collision frequency profile (1 / S)
Field3D vth_i, vth_e;       // Electron/ion Thermal Velocity profile (M / S)
Field3D kappa_par_i;        // Ion Thermal Conductivity profile (kg*M / S^2)
Field3D kappa_par_e;        // Electron Thermal Conductivity profile (kg*M / S^2)
BoutReal kappa_par_i_const, kappa_par_e_const;
Field2D omega_ci, omega_ce; // cyclotron frequency
Field3D kappa_perp_i;       // Ion perpendicular Thermal Conductivity profile (kg*M / S^2)
Field3D kappa_perp_e;       // Electron perpendicular Thermal Conductivity profile (kg*M / S^2)
Field3D kappa_par_i_fl, kappa_par_e_fl;  // flux-limited conductivity
Field3D kappa_perp_i_fl, kappa_perp_e_fl;
Field3D kappa_par_i_sp, kappa_par_e_sp; // Spitzer-Harm conductivity
Field3D q_par_i, q_par_e;
Field3D q_par_fl, q_par_landau;

/********** Model options and additional variables ***************/
bool evolve_psi;
bool emass;
static Field2D acoeff; // for emass
static bool aset = false;
Field3D Ajpar; // Parallel current, electric potential
BoutReal emass_inv; // inverse of electron mass
BoutReal coef_jpar;
BoutReal delta_e;     // Normalized electron skin depth
BoutReal delta_e_inv; // inverse normalized electron skin depth
BoutReal gyroAlv;     // Normalized ion current coef

bool diamag;
BoutReal dia_fact; // Multiply diamagnetic term by this

bool energy_flux, energy_exch; // energy flux term
bool diamag_phi0;              // Include the diamagnetic equilibrium phi0
bool thermal_force;            // Include the thermal flux term in Ohm's law
bool eHall;
bool diff_par_flutter;

bool experiment_Er, KH_term; // read in phi_0 from experiment
Field2D phi0_net, U0_net;    // calculate Kelvin-Helmholtz term
Field2D V0, Dphi0;           // net flow amplitude, differential potential to flux
bool diamag_er;              // switch phi0 to Er


int phi_flags, apar_flags;
bool nonlinear;
bool evolve_jpar;
BoutReal g; // Only if compressible
bool phi_curv;


bool include_curvature, include_jpar0, compress0, compresse;
bool include_vipar;
bool evolve_pressure, continuity;
bool parallel_viscous;
Field3D eta_i0, pi_ci;
bool gyroviscous;
Field3D Dperp2Phi0, Dperp2Phi, GradPhi02, GradPhi2; // Temporary variables for gyroviscous
Field3D GradparPhi02, GradparPhi2, GradcPhi, GradcparPhi;
Field3D Dperp2Pi0, Dperp2Pi, bracketPhi0P, bracketPhiP0, bracketPhiP;

BoutReal laplace_alpha; // test the effect of first order term of invert Laplace function

// Bootsctrap current
bool BScurrent;
bool radial_diffusion;
Field3D Jpar_BS0, nu_estar, nu_istar, ft;
Field3D f33;
Field3D cond_neo;
Field3D diff_radial, ddx_ni, ddx_n0;
BoutReal diffusion_coef_Hmode0, diffusion_coef_Hmode1;

// neoclassical effects
bool neoclassic_i, neoclassic_e;
bool neo_resist, ft_simple;
BoutReal major_radius, minor_radius, epsilon;
Field3D xii_neo, xie_neo, Dri_neo, rho_i, rho_e, tmpddx2;
Field3D heatf_neo_i, heatf_neo_e, partf_neo_i;

// resistivity
bool spitzer_resist;              // Use Spitzer formula for resistivity
BoutReal vac_lund, core_lund;     // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
BoutReal vac_resist, core_resist; // The resistivities (just 1 / S)
Field3D eta;                      // Resistivity profile (1 / S)
BoutReal FZ;                      // correction coefficient in Spitzer model

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

// parallel heat flux
// case1: flux limited expression
bool fluxlimit;        // flux limited condition
BoutReal q_alpha;      // flux limiting coefficient
// case2: Landau damping closure
bool Landau;           // Use Gyro Landau Fluid closure instead of thermal conductivity
bool Landau_coll;      // collisional Landau Damping
BoutReal Landau_coeff; // Coefficient for Landau Damping
int nLorentzian;       // number of Lorentzians, collisional: [3, 7, 12], collisionless: >=7
BoutReal kappa_0;      // collisionless Landau damping coefficient
Field3D kappa_i;       // ion collisional Landau damping coefficient (0.5*nu_i/vth_i)
Field3D kappa_e;       // electron collisional Landau damping coefficient (0.5*nu_e/vthe_e)
Field3D SBC_value_i,SBC_value_e;
bool full_sbc;         // full-f version of sheath BC for ion parallel velocity

// impurity
bool impurity_prof, load_impurity, impurity_gyro;
BoutReal Z_imp, A_imp; // impurity ion info
Field3D Dperp2Pimp0, bracketPhiPimp0;
Field2D Upara_imp;

// source
bool source;
Field3D NiSource,TeSource,TiSource;  // Axisymmetric 2D/3D sources
BoutReal NiAmp,TeAmp,TiAmp; // Amplitude of the Gaussian shape sources
int NiLoc,TeLoc,TiLoc,NiSig,TeSig,TiSig;     // Center locations and standard deviation of the sources

// neutral
bool neutral;
Field3D Nn,lNn,Vn,Pn;
Field3D Dn,etan,Sn,Sn_ext,Sv,S_tmp;
Field3D nu_iz,nu_cx,nu_rc,sigma_cx;
BoutReal NnAmp,NnLoc,NnSig;  // used for initialize neutral profile
BoutReal Rcyc_Nn, Rcyc_Vn;  // recycle coefficients

// parallel and perpendicular hyperdiffusion
BoutReal hyperdiff_par_n4, hyperdiff_par_ti4, hyperdiff_par_te4; // M: 4th Parallel density diffusion
BoutReal hyperdiff_par_v4, hyperdiff_par_apar4, hyperdiff_par_u4; // xqx: parallel hyper-viscous diffusion for vorticity
BoutReal hyperdiff_perp_n4, hyperdiff_perp_ti4, hyperdiff_perp_te4; // M: 4th Perpendicular density diffusion
BoutReal hyperdiff_perp_v4, hyperdiff_perp_apar4, hyperdiff_perp_u4;
BoutReal viscos_par;  // Parallel viscosity
BoutReal viscos_perp; // Perpendicular viscosity
BoutReal hyperviscos; // Hyper-viscosity (radial)
Field3D hyper_mu_x;   // Hyper-viscosity coefficient

// position filter
bool pos_filter, pos_filter2, keep_zonalPF;
BoutReal filter_position_ni, filter_position_ti, filter_position_te;
int position_tmpi, position_tmpe, position_tmp;
bool pos_filter_zf;
BoutReal pos_sink_zf, Grid_NX, Grid_NY;
BoutReal pos_filter_width, pos_filter_length, sink_pos_zf;
Field2D pos_sink_ti, pos_sink_te, pos_sink_ni;

// filter low/high-n mode
bool filter_z, filter_z_nonlinear;
int filter_z_mode;
int low_pass_z;
int zonal_flow;
LaplaceXY2 *lapDC;
int zonal_field;
int zonal_bkgd;

/********** Normalization coefficients ************************/
BoutReal density;              // density normalization factor [m^-3]
BoutReal density_unit;         // density unit for grid [m^-3]
BoutReal Bbar, Lbar, Tbar, Va; // Normalization constants
BoutReal Nbar, Tibar, Tebar, Tau_ie;
// coefficients in the equations
BoutReal Psipara1, Upara0, Upara1;
BoutReal Upara2, Upara3, Nipara1;
BoutReal Tipara1, Tipara2, Tipara3;
BoutReal Tepara1, Tepara2, Tepara3, Tepara4;
BoutReal Vepara, Vipara;
// misc
BoutReal Vt0; // equilibrium toroidal flow normalized to Alfven velocity
BoutReal Vp0; // equilibrium poloidal flow normalized to Alfven velocity

/**** Numerical/solver realted variables ************************/
int jx, jy, jz, ncz;                         // index varriable
int xind, indx, indy;
BoutReal dindx;
Field3D fp, fm; // Interpolated on + and - y locations
Field2D F2D_tmp;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm_exb, bm_mag; // Bracket method for advection terms
int bracket_method_exb, bracket_method_mag;

// For preconditioner
Field3D ni_tmp, ti_tmp, te_tmp, vi_tmp, psi_tmp, u_tmp, p_tmp, jpar1, phi_tmp;
Field3D u_tmp1, ti_tmp2, te_tmp2;
Field2D kappa_par_i_lin, kappa_par_e_lin; // for preconditioner
BoutReal Te_tmp1, Ti_tmp1, N_tmp1;

// Fake run related
bool fakerun;
int timestep;
string path;
Vector3D Btilde;

bool limit_jacobi; // limit the infinity value of jacobi at x-point
BoutReal bpxy_constraint, const_bp;
BoutReal hthe_constraint;
bool PF_limit; // filter the instability in PF region
BoutReal PF_limit_range;
BoutReal PF_sink, PFs_width, PFs_length; // sink at inner boundary of PF

bool relax_j_vac;
BoutReal relax_j_tconst; // Time-constant for j relax
Field3D Psitarget;       // The (moving) target to relax to

bool smooth_j_x;              // Smooth Jpar in the x direction
bool mask_j_x, mask_phi_x;    // Mask Jpar at the inner boundary of x
Field3D mask_jx1d, mask_px1d; // the variable of mask function, normalized to 1.
int mask_flag_j, mask_flag_phi;
BoutReal mask_width, mask_length;
BoutReal filter_nl;

int jpar_bndry_width; // Zero jpar in a boundary region

bool parallel_lr_diff;  // Use left and right shifted stencils for parallel differences
bool parallel_lagrange; // Use (semi-) Lagrangian method for parallel derivatives
bool parallel_project;  // Use Apar to project field-lines

// for debug purpose
Field3D term1, term2, term3, term4, term5;

/****************************************************************/

/*****************************************************************/
/*****************************************************************/

BoutReal n0_height, n0_ave, n0_width, n0_center, n0_bottom_x; // the total height, average width and center of profile of N0
BoutReal Tconst;                                              // the ampitude of congstant temperature

// Vector2D bxgradb;

// Constraint
Field3D C_phi;

// Parameters

BoutReal diffusion_par;  // Parallel thermal conductivity is multiplied by this value. (<0 off)
BoutReal diffusion_perp; // Perpendicular thermal conductivity (>0 open)

// seems like obsolete stuff?
BoutReal heating_P; // heating power in pressure
BoutReal hp_width;  // heating profile radial width in pressure
BoutReal hp_length; // heating radial domain in pressure
BoutReal sink_vp;   // sink in pressure
BoutReal sp_width;  // sink profile radial width in pressure
BoutReal sp_length; // sink radial domain in pressure

BoutReal sink_Ul;    // left edge sink in vorticity
BoutReal su_widthl;  // left edge sink profile radial width in vorticity
BoutReal su_lengthl; // left edge sink radial domain in vorticity

BoutReal sink_Ur;    // right edge sink in vorticity
BoutReal su_widthr;  // right edge sink profile radial width in vorticity
BoutReal su_lengthr; // right edge sink radial domain in vorticity

BoutReal sink_Tel;    // left edge sink in Te
BoutReal ste_widthl;  // left edge sink profile radial width in Te
BoutReal ste_lengthl; // left edge sink radial domain in Te

BoutReal sink_Ter;    // right edge sink in Te
BoutReal ste_widthr;  // right edge sink profile radial width in Te
BoutReal ste_lengthr; // right edge sink radial domain in Te

BoutReal sink_Psir;   // right edge sink in Psi
BoutReal spsi_widthr; // right edge sink profile radial width in Psi/Apar
BoutReal spsi_lengthr;// right edge sink radial domain in Psi/Apar

BoutReal Low_limit; // To limit the negative value of total density and temperatures

Field3D Te_tmp, Ti_tmp, N_tmp, Ne_tmp; // to avoid the negative value of total value
BoutReal gamma_i_BC, gamma_e_BC;       // sheath energy transmission factors
int Sheath_width;
bool SBC_phi;
Field3D c_se, Jpar_sh, q_se, q_si, vth_et, c_set, phi_sh; // variables for sheath boundary conditions
Field2D vth_e0, c_se0, Jpar_sh0, phi_sh0;
BoutReal const_cse;

//********************

Field3D Xip_x, Xip_z; // Displacement of y+1 (in cell index space)

Field3D Xim_x, Xim_z; // Displacement of y-1 (in cell index space)

bool phi_constraint; // Solver for phi using a solver constraint

bool output_transfer; // output the results of energy transfer
bool output_ohm;      // output the results of the terms in Ohm's law
bool output_flux_par; // output the results of parallel particle and heat flux
bool output_vradial;  // output the results of radial velocity, induced by ExB and magnetic flutter
bool output_Teterms, output_Titerms, output_Tevegradte, output_qparcompare;
bool output_P, output_Vepar, output_SBC;
bool output_eta, output_kappa_par;  // if false, only the initial value is written.

Field3D T_M, T_R, T_ID, T_C, T_G; // Maxwell stress, Reynolds stress, ion diamagbetic and curvature term
Field3D ohm_phi, ohm_hall, ohm_thermal;
Field3D gamma_par_i, heatf_par_i, heatf_par_e; // particle flux, ion and elelctron heat flux
Field3D heatf_par_flutter_i, heatf_par_flutter_e;
Field3D bracket1i, gradpar_ti, bracket1e, gradpar_te; // temp variable for perturbed parallel thermal conduction
Field3D Vbti_par, Vbte_par;

BoutReal hyperresist;  // Hyper-resistivity coefficient (in core only)
BoutReal ehyperviscos; // electron Hyper-viscosity coefficient
Field3D hyper_eta_x;   // Radial resistivity profile
Field3D hyper_eta_z;   // Toroidal resistivity profile

int damp_width;        // Width of inner damped region
BoutReal damp_t_const; // Timescale of damping

// const BoutReal Low_limit = 1.e-10;   // limit of the profile to prevent minus total value

// Communication objects
FieldGroup comms;

// additional functions
int precon(BoutReal t, BoutReal cj, BoutReal delta);

// Functions for sheath boundary conditions
void SBC_Dirichlet(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range);
void SBC_Gradpar(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range);
void SBC_yup_eq(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range);
void SBC_ydown_eq(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range);
void SBC_yup_Grad_par(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range);
void SBC_ydown_Grad_par(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range);

void advect_tracer(const Field3D &p,                               // phi (input)
                   const Field3D &delta_x, const Field3D &delta_z, // Current location (input)
                   Field3D &F_dx, Field3D &F_dz);                  // Time-derivative of location

const Field3D Grad2_par2new(const Field3D &f); // for 4th order diffusion

const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x);

const Field3D field_larger(const Field3D &f, const BoutReal limit);
const Field2D field_larger(const Field2D &f, const BoutReal limit);

const Field2D Invert_laplace2(const Field2D &f, int flags);

const Field3D BS_ft(const int index);
const Field3D F31(const Field3D input);
const Field3D F32ee(const Field3D input);
const Field3D F32ei(const Field3D input);
const Field3D F33(const Field3D input);

// const Field2D smooth_xy(const Field2D &f, bool BoutRealspace);

const Field3D PF_filter(const Field3D &input, const BoutReal PF_limit_range);
const Field3D sink_PF(const Field2D &f0, const Field3D &f, const BoutReal width, const BoutReal length);

// reload the mask function on x direction, return a normalized function
const Field3D mask_x_1d(bool BoutRealspace, int mask_flag, BoutReal mask_width, BoutReal mask_length);

// reload the filter function for the nonlinear running with one mode and zonal modes
const Field3D filter(const Field3D &var, int N0, int N1);
BoutReal TanH(BoutReal a);

const Field3D lowPass_pos(const Field3D &var, int filter_index);
const Field3D lowPass_pos2(const Field3D &var, const Field3D &prof); // filter the zonal component where is the negative value.
const Field3D sink_zonal_core(const Field3D &f, int filter_index);
const Field2D tanhxl_core(const int filter_index);

// functions
const Field2D tanhxl_core(const int filter_index) {
    Field2D result;
    result.allocate();

#ifdef CHECK
    msg_stack.push("tanhxl_core(int)", filter_index);
#endif

    result = Field2D(0.);
    BoutReal xpos, width, length, tanh_tmp;
    int indy;
    xpos = filter_index;
    width = xpos * pos_filter_width;   // the width of the tanh filter function
    length = xpos * pos_filter_length; // the middle point of the tanh filter function

    // output.write("Y index: %i, %i, %i, %i", int(jysep1), int(jysep2_1), int(jysep1_2), int(jysep2));

    for (jx = 0; jx < mesh->ngx; jx++) {
        indx = mesh->XGLOBAL(jx);
        tanh_tmp = (1. - TanH((indx - length) / width)) / 2.;
        if (tanh_tmp < 0.005)
            tanh_tmp = 0.;
        for (jy = 0; jy < mesh->ngy; jy++) {
            indy = mesh->YGLOBAL(jy);
            if (((indy > int(jysep1)) && (indy <= int(jysep2_1))) || ((indy > int(jysep1_2)) && (indy <= int(jysep2)))) {
                result[jx][jy] = tanh_tmp;
                // output.write("yindex= %i",indy);
                // output.write("%i, %i, %i, %i",jysep1, jysep2_1, jysep1_2, jysep2);
            } else
                result[jx][jy] = 0.;
            // output.write("jx= %d, jy= %d\n",jx,jy);
        }
    }

#ifdef CHECK
    msg_stack.pop();
#endif

    return result;
}

const Field3D lowPass_pos2(const Field3D &var, const Field3D &prof) {
    Field3D tmp_result, result;
    result.allocate();
    tmp_result.allocate();
    Field2D prof2d = prof.DC();
    BoutReal y_ind;

    result = var;

    tmp_result = lowPass(var, low_pass_z, 0);

    for (int jx = 0; jx < mesh->ngx; jx++) {
        for (int jy = 0; jy < mesh->ngy; jy++) {
            y_ind = mesh->YGLOBAL(jy);
            if (keep_zonalPF && (((y_ind > int(jysep1)) && (y_ind <= int(jysep2_1))) || ((y_ind > int(jysep1_2)) && (y_ind <= int(jysep2))))) {
                if (prof2d[jx][jy] < 0.) {
                    for (int jz = 0; jz < mesh->ngz; jz++) {
                        result[jx][jy][jz] = tmp_result[jx][jy][jz];
                    }
                }
            }
        }
    }

    return result;
}

const Field3D lowPass_pos(const Field3D &var, int filter_index) {
    Field3D tmp_result, result;
    result.allocate();
    tmp_result.allocate();
    int y_ind;

    result = var;

    tmp_result = lowPass(var, low_pass_z, 0);

    for (int jx = 0; jx < mesh->ngx; jx++) {
        if (mesh->XGLOBAL(jx) <= filter_index) {
            for (int jy = 0; jy < mesh->ngy; jy++) {
                y_ind = mesh->YGLOBAL(jy);
                if (keep_zonalPF && (((y_ind > int(jysep1)) && (y_ind <= int(jysep2_1))) || ((y_ind > int(jysep1_2)) && (y_ind <= int(jysep2))))) {
                    for (int jz = 0; jz < mesh->ngz; jz++) {
                        result[jx][jy][jz] = tmp_result[jx][jy][jz];
                    }
                }
            }
        }
    }

    return result;
}

const Field3D filter(const Field3D &var, int N0, int N1) {
    // ASSERT1(var.isAllocated());

    static dcomplex *f = (dcomplex *)NULL;

    ncz = mesh->ngz - 1;

    if (f == (dcomplex *)NULL) {
        // Allocate memory
        f = new dcomplex[ncz / 2 + 1];
    }

    Field3D result;
    result.allocate();

    for (jx = 0; jx < mesh->ngx; jx++) {
        for (jy = 0; jy < mesh->ngy; jy++) {

            rfft(var[jx][jy], ncz, f); // Forward FFT

            for (jz = 0; jz <= ncz / 2; jz++) {

                if ((jz != N0) && (jz != N1)) {
                    // Zero this component
                    f[jz] = 0.0;
                }
            }

            irfft(f, ncz, result[jx][jy]); // Reverse FFT

            result[jx][jy][ncz] = result[jx][jy][0];
        }
    }

#ifdef TRACK
    result.name = "filter(" + var.name + ")";
#endif

    // result.location = var.location;

    return result;
}

/*BoutReal TanH(BoutReal a)
{
  BoutReal temp = exp(a);
  return ((temp - 1.0 / temp) / (temp + 1.0 / temp));
}*/

const Field3D mask_x_1d(bool BoutRealspace, int mask_flag, BoutReal mask_width, BoutReal mask_length) {
    Field3D result;
    result.allocate();

    BoutReal Grid_NX; // the grid number on x, and the
    mesh->get(Grid_NX, "nx");

    // create a radial buffer zone to set jpar zero near radial boundary

    BoutReal min_tmp = (TanH((4. / Grid_NX - mask_length) / mask_width) + 1.) / 2.;
    BoutReal max_tmp = (TanH(((1. - (Grid_NX - 5.) / Grid_NX) - mask_length) / mask_width) + 1.) / 2.;

    for (jx = 0; jx < mesh->ngx; jx++)
        for (jy = 0; jy < mesh->ngy; jy++)
            for (jz = 0; jz < mesh->ngz; jz++) {
                BoutReal lx = mesh->GlobalX(jx);
                BoutReal dampl = (TanH((lx - mask_length) / mask_width) + 1.) / 2. - min_tmp;
                BoutReal dampr = (TanH(((1. - lx) - mask_length) / mask_width) + 1.) / 2. - max_tmp;
                if (mask_flag == 0) // left mask
                    result[jx][jy][jz] = dampl;
                else if (mask_flag == 1) // right mask
                    result[jx][jy][jz] = dampr;
                else // mask on both boundary
                    result[jx][jy][jz] = dampl * dampr;
                if (result[jx][jy][jz] < 0)
                    result[jx][jy][jz] = 0.;
            }

    result /= max(result, true);
    if (BoutRealspace)
        result = result.shiftZ(false); // Shift back

    // Need to communicate boundaries
    mesh->communicate(result);

    return result;
}

const Field2D Invert_laplace2(const Field2D &f, int flags) {
    Field3D f_tmp, result_tmp;
    Field2D result;
    f_tmp.allocate();
    result_tmp.allocate();
    result.allocate();

    f_tmp = f;

    result_tmp = invert_laplace(f_tmp, flags, NULL);
    mesh->communicate(result_tmp);
    result_tmp = smooth_x(result_tmp);
    result_tmp = nl_filter_y(result_tmp, 1);

    for (jx = 0; jx < mesh->ngx; jx++)
        for (jy = 0; jy < mesh->ngy; jy++)
            result[jx][jy] = result_tmp[jx][jy][0];

    return result;
}

const Field3D field_larger(const Field3D &f, const BoutReal limit) {
    Field3D result;
    result.allocate();

    //  #pragma omp parallel for
    for (jx = 0; jx < mesh->ngx; jx++)
        for (jy = 0; jy < mesh->ngy; jy++)
            for (jz = 0; jz < mesh->ngz; jz++) {
                if (f[jx][jy][jz] >= limit)
                    result[jx][jy][jz] = f[jx][jy][jz];
                else
                    result[jx][jy][jz] = limit;
            }
    mesh->communicate(result);
    return result;
}

const Field2D field_larger(const Field2D &f, const BoutReal limit) {
    Field2D result;
    result.allocate();

    //  #pragma omp parallel for
    for (jx = 0; jx < mesh->ngx; jx++)
        for (jy = 0; jy < mesh->ngy; jy++) {
            if (f[jx][jy] >= limit)
                result[jx][jy] = f[jx][jy];
            else {
                result[jx][jy] = 0.9 * limit + 0.1 * f[jx][jy];
            }
        }

    for (jx = 1; jx < mesh->ngx - 1; jx++)
        for (jy = 1; jy < mesh->ngy - 1; jy++) {
            {
                if (f[jx][jy] <= 1.2 * limit) {
                    result[jx][jy] = 0.5 * result[jx][jy] + 0.25 * (result[jx - 1][jy] + result[jx + 1][jy]);
                    result[jx][jy] = 0.5 * result[jx][jy] + 0.25 * (result[jx - 1][jy] + result[jx + 1][jy]);
                    // result[jx][jy] = 0.5*result[jx][jy] + 0.25*( result[jx-1][jy] + result[jx+1][jy] );
                    // result[jx][jy] = 0.5*result[jx][jy] + 0.25*( result[jx-1][jy] + result[jx+1][jy] );
                    // result[jx][jy] = 0.5*result[jx][jy] + 0.25*( result[jx-1][jy] + result[jx+1][jy] );
                    // result[jx][jy] = 0.5*result[jx][jy] + 0.25*( result[jx-1][jy] + result[jx+1][jy] );
                }
            }
        }
    mesh->communicate(result);
    return result;
}

/*const Field2D smooth_xy(const Field2D &f);
{
Field2D fs, result;

fs = f;

result.allocate();

  // Copy boundary region
  for(int jy=0;jy<mesh->ngy;jy++)
  {
  result[0][jy] = fs[0][jy];
  result[mesh->ngx-1][jy] = fs[mesh->ngx-1][jy];
      }

  // Smooth using simple 1-2-1 filter

  for(int jx=1;jx<mesh->ngx-1;jx++)
  for(int jy=0;jy<mesh->ngy;jy++)
  {
  result[jx][jy] = 0.5*fs[jx][jy] + 0.25*( fs[jx-1][jy] + fs[jx+1][jy] );
  result[jx][jy] = 0.5*fs[jx][jy] + 0.25*( fs[jx][jy-1] + fs[jx][jy+1] );
  }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
  }*/

const Field3D PF_filter(const Field3D &input, const BoutReal PF_limit_range) {
    Field3D result;
    result.allocate();

    // BoutReal jysep1, jysep2;
    // mesh->get(jysep1, "jyseps1_1");
    // mesh->get(jysep2, "jyseps2_2");

    for (jx = 0; jx < mesh->ngx; jx++) {
        indx = mesh->XGLOBAL(jx);
        dindx = indx / ixsep;
        for (jy = 0; jy < mesh->ngy; jy++) {
            indy = mesh->YGLOBAL(jy);
            // output.write("dinx: %e   indy: %e  jysep1: %e  jysep2: %e\n", dindx, indy, jysep1, jysep2);
            if ((dindx < PF_limit_range) && ((indy <= int(jysep1)) || ((indy > int(jysep2_1)) && (indy <= int(jysep1_2))) || (indy > int(jysep2)))) {
                for (jz = 0; jz < mesh->ngz; jz++)
                    result[jx][jy][jz] = 0.;
            } else {
                for (jz = 0; jz < mesh->ngz; jz++)
                    result[jx][jy][jz] = input[jx][jy][jz];
            }
        }
    }

    mesh->communicate(result);
    return result;
}

const Field3D sink_zonal_core(const Field3D &var, int filter_index) {
    Field3D result;
    result.allocate();
    static dcomplex *f = NULL, *f2 = NULL;
    int indx, indy;

#ifdef CHECK
    msg_stack.push("sink_zonal_core(Field3D, int)", filter_index);
#endif

    BoutReal xpos, width, length, tanh_tmp;
    xpos = filter_index;
    width = xpos * pos_filter_width;   // the width of the tanh filter function
    length = xpos * pos_filter_length; // the middle point of the tanh filter function

    // output.write("Error zonal_part~\n");
    if (!var.isAllocated()) {
        return var.DC();
    }

    int ncz = mesh->ngz - 1;

    if (f == NULL)
        f = new dcomplex[ncz / 2 + 1];

    if (f2 == NULL)
        f2 = new dcomplex[ncz / 2 + 1];

    for (jx = 0; jx < mesh->ngx; jx++) {
        indx = mesh->XGLOBAL(jx);
        tanh_tmp = (1. + TanH((indx - length) / width)) / 2.;
        if (tanh_tmp < 0.005)
            tanh_tmp = 0.;
        for (jy = 0; jy < mesh->ngy; jy++) {
            indy = mesh->YGLOBAL(jy);
            //	  if ( ((indy > int(jysep1)) && (indy <= int(jysep2_1))) || ((indy > int(jysep1_2)) && (indy <= int(jysep2))) )
            {
                // Take FFT in the Z direction
                rfft(var[jx][jy], ncz, f);
                // Filter the zonal component based on the filter_index
                //	      f[0] *= (1.+TanH( (indx-length)/width ))/2.;
                f[0] *= tanh_tmp;
            }

            irfft(f, ncz, result[jx][jy]); // Reverse FFT
            result[jx][jy][ncz] = result[jx][jy][0];
        }
    }

#ifdef CHECK
    msg_stack.pop();
#endif
    // mesh->communicate(result);
    return result;
}

const Field3D sink_PF(const Field2D &f0, const Field3D &f, const BoutReal width, const BoutReal length) {
    Field3D result;
    result.allocate();

    // BoutReal jysep1, jysep2;
    // mesh->get(jysep1, "jyseps1_1");
    // mesh->get(jysep2, "jyseps2_2");

    result = sink_tanhxl(f0, f, width, length);
    for (jy = 0; jy < mesh->ngy; jy++) {
        indy = mesh->YGLOBAL(jy);
        if (((indy > int(jysep1)) && (indy <= int(jysep2_1))) || ((indy > int(jysep1_2)) && (indy <= int(jysep2)))) {
            for (jx = 0; jx < mesh->ngx; jx++)
                for (jz = 0; jz < mesh->ngz; jz++)
                    result[jx][jy][jz] = 0.;
        } else {
            for (jx = 0; jx < mesh->ngx; jx++) {
                if (mesh->XGLOBAL(jx) >= ixsep)
                    for (jz = 0; jz < mesh->ngz; jz++)
                        result[jx][jy][jz] = 0.;
            }
        }
    }

    mesh->communicate(result);
    return result;
}

const Field3D Grad2_par2new(const Field3D &f) {
/*
 * This function implements d2/dy2 where y is the poloidal coordinate theta
 */

#ifdef CHECK
    int msg_pos = msg_stack.push("Grad2_par2new( Field3D )");
#endif

    Field3D result = D2DY2(f);

#ifdef TRACK
    result.name = "Grad2_par2new(" + f.name + ")";
#endif
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return result;
}

const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x) {
    Field2D result;
    result.allocate();

    BoutReal Grid_NX, Grid_NXlimit; // the grid number on x, and the
    mesh->get(Grid_NX, "nx");

    Grid_NXlimit = n0_bottom_x * Grid_NX;
    output.write("Jysep1_1 = %i   Grid number = %e\n", int(jysep1), Grid_NX);

    if (jysep1 > 0.) // for single null geometry
    {
        // BoutReal Jysep2;
        // mesh->get(Jysep2, "jyseps2_2");
        // output.write("Jysep2_2 = %i   Ixsep1 = %i\n", int(Jysep2), int(Jxsep));

        for (jx = 0; jx < mesh->ngx; jx++) {
            BoutReal mgx = mesh->GlobalX(jx);
            BoutReal xgrid_num = (ixsep + 1.) / Grid_NX;
            // output.write("mgx = %e xgrid_num = %e\n", mgx);
            for (jy = 0; jy < mesh->ngy; jy++) {
                int globaly = mesh->YGLOBAL(jy);
                // output.write("local y = %i;   global y: %i\n", jy, globaly);
                if (mgx > xgrid_num || (globaly <= int(jysep1)) || ((globaly > int(jysep2_1)) && (globaly <= int(jysep1_2))) || (globaly > int(jysep2)))
                    mgx = xgrid_num;
                BoutReal rlx = mgx - n0_center;
                BoutReal temp = exp(rlx / n0_width);
                BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
                result[jx][jy] = 0.5 * (1.0 - dampr) * n0_height + n0_ave;
            }
        }
    } else // circular geometry
    {
        for (jx = 0; jx < mesh->ngx; jx++) {
            BoutReal mgx = mesh->GlobalX(jx);
            BoutReal xgrid_num = Grid_NXlimit / Grid_NX;
            if (mgx > xgrid_num)
                mgx = xgrid_num;
            BoutReal rlx = mgx - n0_center;
            BoutReal temp = exp(rlx / n0_width);
            BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
            for (jy = 0; jy < mesh->ngy; jy++)
                result[jx][jy] = 0.5 * (1.0 - dampr) * n0_height + n0_ave;
        }
    }

    mesh->communicate(result);

    return result;
}

int physics_init(bool restarting) {

    bool noshear;

    output.write("Solving high-beta flute reduced equations\n");
    output.write("\tFile    : %s\n", __FILE__);
    output.write("\tCXXVERSION : %s\n", CXXVERSION);
    output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

    //////////////////////////////////////////////////////////////
    // Load data from the grid

    mesh->get(ixsep, "ixseps1");
    mesh->get(ixsep2, "ixseps2");
    mesh->get(jysep1, "jyseps1_1");
    mesh->get(jysep2, "jyseps2_2");
    mesh->get(jysep1_2, "jyseps1_2");
    mesh->get(jysep2_1, "jyseps2_1");

    mesh->get(Grid_NX, "nx");
    mesh->get(Grid_NY, "ny");
    // output.write("%i, %i, %i, %i,\n",jysep1, jysep2_1, jysep1_2, jysep2);

    if (ixsep==Grid_NX) {
        output.write("Cicular geometry without limiter!\n");
        mag_config = 1;
    } else if (ixsep2<Grid_NX) {
        output.write("Double null geometry!\n");
        mag_config = 4;
    } else if (jysep1<0) {
        output.write("Circular geometry with limiter!\n");
        mag_config = 2;
    } else if (jysep1>0) {
        output.write("Single null geometry!\n");
        mag_config = 3;
    } else {
        output.write("CAUTION: magnetic configuration cannot be determined!\n");
    }

    // Load 2D profiles
    mesh->get(J0, "Jpar0");          // A / m^2
    if (mesh->get(P0, "pressure_s")) // Pascals
    {
        mesh->get(P0, "pressure");
        output.write("Using pressure as P0.\n");
    } else
        output.write("Using pressure_s as P0.\n");

    // Load curvature term
    b0xcv.covariant = false;  // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

    // Load metrics
    if (mesh->get(Rxy, "Rxy")) { // m
        output.write("Error: Cannot read Rxy from grid\n");
        return 1;
    }
    if (mesh->get(Bpxy, "Bpxy")) { // T
        output.write("Error: Cannot read Bpxy from grid\n");
        return 1;
    }
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(B0, "Bxy");    // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I, "sinty");   // m^-2 T^-1

    //////////////////////////////////////////////////////////////
    // Read parameters from the options file
    //
    // Options.get ( NAME,    VARIABLE,    DEFAULT VALUE)
    //
    // or if NAME = "VARIABLE" then just
    //
    // OPTION(VARIABLE, DEFAULT VALUE)
    //
    // Prints out what values are assigned
    /////////////////////////////////////////////////////////////

    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("highbeta");

    OPTION(options, n0_fake_prof,          false);  // use the hyperbolic profile of n0. If both  n0_fake_prof and T0_fake_prof are false, use the profiles from grid file
    OPTION(options, n0_p0_0p3,             false);  // use n0 ~ P0^0.3
    OPTION(options, n0_height,               0.4);  // the total height of profile of N0, in percentage of Ni_x
    OPTION(options, n0_ave,                 0.01);  // the center or average of N0, in percentage of Ni_x
    OPTION(options, n0_width,                0.1);  // the width of the gradient of N0,in percentage of x
    OPTION(options, n0_center,             0.633);  // the grid number of the center of N0, in percentage of x
    OPTION(options, n0_bottom_x,            0.81);  // the start of flat region of N0 on SOL side, in percentage of x
    OPTION(options, T0_fake_prof,          false);
    OPTION(options, Tconst,                 -1.0);  // the amplitude of constant temperature, in percentage
    OPTION(options, impurity_prof,         false);  // Include the profile of impurity in the equilibrium
    OPTION(options, impurity_gyro,         false);  // Include the gyro viscous of impurity
    OPTION(options, load_impurity,         false);  // if load impurity from grid
    OPTION(options, Z_imp,                    6.);  // The charge number of impurity
    OPTION(options, A_imp,                   12.);  // The mass number of impurity
    OPTION(options, Nimp_lowlimit,         false);  // The switch to limit the lowest value of impurity density
    OPTION(options, quasi_neutral_Ni,      false);  // The switch to use quasi neutral condition to calculate N0 instead of reading from grid

    OPTION(options, experiment_Er,         false);
    OPTION(options, Er0_factor,              1.0);  // change Er0 *= Er0_factor
    OPTION(options, KH_term,               false);  // switch to Kelvin-Helmholtz term
    OPTION(options, J0_factor,               1.0);
    OPTION(options, P0_factor,               1.0);

    OPTION(options, laplace_alpha,           1.0);  // test parameter for the cross term of invert Lapalace
    OPTION(options, Low_limit,           1.0e-10);  // limit the negative value of total quantities
    OPTION(options, q95_input,               5.0);  // input q95 as a constant, if <0 use profile from grid
    OPTION(options, local_q,               false);  // using magnetic field to calculate q profile

    OPTION(options, gamma_i_BC,             -1.0);  // sheath energy transmission factor for ion
    OPTION(options, gamma_e_BC,             -1.0);  // sheath energy transmission factor for electron
    OPTION(options, Sheath_width,              1);  // Sheath boundary width in grid number
    OPTION(options, SBC_phi,               false);  // use sheath boundary on phi instead of Jpar

    OPTION(options, density,              1.0e20);  // number density normalization factor [m^-3]
    OPTION(options, density_unit,         1.0e20);  // Number density unit for grid [m^-3]
    OPTION(options, Zi,                        1);  // ion charge number
    OPTION(options, Zeff,                      1);  // Zeff used in resistivity calculation

    OPTION(options, evolve_jpar,           false);  // If true, evolve J raher than Psi
    OPTION(options, phi_constraint,        false);  // Use solver constraint for phi

    // Effects to include/exclude
    OPTION(options, nonlinear,             false);
    OPTION(options, include_curvature,      true);
    OPTION(options, include_jpar0,          true);
    OPTION(options, evolve_pressure,        true);
    OPTION(options, evolve_psi,             true);

    OPTION(options, continuity,            false);  // use continuity equation
    OPTION(options, compress0,             false);
    OPTION(options, compresse,                 false);
    OPTION(options, gyroviscous,           false);
    OPTION(options, parallel_viscous,      false);

    OPTION(options, BScurrent,             false);

    OPTION(options, radial_diffusion,      false);
    OPTION(options, diffusion_coef_Hmode0,   1.0);  // default value of radial diffusion coefficient
    OPTION(options, diffusion_coef_Hmode1,  10.0);  // upper limit of radial diffusion coefficient

    OPTION(options, fakerun,               false);  // Only calculate the transport equation, work if Hmode_rc2 = true
    OPTION(options, path,                   "./");  // The path of the original Vexb data
    // OPTION(options, timestep,                  0);  // The timestep of previous data used as the initial of transport equation

    OPTION(options, pos_filter,            false);  // switch to turn on the filter of the negative value of zonal background
    OPTION(options, pos_filter2,           false);  // switch to turn on the filter inside certain position
    OPTION(options, pos_filter_zf,         false);  // switch to turn on the filter of the dc profiles inside certain postion with tanh function
    OPTION(options, keep_zonalPF,          false);  // keep the zonal component in PF region when zonal filter is turned on
    OPTION(options, filter_position_ni,      100);  // radial index of the filter. Zonal component of Ni in the x range of 0 - filter_position will be filtered.
    OPTION(options, filter_position_ti,      100);  // radial index of the filter. Zonal component of Ti in the x range of 0 - filter_position will be filtered.
    OPTION(options, filter_position_te,      100);  // radial index of the filter. Zonal component of Te in the x range of 0 - filter_position will be filtered.
    OPTION(options, pos_filter_width,       0.02);  // width of po_filter_zf with the normalization of filter_position
    OPTION(options, pos_filter_length,      0.95);  // length of po_filter_zf with the normalization of filter_position
    OPTION(options, pos_sink_zf,            -10.);  // switch to turn on the sink of the dc profiles inside certain postion with tanh function

    if (pos_filter2 || pos_filter_zf || (pos_sink_zf > 0.)) {
        if (filter_position_ni <= 0) {
            filter_position_ni = ixsep;
            output << "\tWarning: filter_position_ni has a negtive value!\n";
            output << "\tfilter_position_ni is forced to be isxeps1=" << ixsep << ".\n";
        }
        if (filter_position_ti <= 0) {
            filter_position_ti = ixsep;
            output << "\tWarning: filter_position_ti has a negtive value!\n";
            output << "\tfilter_position_ti is forced to be isxeps1=" << ixsep << ".\n";
        }
        if (filter_position_te <= 0) {
            filter_position_te = ixsep;
            output << "\tWarning: filter_position_te has a negtive value!\n";
            output << "\tfilter_position_te is forced to be isxeps1=" << ixsep << ".\n";
        }

        position_tmpi = min(filter_position_ni, filter_position_ti);
        position_tmpe = min(filter_position_ni, filter_position_te);
        position_tmp = min(position_tmpi, position_tmpe);
    }

    OPTION(options, NiAmp,             -1.0);  // Amplitude of the explicit particle sourcing
    OPTION(options, TeAmp,             -1.0);  // Amplitude of the explicit electron energy sourcing
    OPTION(options, TiAmp,             -1.0);  // Amplitude of the explicit ion energy sourcing
    OPTION(options, NiLoc,  floor(ixsep/4.));
    OPTION(options, TeLoc,  floor(ixsep/4.));
    OPTION(options, TiLoc,  floor(ixsep/4.));
    OPTION(options, NiSig, floor(ixsep/12.));
    OPTION(options, TeSig, floor(ixsep/12.));
    OPTION(options, TiSig, floor(ixsep/12.));

    OPTION(options, neutral,           false);
    OPTION(options, NnAmp,               0.2);
    OPTION(options, NnLoc,             ixsep);
    OPTION(options, NnSig, floor(Grid_NX/4.));
    OPTION(options, Rcyc_Nn,             1.0);
    OPTION(options, Rcyc_Vn,             1.0);

    // int bracket_method;
    OPTION(options, bracket_method_exb,        0);
    switch (bracket_method_exb) {
        case 0: {
            bm_exb = BRACKET_STD;
            output << "\tBrackets for ExB: default differencing\n";
            break;
        }
        case 1: {
            bm_exb = BRACKET_SIMPLE;
            output << "\tBrackets for ExB: simplified operator\n";
            break;
        }
        case 2: {
            bm_exb = BRACKET_ARAKAWA;
            output << "\tBrackets for ExB: Arakawa scheme\n";
            break;
        }
        case 3: {
            bm_exb = BRACKET_CTU;
            output << "\tBrackets for ExB: Corner Transport Upwind method\n";
            break;
        }
        default:
            output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
            return 1;
    }

    // int bracket_method;
    OPTION(options, bracket_method_mag,        2);
    switch (bracket_method_mag) {
        case 0: {
            bm_mag = BRACKET_STD;
            output << "\tBrackets: default differencing\n";
            break;
        }
        case 1: {
            bm_mag = BRACKET_SIMPLE;
            output << "\tBrackets: simplified operator\n";
            break;
        }
        case 2: {
            bm_mag = BRACKET_ARAKAWA;
            output << "\tBrackets: Arakawa scheme\n";
            break;
        }
        case 3: {
            bm_mag = BRACKET_CTU;
            output << "\tBrackets: Corner Transport Upwind method\n";
            break;
        }
        default:
            output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
            return 1;
    }

    OPTION(options, eHall,                 false);  // electron Hall or electron parallel pressue gradient effects?
    OPTION(options, thermal_force,         false);  // Thermal flux in Ohm's Law
    OPTION(options, AA,                      2.0);  // ion mass in units of proton mass
    Mi = Mp * AA;

    OPTION(options, emass,                 false);  // including electron inertial, electron mass
    OPTION(options, emass_inv,               1.0);  // inverse of electron mass

    OPTION(options, diamag,                false);  // Diamagnetic effects?
    OPTION(options, energy_flux,           false);  // energy flux
    OPTION(options, energy_exch,           false);  // energy exchange
    OPTION(options, diamag_phi0,           false);  // Include equilibrium phi0
    OPTION(options, dia_fact,                1.0);  // Scale diamagnetic effects by this factor
    OPTION(options, diamag_er,             false);  // switch from phi0 to Er0

    OPTION(options, noshear,               false);

    OPTION(options, relax_j_vac,           false);  // Relax vacuum current to zero
    OPTION(options, relax_j_tconst,          0.1);

    // Toroidal filtering
    OPTION(options, filter_z,              false);  // Filter a single n
    OPTION(options, filter_z_mode,             1);
    OPTION(options, filter_z_nonlinear,    false);  // Filter a single n and zonal
    OPTION(options, low_pass_z,               -1);  // Low-pass filter
    OPTION(options, zonal_flow,               -1);  // zonal flow filter
    OPTION(options, zonal_field,              -1);  // zonal field filter
    OPTION(options, zonal_bkgd,               -1);  // zonal background P filter

    OPTION(options, filter_nl,                -1);  // zonal background P filter

    OPTION(options, limit_jacobi,          false);  // limit the value of jacobi at x-point region
    OPTION(options, bpxy_constraint,        0.04);
    OPTION(options, hthe_constraint,        0.04);
    OPTION(options, const_bp,                -1.);
    OPTION(options, PF_limit,              false);  // filter the instability in PF
    OPTION(options, PF_limit_range,          0.1);  // range of filter in PF
    OPTION(options, PF_sink,                 -1.);  // the coefficents of PF sink
    OPTION(options, PFs_width,               0.2);  // the percentage of radial grid points for sink profile radial width in PF
    OPTION(options, PFs_length,              0.4);  // the percentage of radial grid points for sink profile radial domain in PF

    // Radial smoothing
    OPTION(options, smooth_j_x,            false);  // Smooth Jpar in x
    OPTION(options, mask_j_x,              false);  // mask Jpar in x at boundary with tanh function
    OPTION(options, mask_flag_j,               1);  // mask flag, 0: mask on left boundary; 1: right; others: both on left and right
    OPTION(options, mask_phi_x,            false);  // mask phi in x at boundary with tanh function
    OPTION(options, mask_flag_phi,             0);  // mask flag, 0: mask on left boundary; 1: right; others: both on left and right
    OPTION(options, mask_length,             0.1);  // the center of tanh function
    OPTION(options, mask_width,             0.03);  // the width of tanh function

    // Jpar boundary region
    OPTION(options, jpar_bndry_width,         -1);

    // Parallel differencing
    OPTION(options, parallel_lr_diff,      false);
    OPTION(options, parallel_lagrange,     false);  // Use a (semi-) Lagrangian method for Grad_parP
    OPTION(options, parallel_project,      false);

    // Vacuum region control
    OPTION(options, vacuum_pressure,        0.02);  // Fraction of peak pressure
    OPTION(options, vacuum_trans,          0.005);  // Transition width in pressure

    // Resistivity and hyper-resistivity options
    OPTION(options, vac_lund,                0.0);  // Lundquist number in vacuum region
    OPTION(options, core_lund,               0.0);  // Lundquist number in core region
    OPTION(options, hyperresist,            -1.0);
    OPTION(options, ehyperviscos,           -1.0);
    OPTION(options, spitzer_resist,        false);  // Use Spitzer resistivity

    // Inner boundary damping
    OPTION(options, damp_width,                0);
    OPTION(options, damp_t_const,            0.1);

    // Viscosity and hyper-viscosity
    OPTION(options, viscos_par,             -1.0);  // Parallel viscosity
    OPTION(options, viscos_perp,            -1.0);  // Perpendicular viscosity
    OPTION(options, hyperviscos,            -1.0);  // Radial hyperviscosity

    OPTION(options, diffusion_par,          -1.0);  // parallel thermal conductivity is multiplied by this coefficient (<0 off)
    OPTION(options, kappa_par_i_const,      -1.0);  // add a constant heat conductivity. In 10^20m^-3m^2/s
    OPTION(options, kappa_par_e_const,      -1.0);
    OPTION(options, diffusion_perp,         -1.0);  // Perpendicular temperature diffusion
    OPTION(options, diff_par_flutter,      false);  // add magnetic flutter terms
    OPTION(options, full_sbc,              false);
    OPTION(options, fluxlimit,             false);
    OPTION(options, q_alpha,                 1.0);  // flux-limiting coefficient, typical value is [0.03, 3]
    OPTION(options, Landau,                false);  // Using Landau Fluid closure instead of thermal conductivity
    OPTION(options, Landau_coeff,            1.0);  // Coefficent for Landau Damping
    OPTION(options, Landau_coll,           false);  // Using collisional Landau Fluid closure instead of thermal conductivity
    OPTION(options, nLorentzian,               7);  // Number of Lorentzians
    OPTION(options, hyperdiff_par_n4,       -1.0);  // M: 4th Parallel density diffusion
    OPTION(options, hyperdiff_par_ti4,      -1.0);  // M: 4th Parallel ion temperature diffusion
    OPTION(options, hyperdiff_par_te4,      -1.0);  // M: 4th Parallel electron temperature diffusion
    OPTION(options, hyperdiff_par_v4,       -1.0);  // M: 4th Parallel ion parallel velocity diffusion
    OPTION(options, hyperdiff_par_apar4,    -1.0);
    OPTION(options, hyperdiff_par_u4,       -1.0);  // xqx: parallel hyper-viscous diffusion for vorticity
    OPTION(options, hyperdiff_perp_n4,      -1.0);  // M: 4th Perpendicular density diffusion
    OPTION(options, hyperdiff_perp_ti4,     -1.0);  // M: 4th Perpendicular ion temperature diffusion
    OPTION(options, hyperdiff_perp_te4,     -1.0);  // M: 4th Perpendicular electron temperature diffusion
    OPTION(options, hyperdiff_perp_v4,      -1.0);  // M: 4th Perpendicular ion parallel velocity diffusion
    OPTION(options, hyperdiff_perp_apar4,   -1.0);
    OPTION(options, hyperdiff_perp_u4,      -1.0);
    OPTION(options, neoclassic_i,          false);  // switch for ion neoclassical transport
    OPTION(options, neoclassic_e,          false);  // switch for electron neoclassical transport
    OPTION(options, neo_resist,            false);  // Include neoclassical effect in Spitzer resistivity
    OPTION(options, ft_simple,              true);  // use a simpler formula to calculate ft (the fraction of trapped electron) instead of Sauter's
    OPTION(options, major_radius,           0.68);  // R, in meter
    OPTION(options, minor_radius,           0.21);  // a, in meter

    OPTION(options, output_Teterms,        false);
    OPTION(options, output_Titerms,        false);
    OPTION(options, output_Tevegradte,     false);
    OPTION(options, output_transfer,       false);
    OPTION(options, output_ohm,            false);
    OPTION(options, output_flux_par,       false);
    OPTION(options, output_vradial,        false);
    OPTION(options, output_Tevegradte,     false);
    OPTION(options, output_qparcompare,    false);
    OPTION(options, output_P,                  true);
    OPTION(options, output_Vepar,              true);
    OPTION(options, output_SBC,                true);
    OPTION(options, output_eta,                true);
    OPTION(options, output_kappa_par,          false);

    // heating factor in pressure
    OPTION(options, heating_P,              -1.0);  // heating power in pressure
    OPTION(options, hp_width,                0.1);  // the percentage of radial grid points for heating profile radial width in pressure
    OPTION(options, hp_length,              0.04);  // the percentage of radial grid points for heating profile radial domain in pressure

    // sink factor in Vipar
    OPTION(options, sink_vp,                -1.0);  // sink in Vipar
    OPTION(options, sp_width,               0.05);  // the percentage of radial grid points for sink profile radial width in Vipar
    OPTION(options, sp_length,               0.1);  // the percentage of radial grid points for sink profile radial domain in Vipar

    // left edge sink factor in vorticity
    OPTION(options, sink_Ul,                -1.0);  // left edge sink in vorticity
    OPTION(options, su_widthl,              0.06);  // the percentage of left edge radial grid points for sink profile radial width in vorticity
    OPTION(options, su_lengthl,             0.15);  // the percentage of left edge radial grid points for sink profile radial domain in vorticity

    // right edge sink factor in vorticity
    OPTION(options, sink_Ur,                -1.0);  // right edge sink in vorticity
    OPTION(options, su_widthr,              0.06);  // the percentage of right edge radial grid points for sink profile radial width in vorticity
    OPTION(options, su_lengthr,             0.15);  // the percentage of right edge radial grid points for sink profile radial domain in vorticity

    // left edge sink factor in Te
    OPTION(options, sink_Tel,                   -1.0);   // left edge sink in Te
    OPTION(options, ste_widthl,                  0.06);  // the percentage of left edge radial grid points for sink profile radial width in Te
    OPTION(options, ste_lengthl,                 0.15);  // the percentage of left edge radial grid points for sink profile radial domain in Te
    // right edge sink factor in Te
    OPTION(options, sink_Ter,               -1.0);  // right edge sink in Te
    OPTION(options, ste_widthr,             0.06);  // the percentage of right edge radial grid points for sink profile radial width in Te
    OPTION(options, ste_lengthr,            0.15);  // the percentage of right edge radial grid points for sink profile radial domain in Te

    // right edge sink factor in Psi/Apar
    OPTION(options, sink_Psir,              -1.0);
    OPTION(options, spsi_widthr,                0.06);
    OPTION(options, spsi_lengthr,               0.15);

    // Compressional terms
    OPTION(options, phi_curv,               true);
    options->get("gamma", g, 5.0 / 3.0);

    // Field inversion flags
    OPTION(options, phi_flags,                 0);
    OPTION(options, apar_flags,                0);

    if (diffusion_par < 0. && output_flux_par) {
        output_flux_par = false;
        output.write("No parallel thermal conduction. Set 'output_flux_par' to be false.\n");
        if (diff_par_flutter) {
            diff_par_flutter = false;
            output.write("No parallel thermal conduction. Set 'diff_par_flutter' to be false.\n");
        }
    }

    if (!nonlinear) {
        if (output_transfer) {
            output_transfer = false;
            output.write("Linear simulation! Set output_transfer to false.\n");
        }
        if (output_ohm) {
            output_ohm = false;
            output.write("Linear simulation! Set output_ohm to false.\n");
        }
        if (output_flux_par) {
            output_flux_par = false;
            output.write("Linear simulation! Set output_flux_par to false.\n");
        }
        if (output_vradial) {
            output_vradial = false;
            output.write("Linear simulation! Set output_flux_par to false.\n");
        }
    }

    if (!include_curvature)
        b0xcv = 0.0;

    if (!include_jpar0)
        J0 = 0.0;

    if (noshear) {
        if (include_curvature)
            b0xcv.z += I * b0xcv.x;
        mesh->ShiftXderivs = false;
        I = 0.0;
    }

    //////////////////////////////////////////////////////////////
    // SHIFTED RADIAL COORDINATES

    if (mesh->ShiftXderivs) {
        if (mesh->IncIntShear) {
            // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
            mesh->IntShiftTorsion = I;

        } else {
            // Dimits style, using local coordinate system
            if (include_curvature)
                b0xcv.z += I * b0xcv.x;
            I = 0.0; // I disappears from metric
        }
    }

    //////////////////////////////////////////////////////////////
    // NORMALISE QUANTITIES

    if (mesh->get(Bbar, "bmag")) // Typical magnetic field
        Bbar = 1.0;

    if (mesh->get(Lbar, "rmag")) // Typical length scale
        Lbar = 1.0;

    if (mesh->get(Tibar, "Ti_x")) // Typical ion temperature scale
        Tibar = 1.0;

    if (mesh->get(Tebar, "Te_x")) // Typical electron temperature scale
        Tebar = 1.0;

    if (mesh->get(Nbar, "Nixexp")) // Typical ion density scale
        Nbar = 1.0;
    Nbar *= density_unit / density;

    Tau_ie = Tibar / Tebar;

    Va = sqrt(Bbar * Bbar / (MU0 * Mi * Nbar * density));

    Tbar = Lbar / Va;

    output.write("Normalisations:\n");
    output.write("\tBbar = %e T   Lbar = %e m\n", Bbar, Lbar);
    output.write("\tVa = %e m/s   Tbar = %e s\n", Va, Tbar);
    output.write("\tNbar = %e * %e m^-3\n", Nbar, density);
    output.write("\tPibar = %e Pa   Pebar = %e Pa\n", ee * Tibar * Nbar * density, ee * Tebar * Nbar * density);
    output.write("\tTibar = %e eV   Tebar = %e eV    Ti/Te = %e\n", Tibar, Tebar, Tau_ie);
    output.write("\tetabar = %e [Ohm m]\n", MU0 * Va * Lbar);

    if (emass) {
        delta_e = 5.31e5 / sqrt(Nbar * density / 1e6) / (Lbar * 100.0) * emass_inv;
        delta_e_inv = 1.e0 / delta_e / delta_e;
        gyroAlv = 1.602e-19 * Bbar * Tbar / Mi;
        output.write("                delta_e = %e    wci*T_A = %e\n", delta_e, gyroAlv);
    }

    if (thermal_force || eHall) {
        Psipara1 = KB * Tebar * eV_K / ee / Bbar / Lbar / Va;
        output.write("                Psipara1 = %e   AA = %e\n", Psipara1, AA);
    }

    Upara0 = KB * Tibar * eV_K / (Zi * ee * Bbar * Va * Lbar);
    Upara1 = KB * Tebar * eV_K / (Mi * Va * Va);
    output.write("vorticity constant: Upara0 = %e     Upara1 = %e\n", Upara0, Upara1);

    if (gyroviscous) {
        Upara2 = KB * Tibar * eV_K / (Zi * ee * Bbar * Lbar * Va);
        Upara3 = 1.0;
        output.write("Upara2 = %e     Upara3 = %e\n", Upara2, Upara3);
    }

    if ((diamag && continuity) || energy_flux) {
        Nipara1 = KB * Tibar * eV_K / (Zi * ee * Bbar * Lbar * Va);
        Tipara2 = Nipara1;
        Tipara3 = Mi / Zi / ee / Bbar / Tbar;
        Tepara2 = KB * Tebar * eV_K / (ee * Bbar * Lbar * Va);
        Tepara3 = Bbar / (ee * MU0 * Nbar * density * Lbar * Va);
        output.write("Nipara1 = %e     Tipara2 = %e\n", Nipara1, Tipara2);
        output.write("Tepara2 = %e     Tepara3 = %e\n", Tepara2, Tepara3);
    }

    if (energy_exch) {
        Tepara4 = Bbar * Bbar / (MU0 * KB * Nbar * density * Tebar * eV_K);
        output.write("energy exchange constant:   Tepara4 = %e\n", Tepara4);
    }

    if (compress0) {
        output.write("Including compression (Vipar) effects\n");
        if (include_vipar)
            output.write("Including all compresssion (Vipar) effects\n");
        Vipara = MU0 * KB * Nbar * density * Tebar * eV_K / (Bbar * Bbar);
        Vepara = Bbar / (MU0 * ee * Nbar * density * Lbar * Va);
        output.write("Normalized constant for Vipar :   Vipara = %e\n", Vipara);
        output.write("Normalized constant for Vepar :   Vepara = %e\n", Vepara);
    } else
        include_vipar = false;

    if (diffusion_par > 0.0 || diffusion_perp > 0.0) {
        Tipara1 = 2.0 / 3.0 / (Lbar * Va);
        Tepara1 = Tipara1;
    }

    if (hyperresist > 0.0) {
        output.write("    Hyper-resistivity coefficient: %e\n", hyperresist);
        dump.add(hyper_eta_x, "hyper_eta_x", 1);
        dump.add(hyper_eta_z, "hyper_eta_z", 1);
    }

    if (ehyperviscos > 0.0) {
        output.write("    electron Hyper-viscosity coefficient: %e\n", ehyperviscos);
    }

    if (hyperviscos > 0.0) {
        output.write("    Hyper-viscosity coefficient: %e\n", hyperviscos);
        dump.add(hyper_mu_x, "hyper_mu_x", 1);
    }

    // dump.add(sink_PFtmp, "sink_PFtmp",1);

    if (diffusion_par > 0.0) {
        output.write("    diffusion_par: %e\n", diffusion_par);
        // dump.add(diffusion_par, "diffusion_par", 0);
    }

    if (diffusion_perp > 0.0) {
        output.write("    diffusion_perp: %e\n", diffusion_perp);
        dump.add(diffusion_perp, "diffusion_perp", 0);
    }

    // M: 4th order diffusion of p
    if (hyperdiff_par_n4 > 0.0) {
        output.write(" parallel hyperdiffusion_n4: %e\n", hyperdiff_par_n4);
        dump.add(hyperdiff_par_n4, "hyperdiff_par_n4", 0);
    }

    // M: 4th order diffusion of Ti
    if (hyperdiff_par_ti4 > 0.0) {
        output.write(" parallel hyperdiffusion_ti4: %e\n", hyperdiff_par_ti4);
        dump.add(hyperdiff_par_ti4, "hyperdiff_par_ti4", 0);
    }

    // M: 4th order diffusion of Te
    if (hyperdiff_par_te4 > 0.0) {
        output.write(" parallel hyperdiffusion_te4: %e\n", hyperdiff_par_te4);
        dump.add(hyperdiff_par_te4, "hyperdiff_par_te4", 0);
    }

    // M: 4th order diffusion of Vipar
    if (hyperdiff_par_v4 > 0.0) {
        output.write(" parallel hyperdiffusion_v4: %e\n", hyperdiff_par_v4);
        dump.add(hyperdiff_par_v4, "hyperdiff_par_v4", 0);
    }

    if (hyperdiff_par_apar4 > 0.0) {
        output.write(" parallel hyperdiffusion_apar4: %e\n", hyperdiff_par_apar4);
        dump.add(hyperdiff_par_apar4, "hyperdiff_par_apar4", 0);
    }

    // xqx: parallel hyper-viscous diffusion for vorticity
    if (hyperdiff_par_u4 > 0.0) {
        output.write(" parallel hyperdiffusion_u4: %e\n", hyperdiff_par_u4);
        dump.add(hyperdiff_par_u4, "hyperdiff_par_u4", 0);
    }

    if (hyperdiff_perp_n4 > 0.0) {
        output.write(" perpendicular hyperdiffusion_n4: %e\n", hyperdiff_perp_n4);
        dump.add(hyperdiff_perp_n4, "hyperdiff_perp_n4", 0);
    }

    if (hyperdiff_perp_ti4 > 0.0) {
        output.write(" perpendicular hyperdiffusion_ti4: %e\n", hyperdiff_perp_ti4);
        dump.add(hyperdiff_perp_ti4, "hyperdiff_perp_ti4", 0);
    }

    if (hyperdiff_perp_te4 > 0.0) {
        output.write(" perpendicular hyperdiffusion_te4: %e\n", hyperdiff_perp_te4);
        dump.add(hyperdiff_perp_te4, "hyperdiff_perp_te4", 0);
    }

    if (hyperdiff_perp_v4 > 0.0) {
        output.write(" perpendicular hyperdiffusion_v4: %e\n", hyperdiff_perp_v4);
        dump.add(hyperdiff_perp_v4, "hyperdiff_perp_v4", 0);
    }

    if (hyperdiff_perp_apar4 > 0.0) {
        output.write(" perpendicular hyperdiffusion_apar4: %e\n", hyperdiff_perp_apar4);
        dump.add(hyperdiff_perp_apar4, "hyperdiff_perp_apar4", 0);
    }

    if (hyperdiff_perp_u4 > 0.0) {
        output.write(" perpendicular hyperdiffusion_u4: %e\n", hyperdiff_perp_u4);
        dump.add(hyperdiff_perp_u4, "hyperdiff_perp_u4", 0);
    }

    if (sink_vp > 0.0) {
        output.write("    sink_vp(rate): %e\n", sink_vp);
        dump.add(sink_vp, "sink_vp", 1);

        output.write("    sp_width(%): %e\n", sp_width);
        dump.add(sp_width, "sp_width", 1);

        output.write("    sp_length(%): %e\n", sp_length);
        dump.add(sp_length, "sp_length", 1);
    }

    if (limit_jacobi) {
        Bpxy = field_larger(Bpxy, bpxy_constraint);
        // Bpxy = smooth_xy(Bpxy);
        dump.add(Bpxy, "Bpxy", 0);
        if (const_bp > 0.)
            Bpxy = const_bp;
        if (hthe_constraint > 0.) {
            hthe = field_larger(hthe, hthe_constraint);
            dump.add(hthe, "hthe", 0);
        }
    }

    J0 = MU0 * Lbar * J0 / Bbar;
    P0 = P0 / (KB * (Tebar)*eV_K * Nbar * density);

    b0xcv.x /= Bbar;
    b0xcv.y *= Lbar * Lbar;
    b0xcv.z *= Lbar * Lbar;

    Rxy /= Lbar;
    Bpxy /= Bbar;
    Btxy /= Bbar;
    B0 /= Bbar;
    hthe /= Lbar;
    major_radius /= Lbar;
    minor_radius /= Lbar;
    mesh->dx /= Lbar * Lbar * Bbar;
    I *= Lbar * Lbar * Bbar;

    if ((!T0_fake_prof) && n0_fake_prof) {
        if (n0_p0_0p3) {
            output.write("N0 ~ P0^0.3 used!\n");
            // n0 = n0_height*(P0/P00)^0.3*density
            // this relation is used in the bootstrap current calculation
            // when generating the cbm18_dens_ne* serial grids
            // P00: P0 at magnetic axis, for cbm18_dens6, it's 23072.3
            // densn: P00 = 23072.3*n/6.
            BoutReal P00; // P0 at magnetic axis
            if (options->isSet("P00"))
                OPTION(options, P00,                  23072.3);
            else {
                output.write("n0=n0_height*(P0/P00)^0.3*density is used, P00 is required!\n");
                return 1;
            }
            P00 = P00 / (ee * Tebar * Nbar * density);
            N0 = n0_height * ((P0 / P00) ^ 0.3);
        } else
            N0 = N0tanh(n0_height, n0_ave, n0_width, n0_center, n0_bottom_x);

        Ti0 = P0 / N0 / (1.0 + Zi);
        Te0 = Ti0;
        N_imp0 = 0;
        Ne0 = Zi * N0;
    } else if (T0_fake_prof) {
        Ti0 = Tconst;
        Te0 = Ti0;
        N0 = P0 / (Ti0 * Tau_ie + Te0);
        N_imp0 = 0;
        Ne0 = Zi * N0;
    } else {
        if (mesh->get(N0, "Niexp")) { // N_i0
            output.write("Error: Cannot read Ni0 from grid\n");
            return 1;
        }

        if (mesh->get(Ti0, "Tiexp")) { // T_i0
            output.write("Error: Cannot read Ti0 from grid\n");
            return 1;
        }

        if (mesh->get(Te0, "Teexp")) { // T_e0
            output.write("Error: Cannot read Te0 from grid\n");
            return 1;
        }

        N0 /= Nbar;
        Ti0 /= Tibar;
        Te0 /= Tebar;

        if (impurity_prof) {
            if (mesh->get(Ne0, "Neexp")) { // N_e0
                output.write("Error: Cannot read Ne0 from grid\n");
                return 1;
            }
            Ne0 /= Nbar;
            if (load_impurity) {
                if (mesh->get(N_imp0, "N_imp")) {
                    output.write("Error: Cannot read N_imp from grid\n");
                    return 1;
                }
                N_imp0 /= Nbar;
            } else
                N_imp0 = (Ne0 - Zi * N0) / Z_imp;
            if (min(N_imp0, true) <= 0.) {

                output.write("Error: Impurity density has negative value.\n");
                if (Nimp_lowlimit) {
                    BoutReal min_Nimp = max(N_imp0, true) * 0.01;
                    output.write("Use Nimp_lowlimit of %d to modify the smallest value.\n", min_Nimp);
                    for (jx = 0; jx < mesh->ngx; jx++)
                        for (jy = 0; jy < mesh->ngy; jy++)
                            if (N_imp0[jx][jy] <= min_Nimp)
                                N_imp0[jx][jy] = min_Nimp;
                } else
                    return 1;
            }
            T_imp0 = Ti0;
            P_imp0 = N_imp0 * T_imp0;
            // P_imp0 = P0 - N0*Ti0*Tau_ie - Ne0*Te0;
            // T_imp0 = P_imp0/N_imp0;
            if (min(T_imp0, true) <= 0.) {
                output.write("Error: Impurity temperature has negative value.\n");
                if (Nimp_lowlimit) {
                    BoutReal min_Timp = max(T_imp0, true) * 0.01;
                    output.write("Use Timp_lowlimit of %d to modify the smallest value.\n", min_Timp);
                    for (jx = 0; jx < mesh->ngx; jx++)
                        for (jy = 0; jy < mesh->ngy; jy++)
                            if (T_imp0[jx][jy] <= min_Timp)
                                T_imp0[jx][jy] = min_Timp;
                    P_imp0 = N_imp0 * T_imp0;
                } else
                    return 1;
            }
        } else if (quasi_neutral_Ni) {
            N0 = P0 / (Tau_ie * Ti0 + Te0);
            Ne0 = Zi * N0;
            N_imp0 = 0;
        } else {
            Ne0 = Zi * N0; // quasi-neutral condition
            N_imp0 = 0;
        }
    }

    J0 *= J0_factor;
    P0 *= P0_factor;


    Ne0.applyBoundary("neumann");

    Pi0 = N0 * Ti0;
    Pe0 = Ne0 * Te0;
    Pi0.applyBoundary("neumann");
    Pe0.applyBoundary("neumann");


    jpar1.setBoundary("J");
    u_tmp1.setBoundary("U");
    ti_tmp2.setBoundary("Ti");
    te_tmp2.setBoundary("Te");
    phi_tmp.setBoundary("phi");

    nu_e.setLocation(CELL_YLOW);
    nu_e.setBoundary("kappa");
    if (spitzer_resist) {
        eta.setLocation(CELL_YLOW);
        eta.setBoundary("kappa");
    }
    if (diffusion_par > 0.0 || diffusion_perp > 0.0 || neo_resist) {
        nu_i.setLocation(CELL_YLOW);
        nu_i.setBoundary("kappa");
        vth_i.setLocation(CELL_YLOW);
        vth_e.setLocation(CELL_YLOW);
        vth_i.setBoundary("kappa");
        vth_e.setBoundary("kappa");
        kappa_par_i.setLocation(CELL_YLOW);
        kappa_par_e.setLocation(CELL_YLOW);
        kappa_par_i.setBoundary("kappa");
        kappa_par_e.setBoundary("kappa");
        kappa_perp_i.setLocation(CELL_YLOW);
        kappa_perp_e.setLocation(CELL_YLOW);
        kappa_perp_i.setBoundary("kappa");
        kappa_perp_e.setBoundary("kappa");
        if (diff_par_flutter) {
            bracket1i.setLocation(CELL_CENTRE);
            bracket1i.setBoundary("Ti");
            gradpar_ti.setLocation(CELL_CENTRE);
            gradpar_ti.setBoundary("Ti");
            bracket1e.setLocation(CELL_CENTRE);
            bracket1e.setBoundary("Te");
            gradpar_te.setLocation(CELL_CENTRE);
            gradpar_te.setBoundary("Te");
        }
        if (Landau) {
            q_par_e.setLocation(CELL_CENTRE);
            q_par_i.setLocation(CELL_CENTRE);
            q_par_e.setBoundary("q_par_e");
            q_par_i.setBoundary("q_par_i");
            q_par_e = 0.;
            q_par_i = 0.;
            dump.add(q_par_e, "q_par_e", 1);
            dump.add(q_par_i, "q_par_i", 1);
            if (output_qparcompare) {
                q_par_fl.setLocation(CELL_CENTRE);
                q_par_fl.setBoundary("q_par_e");
                dump.add(q_par_fl, "q_par_fl", 1);
                q_par_landau.setLocation(CELL_CENTRE);
                q_par_landau.setBoundary("q_par_e");
                dump.add(q_par_landau, "q_par_landau", 1);
            }
            if (Landau_coll) {
                kappa_i.setLocation(CELL_YLOW);
                kappa_i.setBoundary("kappa");
                kappa_e.setLocation(CELL_YLOW);
                kappa_e.setBoundary("kappa");
            }
            SBC_value_i.setLocation(CELL_YLOW);
            SBC_value_i.setBoundary("kappa");
            SBC_value_e.setLocation(CELL_YLOW);
            SBC_value_e.setBoundary("kappa");
        }
    }

    if (neoclassic_i || neoclassic_e) {
        xii_neo.setLocation(CELL_YLOW);
        xii_neo.setBoundary("kappa");
        xie_neo.setLocation(CELL_YLOW);
        xie_neo.setBoundary("kappa");
        // Dri_neo.setLocation(CELL_YLOW);
        // Dri_neo.setBoundary("kappa");
        rho_i.setLocation(CELL_YLOW);
        rho_i.setBoundary("kappa");
        rho_e.setLocation(CELL_YLOW);
        rho_e.setBoundary("kappa");
        tmpddx2.setLocation(CELL_YLOW);
        tmpddx2.setBoundary("kappa");
        partf_neo_i.setLocation(CELL_CENTRE);
        partf_neo_i.setBoundary("Ni");
        heatf_neo_i.setLocation(CELL_CENTRE);
        heatf_neo_i.setBoundary("Ti");
        heatf_neo_e.setLocation(CELL_CENTRE);
        heatf_neo_e.setBoundary("Te");
    }

    if (parallel_viscous && compress0) {
        eta_i0.setLocation(CELL_CENTRE);
        eta_i0.setBoundary("Ti");
        pi_ci.setLocation(CELL_CENTRE);
        pi_ci.setBoundary("Ti");

        // dump.add(eta_i0, "eta_i0", 1);
        // dump.add(pi_ci, "pi_ci", 1);
    }

    if (gyroviscous) {
        Dperp2Phi0.setLocation(CELL_CENTRE);
        Dperp2Phi0.setBoundary("phi");
        Dperp2Phi.setLocation(CELL_CENTRE);
        Dperp2Phi.setBoundary("phi");
        GradPhi02.setLocation(CELL_CENTRE);
        GradPhi02.setBoundary("phi");
        // GradparPhi02.setLocation(CELL_CENTRE);
        // GradparPhi02.setBoundary("phi");
        GradcPhi.setLocation(CELL_CENTRE);
        GradcPhi.setBoundary("phi");
        // GradcparPhi.setLocation(CELL_CENTRE);
        // GradcparPhi.setBoundary("phi");
        Dperp2Pi0.setLocation(CELL_CENTRE);
        Dperp2Pi0.setBoundary("P");
        Dperp2Pi.setLocation(CELL_CENTRE);
        Dperp2Pi.setBoundary("P");
        bracketPhi0P.setLocation(CELL_CENTRE);
        bracketPhi0P.setBoundary("P");
        bracketPhiP0.setLocation(CELL_CENTRE);
        bracketPhiP0.setBoundary("P");
        if (impurity_prof && impurity_gyro) {
            Upara_imp = (A_imp * N_imp0) / (AA * N0);
            output.write("Max Upara_imp = %e, Min Upara_imp = %e\n", max(Upara_imp, true), min(Upara_imp, true));

            Dperp2Pimp0.setLocation(CELL_CENTRE);
            Dperp2Pimp0.setBoundary("P");
            bracketPhiPimp0.setLocation(CELL_CENTRE);
            bracketPhiPimp0.setBoundary("P");
        }
        if (nonlinear) {
            GradPhi2.setLocation(CELL_CENTRE);
            GradPhi2.setBoundary("phi");
            // GradparPhi2.setLocation(CELL_CENTRE);
            // GradparPhi2.setBoundary("phi");
            bracketPhiP.setLocation(CELL_CENTRE);
            bracketPhiP.setBoundary("P");
        }
    }

    if (output_transfer) {
        T_R.setLocation(CELL_YLOW);
        T_R.setBoundary("phi");
        T_M.setLocation(CELL_YLOW);
        T_M.setBoundary("phi");
        T_ID.setLocation(CELL_YLOW);
        T_ID.setBoundary("phi");
        T_C.setLocation(CELL_YLOW);
        T_C.setBoundary("P");
        T_G.setLocation(CELL_YLOW);
        T_G.setBoundary("P");

        dump.add(T_R, "T_R", 1);
        dump.add(T_M, "T_M", 1);
        dump.add(T_ID, "T_ID", 1);
        dump.add(T_C, "T_C", 1);
        dump.add(T_G, "T_G", 1);
    }

    if (output_ohm) {
        ohm_phi.setLocation(CELL_YLOW);
        ohm_phi.setBoundary("Psi");
        ohm_hall.setLocation(CELL_YLOW);
        ohm_hall.setBoundary("Psi");
        ohm_thermal.setLocation(CELL_YLOW);
        ohm_thermal.setBoundary("Psi");

        dump.add(ohm_phi, "ohm_phi", 1);
        dump.add(ohm_hall, "ohm_hall", 1);
        dump.add(ohm_thermal, "ohm_thermal", 1);
    }

    if (output_flux_par && diffusion_par > 0.) {
        // gamma_par_i.setLocation(CELL_YLOW);
        // gamma_par_i.setBoundary("Ni");
        heatf_par_i.setLocation(CELL_CENTRE);
        heatf_par_i.setBoundary("Ti");
        heatf_par_e.setLocation(CELL_CENTRE);
        heatf_par_e.setBoundary("Te");

        // dump.add(gamma_par_i, "gamma_i", 1);
        dump.add(heatf_par_i, "heatflux_par_i", 1);
        dump.add(heatf_par_e, "heatflux_par_e", 1);
        if (diff_par_flutter) {
            heatf_par_flutter_i.setLocation(CELL_CENTRE);
            heatf_par_flutter_i.setBoundary("Ti");
            heatf_par_flutter_e.setLocation(CELL_CENTRE);
            heatf_par_flutter_e.setBoundary("Te");

            // dump.add(gamma_par_i, "gamma_i", 1);
            dump.add(heatf_par_flutter_i, "heatflux_par_flutter_i", 1);
            dump.add(heatf_par_flutter_e, "heatflux_par_flutter_e", 1);
        }
    }

    if (output_vradial) {
        Vexb.covariant = false;
        Vexb.setLocation(CELL_YLOW);
        Vexb.setBoundary("Vipar");
        Vbtilde.covariant = false;
        Vbtilde.setLocation(CELL_YLOW);
        Vbtilde.setBoundary("Vipar");
        // Vbti_par.setLocation(CELL_YLOW);
        // Vbti_par.setBoundary("Vipar");
        // Vbte_par.setLocation(CELL_YLOW);
        // Vbte_par.setBoundary("Vipar");

        dump.add(Vexb, "Vexb", 1);
        dump.add(Vbtilde, "Vbtild", 1);
    }

    if (mask_j_x) {
        mask_jx1d = mask_x_1d(false, mask_flag_j, mask_width, mask_length);
        // dump.add(mask_jx1d, "mask_jx1d", 0);
    }

    if (mask_phi_x) {
        mask_px1d = mask_x_1d(false, mask_flag_phi, mask_width, mask_length);
        // dump.add(mask_px1d, "mask_px1d", 0);
    }

    BoutReal pnorm = max(P0, true); // Maximum over all processors

    vacuum_pressure *= pnorm; // Get pressure from fraction
    vacuum_trans *= pnorm;

    // Transitions from 0 in core to 1 in vacuum
    vac_mask = (1.0 - tanh((P0 - vacuum_pressure) / vacuum_trans)) / 2.0;

    LnLambda = 24.0 - log(sqrt(Zi * Nbar * density * Ne0 / 1.e6) / (Tebar * Te0));
    output.write("\tlog Lambda: %e -> %e \n", min(LnLambda), max(LnLambda));

    if (spitzer_resist || neo_resist) {
        // Use Spitzer resistivity
        output.write("");
        output.write("\tSpizter parameters");
        // output.write("\tTemperature: %e -> %e [eV]\n", min(Te), max(Te));
        FZ = (1. + 1.198 * Zeff + 0.222 * Zeff * Zeff) / (1. + 2.996 * Zeff + 0.753 * Zeff * Zeff);
        eta = FZ * 1.03e-4 * Zeff * LnLambda * ((Te0 * Tebar) ^ (-1.5)); // eta in Ohm-m. NOTE: ln(Lambda) = 20
        output.write("\tSpitzer resistivity: %e -> %e [Ohm m]\n", min(eta), max(eta));
        eta /= MU0 * Va * Lbar;
        // eta.applyBoundary();
        // mesh->communicate(eta);
        output.write("\t -> Lundquist %e -> %e\n", 1.0 / max(eta), 1.0 / min(eta));
        if (output_eta)
            dump.add(eta, "eta", 1);
        else
            dump.add(eta, "eta", 0);
    } else {
        // transition from 0 for large P0 to resistivity for small P0
        if (vac_lund > 0.0) {
            output.write("        Vacuum  Tau_R = %e s   eta = %e Ohm m\n", vac_lund * Tbar, MU0 * Lbar * Lbar / (vac_lund * Tbar));
            vac_resist = 1. / vac_lund;
        } else {
            output.write("        Vacuum  - Zero resistivity -\n");
            vac_resist = 0.0;
        }

        if (core_lund > 0.0) {
            output.write("        Core    Tau_R = %e s   eta = %e Ohm m\n", core_lund * Tbar, MU0 * Lbar * Lbar / (core_lund * Tbar));
            core_resist = 1. / core_lund;
        } else {
            output.write("        Core    - Zero resistivity -\n");
            core_resist = 0.0;
        }

        eta = core_resist + (vac_resist - core_resist) * vac_mask;
        dump.add(eta, "eta", 0);
    }

    if (diffusion_par > 0.0 || diffusion_perp > 0.0 || neoclassic_i || neoclassic_e || neo_resist) {
        if (q95_input > 0)
            q95 = q95_input; // use a constant for test
        else {
            if (local_q)
                q95 = abs(hthe * Btxy / (Bpxy));
            else {
                output.write("\tUsing q profile from grid.\n");
                if (mesh->get(q95, "q")) {
                    output.write("Cannot get q profile from grid!\nPlease run addqprofile.pro first\n");
                    return 1;
                }
            }
        }

        output.write("\tlocal max q: %e\n", max(q95));
        output.write("\tlocal min q: %e\n", min(q95));
    }

//  LnLambda = 24.0 - log(pow(Zi * Nbar * density / 1.e6, 0.5) * pow(Tebar, -1.0)); // xia: ln Lambda
//  output.write("\tlog Lambda: %e\n", LnLambda);

    nu_e = 2.91e-6 * LnLambda * ((Ne0)*Nbar * density / 1.e6) * (((Te0)*Tebar) ^ (-1.5)); // nu_e in 1/S.
    output.write("\telectron collision rate: %e -> %e [1/s]\n", min(nu_e), max(nu_e));
    // nu_e.applyBoundary();
    // mesh->communicate(nu_e);

    if (diffusion_par > 0.0 || diffusion_perp > 0.0 || parallel_viscous || neoclassic_i || neoclassic_e || neo_resist) {

        output.write("\tion thermal noramlized constant: Tipara1 = %e\n", Tipara1);
        output.write("\telectron normalized thermal constant: Tepara1 = %e\n", Tepara1);
        // xqx addition, begin
        // Use Spitzer thermal conductivities
        nu_i = 4.80e-8 * (Zi * Zi * Zi * Zi / sqrt(AA)) * LnLambda * ((N0)*Nbar * density / 1.e6) * (((Ti0)*Tibar) ^ (-1.5)); // nu_i in 1/S.
        // output.write("\tCoulomb Logarithm: %e\n", max(LnLambda));
        output.write("\tion collision rate: %e -> %e [1/s]\n", min(nu_i), max(nu_i));

        // nu_i.applyBoundary();
        // mesh->communicate(nu_i);

        vth_i = 9.79e3 * sqrt((Ti0)*Tibar / AA); // vth_i in m/S.
        output.write("\tion thermal velocity: %e -> %e [m/s]\n", min(vth_i), max(vth_i));
        // vth_i.applyBoundary();
        // mesh->communicate(vth_i);
        vth_e = 4.19e5 * sqrt((Te0)*Tebar); // vth_e in m/S.
        output.write("\telectron thermal velocity: %e -> %e [m/s]\n", min(vth_e), max(vth_e));
        // vth_e.applyBoundary();
        // mesh->communicate(vth_e);
    }

    if (parallel_viscous && compress0) {
        //eta_i0 = 0.96 * Pi0 * Tau_ie * nu_i * Tbar;
        eta_i0 = 0.96 * Pi0 * Tau_ie / nu_i / Tbar;
        output.write("\tCoefficients of parallel viscocity: %e -> %e [kg/(m s)]\n", min(eta_i0), max(eta_i0));
    }

    if (diffusion_par > 0.0) {
//        kappa_par_i = 3.9 * vth_i * vth_i / nu_i; // * 1.e4;
//        kappa_par_e = 3.2 * vth_e * vth_e / nu_e; // * 1.e4;
//
//        output.write("\tion thermal conductivity: %e -> %e [m^2/s]\n", min(kappa_par_i), max(kappa_par_i));
//        output.write("\telectron thermal conductivity: %e -> %e [m^2/s]\n", min(kappa_par_e), max(kappa_par_e));
//
//        output.write("\tnormalized ion thermal conductivity: %e -> %e\n", min(kappa_par_i * Tipara1), max(kappa_par_i * Tipara1));
//        output.write("\tnormalized electron thermal conductivity: %e -> %e\n", min(kappa_par_e * Tepara1), max(kappa_par_e * Tepara1));
//
//        kappa_par_i_fl = q_alpha * vth_i * q95 * Lbar; // * 1.e2;
//        kappa_par_e_fl = q_alpha * vth_e * q95 * Lbar; // * 1.e2;
//
//        if (fluxlimit) {
//            kappa_par_i *= kappa_par_i_fl / (kappa_par_i + kappa_par_i_fl);
//            kappa_par_e *= kappa_par_e_fl / (kappa_par_e + kappa_par_e_fl);
//        }
//        kappa_par_i *= diffusion_par * Tipara1 * N0;
//        output.write("\tUsed normalized ion thermal conductivity: %e -> %e\n", min(kappa_par_i), max(kappa_par_i));
//        // kappa_par_i.applyBoundary();
//        // mesh->communicate(kappa_par_i);
//        kappa_par_e *= diffusion_par * Tepara1 * Ne0;
//        output.write("\tUsed normalized electron thermal conductivity: %e -> %e\n", min(kappa_par_e), max(kappa_par_e));
//        // kappa_par_e.applyBoundary();
//        // mesh->communicate(kappa_par_e);
//
//        dump.add(kappa_par_i, "kappa_par_i", 1);
//        dump.add(kappa_par_e, "kappa_par_e", 1);

        kappa_par_i_sp = 3.9 * vth_i * vth_i / nu_i;  // in SI units
        kappa_par_e_sp = 3.2 * vth_e * vth_e / nu_e;
        kappa_par_i_sp *= Tipara1 * N0;  // normalized 2/3*Spitzer-Harm
        kappa_par_e_sp *= Tepara1 * Ne0;
        mesh->communicate(kappa_par_i_sp);
        mesh->communicate(kappa_par_e_sp);
        if (output_kappa_par) {
            dump.add(kappa_par_i_sp, "kappa_par_i_sp", 1);
            dump.add(kappa_par_e_sp, "kappa_par_e_sp", 1);
        } else {
            dump.add(kappa_par_i_sp, "kappa_par_i_sp", 0);
            dump.add(kappa_par_e_sp, "kappa_par_e_sp", 0);
        }
        if (fluxlimit) {
            kappa_par_i_fl = q_alpha * vth_i * q95 * major_radius * Lbar;  // in SI units
            kappa_par_e_fl = q_alpha * vth_e * q95 * major_radius * Lbar;
//            kappa_par_i_fl = q_alpha * vth_i * q95 * Lbar;
//            kappa_par_e_fl = q_alpha * vth_e * q95 * Lbar;
            kappa_par_i_fl *= Tipara1 * N0;  // normalized 2/3*flux-limited
            kappa_par_e_fl *= Tepara1 * Ne0;
            mesh->communicate(kappa_par_i_fl);
            mesh->communicate(kappa_par_e_fl);
            if (output_kappa_par) {
                dump.add(kappa_par_i_fl, "kappa_par_i_fl", 1);
                dump.add(kappa_par_e_fl, "kappa_par_e_fl", 1);
            } else {
                dump.add(kappa_par_i_fl, "kappa_par_i_fl", 0);
                dump.add(kappa_par_e_fl, "kappa_par_e_fl", 0);
            }
            kappa_par_i = diffusion_par * (kappa_par_i_sp * kappa_par_i_fl);
            kappa_par_i /= (kappa_par_i_sp + kappa_par_i_fl);
            kappa_par_e = diffusion_par * (kappa_par_e_sp * kappa_par_e_fl);
            kappa_par_e /= (kappa_par_e_sp + kappa_par_e_fl);
        } else {
            kappa_par_i = diffusion_par * kappa_par_i_sp;
            kappa_par_e = diffusion_par * kappa_par_e_sp;
        }
        if (kappa_par_i_const > 0.0) {
            kappa_par_i_const *= Tipara1 * 1e20 / (Nbar * density);
            dump.add(kappa_par_i_const, "kappa_par_i_const", 0);
            kappa_par_i += kappa_par_i_const;
        }
        if (kappa_par_e_const > 0.0) {
            kappa_par_e_const *= Tepara1 * 1e20 / (Nbar * density);
            dump.add(kappa_par_e_const, "kappa_par_e_const", 0);
            kappa_par_e += kappa_par_e_const;
        }
        mesh->communicate(kappa_par_i);
        mesh->communicate(kappa_par_e);
        if (output_kappa_par) {
            dump.add(kappa_par_i, "kappa_par_i", 1);
            dump.add(kappa_par_e, "kappa_par_e", 1);
        } else {
            dump.add(kappa_par_i, "kappa_par_i", 0);
            dump.add(kappa_par_e, "kappa_par_e", 0);
        }

        kappa_par_i_lin = B0;
        kappa_par_e_lin = B0;
        for (jx = 0; jx < mesh->ngx; jx++)
            for (jy = 0; jy < mesh->ngy; jy++) {
                kappa_par_i_lin[jx][jy] = kappa_par_i[jx][jy][0];
                kappa_par_e_lin[jx][jy] = kappa_par_e[jx][jy][0];
            }

        if (Landau) {
            if (Landau_coll) {
                output.write("Collisional Landau damping closure used!\n");
                kappa_i = 0.5 * nu_i / vth_i;
                kappa_i *= Lbar; // normalization
                kappa_e = 0.5 * nu_e / vth_e;
                kappa_e *= Lbar; // normalization
            } else {
                output.write("Collisionless Landau damping closure used!\n");
                kappa_0 = 2. / (20. * (2 * PI * q95_input));
                output.write("\tkappa_0 = %e\n", kappa_0);
            }
        }
    } else {
        if (kappa_par_i_const > 0.0) {
            kappa_par_i_const *= Tipara1 * 1e20 / (Nbar * density);
            dump.add(kappa_par_i_const, "kappa_par_i_const", 0);
            kappa_par_i = kappa_par_i_const;
            if (output_kappa_par) {
                dump.add(kappa_par_i, "kappa_par_i", 1);
            } else {
                dump.add(kappa_par_i, "kappa_par_i", 0);
            }
        }
        if (kappa_par_e_const > 0.0) {
            kappa_par_e_const *= Tepara1 * 1e20 / (Nbar * density);
            dump.add(kappa_par_e_const, "kappa_par_e_const", 0);
            kappa_par_e = kappa_par_e_const;
            if (output_kappa_par) {
                dump.add(kappa_par_e, "kappa_par_e", 1);
            } else {
                dump.add(kappa_par_e, "kappa_par_e", 0);
            }
        }

    }

    if (diffusion_perp > 0.0) {
        omega_ci = Zi * ee * Bbar * B0 / Mi;
        omega_ce = ratio_pe * AA * ee * Bbar * B0 / Mi;

        kappa_perp_i = 2.0 * vth_i * vth_i * nu_i / (omega_ci * omega_ci); // * 1.e4;
        kappa_perp_e = 4.7 * vth_e * vth_e * nu_e / (omega_ce * omega_ce); // * 1.e4;

        output.write("\tion perp thermal conductivity: %e -> %e [m^2/s]\n", min(kappa_perp_i), max(kappa_perp_i));
        output.write("\telectron perp thermal conductivity: %e -> %e [m^2/s]\n", min(kappa_perp_e), max(kappa_perp_e));

        output.write("\tnormalized perp ion thermal conductivity: %e -> %e\n", min(kappa_perp_i * Tipara1), max(kappa_perp_i * Tipara1));
        output.write("\tnormalized perp electron thermal conductivity: %e -> %e\n", min(kappa_perp_e * Tepara1), max(kappa_perp_e * Tepara1));

        kappa_perp_i_fl = q_alpha * vth_i * q95 * Lbar; // * 1.e4;
        kappa_perp_e_fl = q_alpha * vth_e * q95 * Lbar; // * 1.e4;

        if (fluxlimit) {
            kappa_perp_i *= kappa_perp_i_fl / (kappa_perp_i + kappa_perp_i_fl);
            kappa_perp_e *= kappa_perp_e_fl / (kappa_perp_e + kappa_perp_e_fl);
        }
        kappa_perp_i *= Tipara1 * N0;
        output.write("\tUsed normalized ion perp thermal conductivity: %e -> %e\n", min(kappa_perp_i), max(kappa_perp_i));
        // kappa_perp_i.applyBoundary();
        // mesh->communicate(kappa_perp_i);
        kappa_perp_e *= Tepara1 * Ne0;
        output.write("\tUsed normalized electron perp thermal conductivity: %e -> %e\n", min(kappa_perp_e), max(kappa_perp_e));
        // kappa_perp_e.applyBoundary();
        // mesh->communicate(kappa_perp_e);

        dump.add(kappa_perp_i, "kappa_perp_i", 1);
        dump.add(kappa_perp_e, "kappa_perp_e", 1);
    }

    epsilon = minor_radius / major_radius;
    if (neoclassic_i) {
        rho_i = 1.02e-4 * sqrt(AA * Ti0 * Tibar) / B0 / Bbar / Zi;
        // Dri_neo = (1.+1.6*q95)*(1.+Tau_ie)*nu_i*rho_i*rho_i;
        xii_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_i * (rho_i ^ 2.);

        // output.write("\tion neoclassic particle transport: %e -> %e [m^2/s]\n", min(Dri_neo), max(Dri_neo));
        output.write("\tion neoclassic heat transport: %e -> %e [m^2/s]\n", min(xii_neo), max(xii_neo));

        // Dri_neo *= 3./2.*Tipara1;
        xii_neo *= Tipara1;

        // Dri_neo.applyBoundary();
        // xii_neo.applyBoundary();

        // output.write("\tNormalized ion neoclassic particle transport: %e -> %e [m^2/s]\n", min(Dri_neo), max(Dri_neo));
        output.write("\tNormalized ion neoclassic heat transport: %e -> %e [m^2/s]\n", min(xii_neo), max(xii_neo));

        // dump.add(Dri_neo, "Dri_neo", 1);
        dump.add(xii_neo, "xii_neo", 1);
        dump.add(heatf_neo_i, "heatf_neo_i", 1);
    }

    if (neoclassic_e) {
        rho_e = 2.38e-6 * sqrt(Te0 * Tebar) / B0 / Bbar;
        xie_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_e * (rho_e ^ 2.);

        output.write("\telectron neoclassic heat transport: %e -> %e [m^2/s]\n", min(xie_neo), max(xie_neo));

        xie_neo *= Tepara1;
        // xie_neo.applyBoundary();

        output.write("\tNormalized electron neoclassic heat transport: %e -> %e [m^2/s]\n", min(xie_neo), max(xie_neo));

        dump.add(xie_neo, "xie_neo", 1);
        dump.add(heatf_neo_e, "heatf_neo_e", 1);
    }

    /**************** CALCULATE METRICS ******************/

    mesh->g11 = (Rxy * Bpxy) ^ 2;
    mesh->g22 = 1.0 / (hthe ^ 2);
    mesh->g33 = (I ^ 2) * mesh->g11 + (B0 ^ 2) / mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -I * mesh->g11;
    mesh->g23 = -Btxy / (hthe * Bpxy * Rxy);

    mesh->J = hthe / Bpxy;
    mesh->Bxy = B0;

    mesh->g_11 = 1.0 / mesh->g11 + ((I * Rxy) ^ 2);
    mesh->g_22 = (B0 * hthe / Bpxy) ^ 2;
    mesh->g_33 = Rxy * Rxy;
    mesh->g_12 = Btxy * hthe * I * Rxy / Bpxy;
    mesh->g_13 = I * Rxy * Rxy;
    mesh->g_23 = Btxy * hthe * Rxy / Bpxy;

    mesh->geometry(); // Calculate quantities from metric tensor

    // Set B field vector

    B0vec.covariant = false;
    B0vec.x = 0.;
    B0vec.y = Bpxy / hthe;
    B0vec.z = 0.;

    // Set V0vec field vector

    V0vec.covariant = false;
    V0vec.x = 0.;
    V0vec.y = Vp0 / hthe;
    V0vec.z = Vt0 / Rxy;

    // Set V0eff field vector

    V0eff.covariant = false;
    V0eff.x = 0.;
    V0eff.y = -(Btxy / (B0 * B0)) * (Vp0 * Btxy - Vt0 * Bpxy) / hthe;
    V0eff.z = (Bpxy / (B0 * B0)) * (Vp0 * Btxy - Vt0 * Bpxy) / Rxy;

    /**************** SET VARIABLE LOCATIONS *************/

    P.setLocation(CELL_CENTRE);
    U.setLocation(CELL_CENTRE);
    phi.setLocation(CELL_CENTRE);
    if (zonal_flow < 0.) {
        phiDC.setLocation(CELL_CENTER);
        VortDC.setLocation(CELL_CENTER);
        TeDC.setLocation(CELL_CENTER);
        phiDC.setBoundary("phiDC");
    }
    Psi.setLocation(CELL_YLOW);
    Apar.setLocation(CELL_YLOW);
    if (emass)
        Ajpar.setLocation(CELL_YLOW);
    Jpar.setLocation(CELL_YLOW);

    Ni.setLocation(CELL_YLOW);
    Ti.setLocation(CELL_CENTRE);
    Te.setLocation(CELL_CENTRE);

    Vipar.setLocation(CELL_YLOW);
    Vepar.setLocation(CELL_YLOW);
    Pi.setLocation(CELL_CENTRE);
    Pe.setLocation(CELL_CENTRE);

    N_tmp.setLocation(CELL_CENTRE);
    if (impurity_prof)
        Ne_tmp.setLocation(CELL_CENTRE);
    if (nonlinear) {
        Ti_tmp.setLocation(CELL_CENTRE);
        Te_tmp.setLocation(CELL_CENTRE);
    }

    Pe.setBoundary("P");
    Pi.setBoundary("P");  //initialize sourcing profile

    NiSource = 0.0;
    if (NiAmp > 0.) {
        NiSource.setLocation(CELL_YLOW);
        NiSource.setBoundary("Ni");
//  NiSource.applyBoundary();
//  mesh->communicate(NiSource);
        for (jz=0;jz<mesh->ngz;jz++) {
            for (jy=0;jy<mesh->ngy;jy++) {
                indy = mesh->YGLOBAL(jy);
                for (jx=0;jx<mesh->ngx;jx++) {
                    indx = mesh->XGLOBAL(jx);
                    if (mag_config == 1 || mag_config == 2)
                        NiSource[jx][jy][jz]=NiAmp*exp(-((indx-NiLoc)*(indx-NiLoc)/(2.*NiSig*NiSig)));
                    if (mag_config == 3) {
                        if ( (indy>jysep1) && (indy<=jysep2) )
                            NiSource[jx][jy][jz]=NiAmp*exp(-((indx-NiLoc)*(indx-NiLoc)/(2.*NiSig*NiSig)));
                    }
                    if (mag_config == 4) {
//          output.write("%i,%i,%i,%i,%i,%i\n",mag_config,jx,jy,jz,indx,indy);
                        if ( ((indy>jysep1) && (indy<=jysep2_1)) || ((indy>jysep1_2) && (indy<=jysep2)) )
                            NiSource[jx][jy][jz]=NiAmp*exp(-((indx-NiLoc)*(indx-NiLoc)/(2.*NiSig*NiSig)));
                    }
                }
            }
        }
    }

#if DEBUG_6F>0
    output.write("NiSource initialization finished\n");
#endif
    NiSource.applyBoundary();
    mesh->communicate(NiSource);
    SAVE_ONCE(NiSource);

// for simply test purpose
    TeSource=NiSource;
    TiSource=NiSource;

    if (neutral) {
// setup variable location and boundaries
        Nn.setLocation(CELL_YLOW);
        lNn.setLocation(CELL_YLOW);
        Vn.setLocation(CELL_YLOW);
        Sn.setLocation(CELL_YLOW);
        Sv.setLocation(CELL_YLOW);
        Sn_ext.setLocation(CELL_YLOW);

        term1.setLocation(CELL_YLOW);
        term2.setLocation(CELL_YLOW);
        term3.setLocation(CELL_YLOW);
        term4.setLocation(CELL_YLOW);
        term5.setLocation(CELL_YLOW);

        // zhu: ideally Nn and Vn shall have the neumann bc at radial boundaries
        //Nn.setBoundary("Ni");
        //lNn.setBoundary("Ni");
        //Vn.setBoundary("Vipar");
        SOLVE_FOR(lNn);
        SOLVE_FOR(Vn);

        lNn = 0.;
// initial neutral profile with tanh function, could be changed easily
        for (jz=0;jz<mesh->ngz;jz++) {
            for (jy=0;jy<mesh->ngy;jy++) {
                indy = mesh->YGLOBAL(jy);
                for (jx=0;jx<mesh->ngx;jx++) {
                    indx = mesh->XGLOBAL(jx);
                    if (mag_config == 1 || mag_config == 2)
                        lNn[jx][jy][jz]=0.4*(tanh(-(indx-NnLoc)/NnSig)+1.);
                    if (mag_config == 3) {
                        if ( (indy>jysep1) && (indy<=jysep2) )
                            lNn[jx][jy][jz]=0.4*(tanh(-(indx-NnLoc)/NnSig)+1.);
                        //BoutReal SigNy = Grid_NY/2.;
                        //lNn[jx][jy][jz]=(exp(-indy*indy/SigNy/SigNy) \
			    +exp(-(Grid_NY-1-indy)*(Grid_NY-1-indy)/SigNy/SigNy)) \
			    *exp(-(Grid_NX-1-indx)*(Grid_NX-1-indx)/NnSig/NnSig);
                    }
                    if (mag_config == 4) {
                        if ( ((indy>jysep1) && (indy<=jysep2_1)) || ((indy>jysep1_2) && (indy<=jysep2)) )
                            lNn[jx][jy][jz]=0.4*(tanh(-(indx-NnLoc)/NnSig)+1.);
                    }
                }
            }
        }
        lNn.applyBoundary();
        mesh->communicate(lNn);

        Nn = NnAmp*(1.-lNn);
        //Nn = NnAmp;//*(lNn+1.0);
        Nn.applyBoundary();
        mesh->communicate(Nn);
        output.write("Max and min of normalized Nn = %e, %e.\n", max(Nn), min(Nn));
        SAVE_ONCE(Nn);
        lNn = log(Nn);
        output.write("Max and min of normalized ln(Nn) = %e, %e.\n", max(lNn), min(lNn));
        Vn = 0.;
        Sn_ext = 0.;

//    SOLVE_FOR(lNn);
//    SOLVE_FOR(Vn);

#if DEBUG_6F>0
        output.write("Neutral initialization finished\n");
#endif
    }

    /**************** SET EVOLVING VARIABLES *************/

    // Tell BOUT which variables to evolve
    SOLVE_FOR(U);
    SOLVE_FOR(Ni);
    SOLVE_FOR(Ti);
    SOLVE_FOR(Te);

    if (emass) {
        output.write("Solving for Psi, Differentiating to get jpar\n");
        SOLVE_FOR(Ajpar);
    } else {
        if (evolve_psi) {
            output.write("Solving for Psi, Differentiating to get jpar\n");
            SOLVE_FOR(Psi);
        } else {
            output.write("Solving for Apar, Differentiating to get jpar\n");
            SOLVE_FOR(Apar);
        }
    }
    dump.add(Jpar, "jpar", 1);
    if (output_P)
        dump.add(P, "P", 1);
    if (output_Vepar)
        dump.add(Vepar, "Vepar", 1);

    if (parallel_lagrange) {
        // Evolving the distortion of the flux surfaces (Ideal-MHD only!)

        bout_solve(Xip_x, "Xip_x");
        bout_solve(Xip_z, "Xip_z");

        bout_solve(Xim_x, "Xim_x");
        bout_solve(Xim_z, "Xim_z");
    }

    if (parallel_project) {
        // Add Xi to the dump file
        dump.add(Xip_x, "Xip_x", 1);
        dump.add(Xip_z, "Xip_z", 1);

        dump.add(Xim_x, "Xim_x", 1);
        dump.add(Xim_z, "Xim_z", 1);
    }

    if (compress0) {
        SOLVE_FOR(Vipar);
        if (!restarting)
            Vipar = 0.0;
    }

    // solver->setPrecon(precon);

    if (phi_constraint) {
        // Implicit Phi solve using IDA

        if (!bout_constrain(phi, C_phi, "phi")) {
            output.write("ERROR: Cannot constrain. Run again with phi_constraint=false\n");
            bout_error("Aborting.\n");
        }

    } else {
        // Phi solved in RHS (explicitly)
        dump.add(phi, "phi", 1);
    }

    /*  bxgradb = B0vec ^ Grad(B0)/B0/B0;
    bxgradb.toCovariant();
    dump.add(bxgradb.x, "bcurvx",0);
    dump.add(bxgradb.y, "bcurvy",0);
    dump.add(bxgradb.z, "bcurvz",0);
  */

    // Diamagnetic phi0
    if (diamag && diamag_phi0) {
        if (n0_p0_0p3)
            // n0 ~ p0^0.3
            phi0 = -1 / 0.7 * Upara0 * Pi0 / N0;
        else
            // const density
            phi0 = -Upara0 * Pi0 / N0;

        mesh->communicate(phi0);
        // phi0.applyBoundary('phi');

        if (experiment_Er) {
            if (diamag_er) {
                // get Er0 from grid file
                Field2D Er_tmp;
                mesh->get(Er_tmp, "E_r");
                Er0.x = Er0_factor * Er_tmp / (Bbar * Va * Rxy * Bpxy);
                Er0.y = 0.;
                Er0.z = 0.;

                mesh->communicate(Er0);
                Er0.x.applyBoundary();
                Er0.y.applyBoundary();
                Er0.z.applyBoundary();

                Er0_dia = Upara0 * Grad(Pi0) / N0;
                mesh->communicate(Er0_dia);
                Er0_dia.x.applyBoundary();
                Er0_dia.y.applyBoundary();
                Er0_dia.z.applyBoundary();
                mesh->communicate(Er0_dia);

                Er0_net = Er0 - Er0_dia;
            } else { // WARNING: deprecated
                // get phi0 from grid file
                // TODO: U0_net in the end??
                mesh->get(phi0, "Phi_0");
                phi0 /= Bbar * Lbar * Va;
                // BoutReal N0tmp = max(N0,true);
                // phi0_net =  phi0 + Upara0*Pi0/N0tmp;
                // U0_net = Delp2(phi0_net) * N0tmp/B0;
                U0_net = N0 / B0 * (Delp2(phi0) + Upara0 * Delp2(Pi0) / N0);
                mesh->communicate(U0_net);
                U0_net.applyBoundary("dirichlet");
            }
        } else {
            if (diamag_er) {
                // phi0 = -Upara0*Pi0/max(N0,true);
                Er0_dia = Upara0 * Grad(Pi0) / N0;
                mesh->communicate(Er0_dia);
                Er0_dia.x.applyBoundary();
                Er0_dia.y.applyBoundary();
                Er0_dia.z.applyBoundary();
                mesh->communicate(Er0_dia);

                Er0 = Er0_dia;
                mesh->communicate(Er0);
                Er0.x.applyBoundary();
                Er0.y.applyBoundary();
                Er0.z.applyBoundary();

                Er0_net = 0;
                // SAVE_ONCE(Er0);

            } else
                // Stationary equilibrium plasma. ExB velocity balances diamagnetic drift
            {
                Er0_dia = -Grad(phi0);
                mesh->communicate(Er0_dia);
                Er0_dia.x.applyBoundary();
                Er0_dia.y.applyBoundary();
                Er0_dia.z.applyBoundary();

                Er0 = Er0_dia;
                mesh->communicate(Er0);
                (Er0.x).applyBoundary();
                (Er0.y).applyBoundary();
                (Er0.z).applyBoundary();

                Er0_net = 0;
            }
        }
    } else {
        phi0 = 0.;
        Er0 = 0;
        Er0_net = 0;
        Er0_dia = 0;
        Ve0 = 0;
        Ve0_net = 0;
        Ve0_dia = 0;
    }

    // Er0_net, Ve0, Ve0_dia, Ve0_net, U0_net
    // applyBoundary(), communicate
    mesh->communicate(Er0_net);
    Er0_net.x.applyBoundary();
    Er0_net.y.applyBoundary();
    Er0_net.z.applyBoundary();

    Ve0 = (Er0 ^ B0vec) / (B0 * B0);
    Ve0.setBoundary("Vipar");
    mesh->communicate(Ve0);
    Ve0.x.applyBoundary();
    Ve0.y.applyBoundary();
    Ve0.z.applyBoundary();

    Ve0_dia = (Er0_dia ^ B0vec) / (B0 * B0);
    Ve0_dia.setBoundary("Vipar");
    mesh->communicate(Ve0_dia);
    Ve0_dia.x.applyBoundary();
    Ve0_dia.y.applyBoundary();
    Ve0_dia.z.applyBoundary();

    Ve0_net = Ve0 - Ve0_dia;
    Ve0_net.setBoundary("Vipar");
    mesh->communicate(Ve0_net);
    Ve0_net.x.applyBoundary();
    Ve0_net.y.applyBoundary();
    Ve0_net.z.applyBoundary();

    U0_net = B0vec * Curl(N0 * Ve0_net) / B0;
    mesh->communicate(U0_net);
    U0_net.applyBoundary("dirichlet");

    // Add some equilibrium quantities and normalisations
    // everything needed to recover physical units
    SAVE_ONCE2(J0, P0);
    SAVE_ONCE4(density, Lbar, Bbar, Tbar);
    SAVE_ONCE3(Tibar, Tebar, Nbar);
    SAVE_ONCE2(Va, B0);
    SAVE_ONCE4(Ti0, Te0, N0, Ne0);
    if (impurity_prof)
    SAVE_ONCE3(N_imp0, T_imp0, P_imp0);
    if (diamag) {
        SAVE_ONCE3(phi0, Er0_dia, Ve0_dia);
        if (experiment_Er)
        SAVE_ONCE3(Er0, Ve0, U0_net);
    }

    /////////////// CHECK VACUUM ///////////////////////
    // In vacuum region, initial vorticity should equal zero

    ubyn.setLocation(CELL_CENTRE);
    ubyn.setBoundary("U");

    density_tmp = N0 + A_imp / AA * N_imp0;

    if (zonal_flow < 0.) lapDC = new LaplaceXY2(mesh);

    if (!restarting) {
        // Only if not restarting: Check initial perturbation

        // Set U to zero where P0 < vacuum_pressure
        U = where(P0 - vacuum_pressure, U, 0.0);

        // Field2D lap_temp = 0.0;
        // TODO: diamag in U?
        if (impurity_prof) {
            Field2D logn0 = laplace_alpha * density_tmp;
            ubyn = U * B0 / density_tmp;
            // Phi should be consistent with U
            if (laplace_alpha <= 0.0)
                phi = invert_laplace(ubyn, phi_flags, NULL);
            else
                phi = invert_laplace(ubyn, phi_flags, NULL, &logn0, NULL);
        }
            // Ni_tot = N0+A_imp*N_imp0;
            // else
            // Ni_tot = N0;
        else {
            Field2D logn0 = laplace_alpha * N0;
            // Field3D ubyn;
            Field3D Ntemp;
            Ntemp = max(N0, true);
            // ubyn = U*B0/Ntemp;
            ubyn = U * B0 / N0;
            // Phi should be consistent with U
            if (laplace_alpha <= 0.0)
                phi = invert_laplace(ubyn, phi_flags, NULL);
            else
                phi = invert_laplace(ubyn, phi_flags, NULL, &logn0, NULL);
        }
    }

    if ((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
        output.write("Sheath Boundary conditions applied.\n");
        if (output_SBC) {
            dump.add(c_se, "c_se", 1);
            dump.add(q_si, "q_si", 1);
            dump.add(q_se, "q_se", 1);
        }

        const_cse = sqrt(KB * Tebar * eV_K / Mi);
        vth_et.setLocation(CELL_YLOW);
        c_set.setLocation(CELL_YLOW);
        c_se.setLocation(CELL_YLOW);
        q_se.setLocation(CELL_YLOW);
        q_si.setLocation(CELL_YLOW);

        c_se0 = sqrt(abs(Tau_ie * Ti0 + Te0));
        c_se0 *= const_cse;
        vth_e0 = 4.19e5 * sqrt(Te0 * Tebar);

        // output << max(c_se0, true) << "\t" << max(vth_e0, true) << endl;
        if (output_SBC)
            dump.add(Jpar_sh, "Jpar_sh", 1);
        Jpar_sh.setLocation(CELL_YLOW);
        Jpar_sh0 = Ne0 * Nbar * density * ee;
        Jpar_sh0 *= c_se0 - vth_e0 / (2.0 * sqrt(PI)) * exp(-ee * (phi0 * Va * Lbar * Bbar) / (KB * Te0 * Tebar * eV_K));
        // Jpar_sh0 *= c_se0 -  vth_e0/(2.0*sqrt(PI)) * ( 1. - ee*(phi0*Va*Lbar*Bbar)/(KB*Te0*Tebar*eV_K) );
        if (output_SBC)
            dump.add(phi_sh, "phi_sh", 1);
        phi_sh.setLocation(CELL_YLOW);
        phi_sh.setBoundary("phi");
        // phi_sh0 = -Te0*Tebar;
        // phi_sh0 *= log( 2.*sqrt(PI)*(c_se0-J0*B0*Bbar/(MU0*Lbar)/(Ne0*Nbar*density*ee))/vth_e0 );
    }

    if (BScurrent || neo_resist) {
        nu_estar.setLocation(CELL_YLOW);
        nu_estar = nu_e * q95 * major_radius * Lbar / (vth_e) / pow(epsilon, 1.5);
        output.write("Normalized electron collisionality: nu_e* = %e~%e\n", min(nu_estar, true), max(nu_estar, true));
        ft = BS_ft(3000);
        output.write("modified collisional trapped particle fraction: ft = %e~%e\n", min(ft, true), max(ft, true));
    }
    if (neo_resist) {
        f33 = ft / (1. + (0.55 - 0.1 * ft) * sqrt(nu_estar) + 0.45 * (1. - ft) * nu_estar / Zeff / sqrt(Zeff));
        cond_neo = F33(f33);
        output.write("Neoclassical resistivity used\n");
        output.write("\tlocal max eta is %e~%e * eta_Sp\n", 1/max(cond_neo, true), 1/min(cond_neo, true));
        eta /= cond_neo;
    }
    if (BScurrent) {
        // TODO: check
        Field3D L31, L32, L34;
        Field3D f31, f32ee, f32ei, f34;
        Field3D BSal0, BSal;
        Jpar_BS0.setLocation(CELL_YLOW);
        nu_istar.setLocation(CELL_YLOW);
        Jpar_BS0.setBoundary("J");

//        nu_istar = 100. * nu_i * q95 * Lbar / (vth_i) / pow(epsilon, 1.5);
        nu_istar = nu_i * q95 * major_radius * Lbar / (vth_i) / pow(epsilon, 1.5);
        // nu_estar = 0.012 * N0*Nbar*density/1.e20*Zi*Zi*q95*Lbar/(Te0*Tebar/1000. * pow(Aratio, 1.5));
        // nu_istar = 0.012 * N0*Nbar*density/1.e20*Zi*q95*Lbar/(Ti0*Tibar/1000. * pow(Aratio, 1.5));
        output.write("Bootstrap current is included:\n");
        output.write("Normalized ion collisionality: nu_i* = %e\n", max(nu_istar));
        f31 = ft / (1. + (1. - 0.1 * ft) * sqrt(nu_estar) + 0.5 * (1. - ft) * nu_estar / Zi);
        f32ee = ft / (1. + 0.26 * (1. - ft) * sqrt(nu_estar) + 0.18 * (1. - 0.37 * ft) * nu_estar / sqrt(Zi));
        f32ei = ft / (1. + (1. + 0.6 * ft) * sqrt(nu_estar) + 0.85 * (1. - 0.37 * ft) * nu_estar * (1. + Zi));
        f34 = ft / (1. + (1. - 0.1 * ft) * sqrt(nu_estar) + 0.5 * (1. - 0.5 * ft) * nu_estar / Zi);

        L31 = F31(f31);
        L32 = F32ee(f32ee) + F32ei(f32ei);
        L34 = F31(f34);

        BSal0 = -(1.17 * (1. - ft)) / (1. - 0.22 * ft - 0.19 * ft * ft);
        BSal = (BSal0 + 0.25 * (1 - ft * ft) * sqrt(nu_istar)) / (1. + 0.5 * sqrt(nu_istar)) + 0.31 * nu_istar * nu_istar * ft * ft * ft * ft * ft * ft;
        BSal *= 1. / (1. + 0.15 * nu_istar * nu_istar * ft * ft * ft * ft * ft * ft);

        Jpar_BS0 = L31 * DDX(P0) / Pe0 + L32 * DDX(Te0) / Te0 + L34 * DDX(Ti0) / (Zi * Te0) * BSal;
        Jpar_BS0 *= Field3D(-Rxy * Btxy * Pe0 / (B0 * B0) * (MU0 * KB * Nbar * density * Tebar * eV_K) / (Bbar * Bbar));

        mesh->communicate(Jpar_BS0);
        Jpar_BS0.applyBoundary();

        dump.add(Jpar_BS0, "jpar_BS0", 0);
        dump.add(nu_estar, "nu_estar", 0);
        dump.add(nu_istar, "nu_istar", 0);
    }

    if (radial_diffusion) {
        diffusion_coef_Hmode0 /= Lbar * Lbar / Tbar;
        diffusion_coef_Hmode1 /= Lbar * Lbar / Tbar;

        ddx_ni.setLocation(CELL_YLOW);
        diff_radial.setLocation(CELL_YLOW);
        ddx_ni.setBoundary("Ni");
        diff_radial.setBoundary("Ni");

        ddx_n0 = DDX(N0);
        ddx_n0 = -field_larger(-ddx_n0, Low_limit);
        diff_radial = Field3D(diffusion_coef_Hmode0);

        dump.add(diff_radial, "diff_radial", 1);
    }

    /************** SETUP COMMUNICATIONS **************/

    comms.add(U);
    // comms.add(phi);
    comms.add(Ni);
    comms.add(Ti);
    comms.add(Te);
    if (!emass) {
        if (evolve_psi)
            comms.add(Psi);
        else
            comms.add(Apar);
    }
    else
        comms.add(Ajpar);

    if (compress0) {
        comms.add(Vipar);
        Vepar.setBoundary("Vipar");
    }

    if (neutral) {
        comms.add(lNn);
        comms.add(Vn);
    }

    if (hyperdiff_par_u4 > 0.0 || hyperdiff_perp_u4 > 0.0)
        tmpU2.setBoundary("U");

    if (hyperdiff_par_apar4 > 0.0 || hyperdiff_perp_apar4 > 0.0)
        tmpA2.setBoundary("J");

    if (hyperdiff_par_n4 > 0.0 || hyperdiff_perp_n4 > 0.0)
        tmpN2.setBoundary("Ni");

    if (hyperdiff_par_ti4 > 0.0 || hyperdiff_perp_ti4 > 0.0)
        tmpTi2.setBoundary("Ti");

    if (hyperdiff_par_te4 > 0.0 || hyperdiff_perp_te4 > 0.0)
        tmpTe2.setBoundary("Te");

    if (hyperdiff_par_v4 > 0.0 || hyperdiff_perp_v4 > 0.0)
        tmpVp2.setBoundary("Vipar");

    phi.setBoundary("phi"); // Set boundary conditions

    P.setBoundary("P");
    Jpar.setBoundary("J");
    Jpar2.setBoundary("J");

    F2D_tmp = 0.;

    if (fakerun) {
        // Vector3D Btilde;
        Btilde.setBoundary("P");
        Btilde.setLocation(CELL_YLOW);
        dump.add(Btilde, "btilde", 1);
    }

    /*  if (pos_sink_zf)
      {
        pos_sink_ti = tanhxl_core(filter_position_ti);
        pos_sink_te = tanhxl_core(filter_position_te);
        pos_sink_ni = tanhxl_core(filter_position_ni);
  // dump.add(pos_sink_ti, "pos_sink_ti", 0);
      }
  */
    return 0;
}

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
    Field3D result;

    if (parallel_lagrange || parallel_project) {
        // Moving stencil locations

        // Field3D fp, fm; // Interpolated on + and - y locations

        fp = interpolate(f, Xip_x, Xip_z);
        fm = interpolate(f, Xim_x, Xim_z);

        result.allocate();
        for (jx = 0; jx < mesh->ngx; jx++)
            for (jy = 1; jy < mesh->ngy - 1; jy++)
                for (jz = 0; jz < mesh->ngz - 1; jz++) {
                    result[jx][jy][jz] = (fp[jx][jy + 1][jz] - fm[jx][jy - 1][jz]) / (2. * mesh->dy[jx][jy] * sqrt(mesh->g_22[jx][jy]));
                }
    } else {
        if (parallel_lr_diff) {
            // Use left/right biased stencils. NOTE: First order only!
            if (loc == CELL_YLOW) {
                result = Grad_par_CtoL(f);
            } else
                result = Grad_par_LtoC(f);
        } else
            result = Grad_par(f, loc);

        if (nonlinear) {
            if (evolve_psi)
                result -= bracket(Psi, f, bm_mag) * B0;
            else
                result -= bracket(Apar, f, bm_mag);
        }
    }

    return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t) {

    if (fakerun) {
        // output.write("I see you 0.0!\n");//xia
        int rank;
        char filename[200];
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        sprintf(filename, "%sBOUT.dmp.%d.nc", path.c_str(), rank);
        timestep = (int)floor(t);

        if (false) {
            output.write("timestep = %d from process %d\n", timestep, rank);
            output.write("Reading from file %s\n", filename);
        }

        // output.write("ngx=%i, ngy=%i, ngz=%i\n",mesh->ngx, mesh->ngy, mesh->ngz);
        DataFormat *file = data_format(filename);
        file->openr(filename);
        file->setRecord(timestep);
        Psi.allocate();
        file->read_rec(**(Psi.getData()), "Psi", mesh->ngx, mesh->ngy, mesh->ngz);
        Apar.allocate();
        file->read_rec(**(Apar.getData()), "Apar", mesh->ngx, mesh->ngy, mesh->ngz);
        Vipar.allocate();
        file->read_rec(**(Vipar.getData()), "Vipar", mesh->ngx, mesh->ngy, mesh->ngz);
        Ni.allocate();
        file->read_rec(**(Ni.getData()), "Ni", mesh->ngx, mesh->ngy, mesh->ngz);
        Te.allocate();
        file->read_rec(**(Te.getData()), "Te", mesh->ngx, mesh->ngy, mesh->ngz);
        Ti.allocate();
        file->read_rec(**(Ti.getData()), "Ti", mesh->ngx, mesh->ngy, mesh->ngz);
        U.allocate();
        file->read_rec(**(U.getData()), "U", mesh->ngx, mesh->ngy, mesh->ngz);
        // Vexp_ave_tmp.allocate();
        // file->read_rec(**(Vexp_ave_tmp.getData()), "Vexp_ave_tmp", mesh->ngx, mesh->ngy, mesh->ngz);

        // mesh->communicate(Vexp_ave_tmp);
        // Vexp_ave_tmp.applyBoundary();
        // mesh->communicate(Psi);
        file->close();

        // output.write("I see you 0.1!\n");//xia

        /*      Btilde = Grad(Psi) ^ B0vec;
          mesh->communicate(Btilde);
          Btilde.applyBoundary();
          Btilde.toContravariant();*/
        // output.write("I see you 0.2!\n");//xia
    }

    mesh->communicate(comms);

    /*      Btilde = Grad(Psi) ^ B0vec;
      mesh->communicate(Btilde);
      Btilde.applyBoundary();
      Btilde.toContravariant();*/
#if DEBUG_6F>0
    output.write("I see you 0!\n");//xia
#endif

    // Perform communications
    // mesh->communicate(comms);

    if (pos_filter) {
        Ti = lowPass_pos2(Ti, Ti);
        Te = lowPass_pos2(Te, Te);
        Ni = lowPass_pos2(Ni, Ni);
    }

    // Inversion
    Pi = Ni * Ti0 + N0 * Ti;
    if (nonlinear)
        Pi += Ni * Ti;
    mesh->communicate(Pi);

    Pe = Zi * Ni * Te0 + Ne0 * Te;
    if (nonlinear)
        Pe += Zi * Ni * Te;
    mesh->communicate(Pe);

    P = Tau_ie * Pi + Pe;
    mesh->communicate(P);

    if (nonlinear && pos_filter) {
        Pi = lowPass_pos2(Pi, Pi);
        Pe = lowPass_pos2(Pe, Pe);
        P = lowPass_pos2(P, P);
    }
    if (nonlinear && pos_filter2) {
        Pi = lowPass_pos(Pi, position_tmpi);
        Pe = lowPass_pos(Pe, position_tmpe);
        P = lowPass_pos(P, position_tmp);
    }

    if (nonlinear && (pos_filter_zf || (pos_sink_zf > 0.))) {
        Pi = sink_zonal_core(Pi, position_tmpi);
        Pe = sink_zonal_core(Pe, position_tmpe);
        P = sink_zonal_core(P, position_tmp);
    }
    if (low_pass_z > 0.) {
        Pi = lowPass(Pi, low_pass_z, zonal_bkgd);
        Pe = lowPass(Pe, low_pass_z, zonal_bkgd);
        P = lowPass(P, low_pass_z, zonal_bkgd);
    }

    // output.write("I see you 0.1!\n");//xia

    if ((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
        SBC_Gradpar(U, 0.0, PF_limit, PF_limit_range);
    }

    // Field2D lap_temp=0.0;
    // Field2D logn0 = laplace_alpha * N0;
    // ubyn = U*B0/max(N0, true);
    if (impurity_prof) {
        ubyn = U * B0 / density_tmp;
        if (diamag) {
            ubyn -= Upara0 / density_tmp * Delp2(Pi);
            mesh->communicate(ubyn);
            ubyn.applyBoundary();
        }

        // output.write("I see you 0.2!\n");//xia

        // Invert laplacian for phi
        if (laplace_alpha <= 0.0)
            phi = invert_laplace(ubyn, phi_flags, NULL);
        else
            phi = invert_laplace(ubyn, phi_flags, NULL, &density_tmp, NULL);
        // output.write("I see you 0.21!\n");//xia
    } else {
        ubyn = U * B0 / N0;
        if (diamag) {
            ubyn -= Upara0 / N0 * Delp2(Pi);
            mesh->communicate(ubyn);
            ubyn.applyBoundary();
        }

        // output.write("I see you 0.2!\n");//xia

        // Invert laplacian for phi
        if (laplace_alpha <= 0.0)
            phi = invert_laplace(ubyn, phi_flags, NULL);
        else
            phi = invert_laplace(ubyn, phi_flags, NULL, &N0, NULL);
        // output.write("I see you 0.21!\n");//xia
    }

    if (mask_phi_x) {
        phi *= mask_px1d;
    }

    mesh->communicate(phi);

    if (low_pass_z > 0.)
        phi = lowPass(phi, low_pass_z, zonal_flow);

//  if (low_pass_z > 0. || zonal_flow < 0.)
//    phi = lowPass(phi, low_pass_z, 0); // now phi has no DC component

    if (zonal_flow < 0.) {
        VortDC = getDC(ubyn);
        VortDC.applyBoundary("neumann");
        TeDC = getDC(Te);
        phiDC = 2.5*TeDC;
#if DEBUG_6F>0
        output.write("begin invert Laplacian\n");
#endif
        if (laplace_alpha <= 0.0)
            lapDC->setCoefs(0.0, 1.0, 1.0);
        else {
            if (impurity_prof)
                lapDC->setCoefs(0.0, density_tmp, 1.0);
            else
                lapDC->setCoefs(0.0, N0, 1.0);
        }
#if DEBUG_6F>0
        output.write("coefficients set\n");
#endif
        phiDC = lapDC->solve(VortDC, phiDC);
#if DEBUG_6F>0
        output.write("\tphiDC solved successfully!\n");
#endif
        phi += phiDC;
    }

    // Apply a boundary condition on phi for target plates
    phi.applyBoundary();
    mesh->communicate(phi);
    //  }

    if (emass) {
        // static Field2D acoeff;
        // static bool aset = false;

        if (!aset) // calculate Apar coefficient
            acoeff = -delta_e_inv * N0 * N0;
        aset = true;
        if (compress0) {
            Psi = invert_laplace(acoeff * Ajpar - gyroAlv * Vipar, apar_flags, &acoeff);
        } else {
            Psi = invert_laplace(acoeff * Ajpar, apar_flags, &acoeff);
        }
        mesh->communicate(Psi);
    }

    // output.write("I see you 0.3!\n");//xia

    // BoutReal N_tmp1;
    N_tmp1 = Low_limit;
    if (nonlinear) {
        N_tmp = field_larger(N0 + Ni, N_tmp1);
        if (impurity_prof)
            Ne_tmp = field_larger(Ne0 + Zi * Ni, N_tmp1);
        else
            Ne_tmp = Zi * N_tmp;
    }

    // BoutReal Te_tmp1, Ti_tmp1;
    Te_tmp1 = Low_limit;
    Ti_tmp1 = Low_limit;
    if (nonlinear) {
        Ti_tmp = field_larger(Ti0 + Ti, Ti_tmp1);
        Te_tmp = field_larger(Te0 + Te, Te_tmp1);
    }

    // output.write("I see you 1!\n");//xia

    if (!nonlinear && (parallel_viscous && compress0)) {
        pi_ci = -eta_i0 * 2. * (B0 ^ (-0.5)) * Grad_par(((B0 ^ (0.5)) * Vipar), CELL_YLOW);
        pi_ci -= eta_i0 * b0xcv * (Er0_net + Grad(phi)) / B0;
        mesh->communicate(pi_ci);
        pi_ci.applyBoundary();
    }

    // Transitions from 0 in core to 1 in vacuum
    if (nonlinear) {
        vac_mask = (1.0 - tanh(((P0 + P) - vacuum_pressure) / vacuum_trans)) / 2.0;
        // Update resistivity
        if (spitzer_resist || neo_resist) {
            // Use Spitzer formula
            eta = FZ * 1.03e-4 * Zeff * LnLambda * ((Te_tmp * Tebar) ^ (-1.5)); // eta in Ohm-m. ln(Lambda) = 20
            eta /= MU0 * Va * Lbar;
            // eta.applyBoundary();
            // mesh->communicate(eta);
        } else {
            eta = core_resist + (vac_resist - core_resist) * vac_mask;
        }

        if (impurity_prof)
            nu_e = 2.91e-6 * LnLambda * (Ne_tmp * Nbar * density / 1.e6) * ((Te_tmp * Tebar) ^ (-1.5));
        else
            nu_e = 2.91e-6 * LnLambda * (N_tmp * Nbar * density / 1.e6) * ((Te_tmp * Tebar) ^ (-1.5)); // nu_e in 1/S.
        // nu_e.applyBoundary();
        // mesh->communicate(nu_e);

        if (diffusion_par > 0.0 || diffusion_perp > 0.0 || parallel_viscous || neoclassic_i || neoclassic_e || neo_resist) {

            // xqx addition, begin
            // Use Spitzer thermal conductivities

            nu_i = 4.80e-8 * (Zi * Zi * Zi * Zi / sqrt(AA)) * LnLambda * (N_tmp * Nbar * density / 1.e6) * ((Ti_tmp * Tibar) ^ (-1.5)); // nu_i in 1/S.
            // nu_i.applyBoundary();
            // mesh->communicate(nu_i);

            vth_i = 9.79e3 * sqrt(Ti_tmp * Tibar / AA); // vth_i in m/S.
            // vth_i.applyBoundary();
            // mesh->communicate(vth_i);
            vth_e = 4.19e5 * sqrt(Te_tmp * Tebar); // vth_e in m/S.
            // vth_e.applyBoundary();
            // mesh->communicate(vth_e);
        }

        if (neo_resist) {
            nu_estar = nu_e * q95 * Lbar / (vth_e) / pow(epsilon, 1.5);
            f33 = ft / (1. + (0.55 - 0.1 * ft) * sqrt(nu_estar) + 0.45 * (1. - ft) * nu_estar / Zeff / sqrt(Zeff));
            cond_neo = F33(f33);
            eta /= cond_neo;
        }

        if (parallel_viscous && compress0) {
            eta_i0 = 0.96 * (Pi0 + Pi) * Tau_ie * nu_i * Tbar;
            pi_ci = -eta_i0 * 2. / sqrt(B0) * Grad_parP((sqrt(B0) * Vipar), CELL_YLOW);
            pi_ci -= eta_i0 * b0xcv * (Er0_net + Grad(phi)) / B0;
            mesh->communicate(pi_ci);
            pi_ci.applyBoundary();
        }

        if (diffusion_par > 0.0) {
//            kappa_par_i = 3.9 * vth_i * vth_i / nu_i; // * 1.e4;
//            kappa_par_e = 3.2 * vth_e * vth_e / nu_e; // * 1.e4;
//
//            kappa_par_i_fl = q_alpha * vth_i * q95 * Lbar; // * 1.e2;
//            kappa_par_e_fl = q_alpha * vth_e * q95 * Lbar; // * 1.e2;
//
//            if (fluxlimit) {
//                kappa_par_i *= kappa_par_i_fl / (kappa_par_i + kappa_par_i_fl);
//                kappa_par_e *= kappa_par_e_fl / (kappa_par_e + kappa_par_e_fl);
//            }
//            kappa_par_i *= diffusion_par * Tipara1 * N_tmp;
//            // kappa_par_i.applyBoundary();
//            // mesh->communicate(kappa_par_i);
//            kappa_par_e *= diffusion_par * Tepara1 * Ne_tmp;
//            // kappa_par_e.applyBoundary();
//            // mesh->communicate(kappa_par_e);

            kappa_par_i_sp = 3.9 * vth_i * vth_i / nu_i;
            kappa_par_e_sp = 3.2 * vth_e * vth_e / nu_e;
            kappa_par_i_sp *= Tipara1 * N_tmp;
            kappa_par_e_sp *= Tepara1 * Ne_tmp;
            if (fluxlimit) {
                kappa_par_i_fl = q_alpha * vth_i * q95 * major_radius * Lbar;
                kappa_par_e_fl = q_alpha * vth_e * q95 * major_radius * Lbar;
//                kappa_par_i_fl = q_alpha * vth_i * q95 * Lbar;
//                kappa_par_e_fl = q_alpha * vth_e * q95 * Lbar;
                kappa_par_i_fl *= Tipara1 * N_tmp;
                kappa_par_e_fl *= Tepara1 * Ne_tmp;
                kappa_par_i = diffusion_par * (kappa_par_i_sp * kappa_par_i_fl);
                kappa_par_i /= (kappa_par_i_sp + kappa_par_i_fl);
                kappa_par_e = diffusion_par * (kappa_par_e_sp * kappa_par_e_fl);
                kappa_par_e /= (kappa_par_e_sp + kappa_par_e_fl);
            } else {
                kappa_par_i = diffusion_par * kappa_par_i_sp;
                kappa_par_e = diffusion_par * kappa_par_e_sp;
            }
            if (kappa_par_i_const > 0.0) {
                kappa_par_i += kappa_par_i_const;
            }
            if (kappa_par_e_const > 0.0) {
                kappa_par_e += kappa_par_e_const;
            }
        }

        if (diffusion_perp > 0.0) {

            kappa_perp_i = 2.0 * vth_i * vth_i * nu_i / (omega_ci * omega_ci); // * 1.e4;
            kappa_perp_e = 4.7 * vth_e * vth_e * nu_e / (omega_ce * omega_ce); // * 1.e4;

            kappa_perp_i_fl = q_alpha * vth_i * q95 * Lbar; // * 1.e4;
            kappa_perp_e_fl = q_alpha * vth_e * q95 * Lbar; // * 1.e4;

            if (fluxlimit) {
                kappa_perp_i *= kappa_perp_i_fl / (kappa_perp_i + kappa_perp_i_fl);
                kappa_perp_e *= kappa_perp_e_fl / (kappa_perp_e + kappa_perp_e_fl);
            }
            kappa_perp_i *= Tipara1 * N_tmp;
            // kappa_perp_i.applyBoundary();
            // mesh->communicate(kappa_perp_i);
            kappa_perp_e *= Tepara1 * Ne_tmp;
            // kappa_perp_e.applyBoundary();
            // mesh->communicate(kappa_perp_e);
        }

        if (neoclassic_i) {
            rho_i = 1.02e-4 * sqrt(AA * Ti_tmp * Tibar) / B0 / Bbar / Zi;
            // Dri_neo = (1.+1.6*q95)*(1.+Tau_ie)*nu_i*rho_i*rho_i;
            xii_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_i * (rho_i ^ 2.);

            // Dri_neo *= 3./2.*Tipara1;
            xii_neo *= Tipara1;
            // Dri_neo.applyBoundary();
            // xii_neo.applyBoundary();
        }

        if (neoclassic_e) {
            rho_e = 2.38e-6 * sqrt(Te_tmp * Tebar) / B0 / Bbar;
            xie_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_e * rho_e * rho_e;

            xie_neo *= Tepara1;
            // xie_neo.applyBoundary();
        }
    }

    // update Landau parallel heat flux
/*  if (diffusion_par && Landau) {
    SBC_value_i = gamma_i_BC * Pi * c_set / (1.6 * N0  * vth_i);
    SBC_value_e = gamma_e_BC * Pe * c_set / (1.6 * Ne0 * vth_e);
    if (Landau_coll) {
      q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va *eiSign_kpar_wcoll(Ti, kappa_i, 2, nLorentzian, SBC_value_i);
      q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar_wcoll(Te, kappa_e, 2, nLorentzian, SBC_value_e);
    } else {
      q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va * iSign_kpar(Ti, kappa_0, 2, nLorentzian, SBC_value_i);
      q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar(Te, kappa_0, 2, nLorentzian, SBC_value_e);
    }
    mesh->communicate(q_par_i);
    mesh->communicate(q_par_e);
    q_par_i.applyBoundary();
    q_par_e.applyBoundary();
  }*/


    if (radial_diffusion && nonlinear) {
        ddx_ni = DDX(N_tmp);
        mesh->communicate(ddx_ni);
        ddx_ni.applyBoundary();

        for (jx = 0; jx < mesh->ngx; jx++) {
            for (jy = 0; jy < mesh->ngy; jy++) {
                for (jz = 0; jz < mesh->ngz; jz++) {
                    if (ddx_ni[jx][jy][jz] > -Low_limit && ddx_ni[jx][jy][jz] < 0.)
                        ddx_ni[jx][jy][jz] = -Low_limit;
                    else if (ddx_ni[jx][jy][jz] < Low_limit && ddx_ni[jx][jy][jz] >= 0.)
                        ddx_ni[jx][jy][jz] = Low_limit;
                }
            }
        }

        diff_radial = diffusion_coef_Hmode0 * ddx_n0 / ddx_ni;

        for (jx = 0; jx < mesh->ngx; jx++) {
            for (jy = 0; jy < mesh->ngy; jy++) {
                for (jz = 0; jz < mesh->ngz; jz++) {
                    if (diff_radial[jx][jy][jz] > diffusion_coef_Hmode1)
                        diff_radial[jx][jy][jz] = diffusion_coef_Hmode1;
                }
            }
        }
        diff_radial = nl_filter(diff_radial, 1);
    }

#if DEBUG_6F>0
    output.write("I see you 2!\n");//xia
#endif

    if (evolve_psi)
        Jpar = -B0 * Delp2(Psi);
    else
        Jpar = -Delp2(Apar);
    Jpar.applyBoundary();
    mesh->communicate(Jpar);

    if (jpar_bndry_width > 0) {
        // Zero j in boundary regions. Prevents vorticity drive
        // at the boundary

        for (jx = 0; jx < jpar_bndry_width; jx++)
            for (jy = 0; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz - 1; jz++) {
                    if (mesh->firstX())
                        Jpar[jx][jy][jz] = 0.0;
                    if (mesh->lastX())
                        Jpar[mesh->ngx - 1 - jx][jy][jz] = 0.0;
                }
    }

    // Smooth j in x
    if (smooth_j_x && !((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0))) {
        Jpar = smooth_x(Jpar);
        /*Jpar = smooth_x(Jpar);
    Jpar = smooth_x(Jpar);
    Jpar = smooth_x(Jpar);
    Jpar = smooth_x(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);*/
        mesh->communicate(Jpar);
        Jpar.applyBoundary();
    }

    if (mask_j_x) {
        Jpar *= mask_jx1d;
    }

    // TODO: check
    if ((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
        if (nonlinear) {
            c_set = sqrt(abs(Tau_ie * Ti_tmp + Te_tmp));
            c_se = c_set - c_se0 / const_cse;
            c_se *= const_cse / Va; // normalized
            c_set *= const_cse;     // not normalized, with unit
            vth_et = 4.19e5 * sqrt(Te_tmp * Tebar);
        } else {
            c_set = sqrt(abs(Tau_ie * Ti0 + Te0));
            c_se = 0.;
            c_se *= const_cse / Va; // normalized
            c_set *= const_cse;     // not normalized, with unit
            vth_et = 4.19e5 * sqrt(Te0 * Tebar);
        }

        /*phi_sh = -eta * Jpar / B0;

    if (eHall) {
      phi_sh += Psipara1 * Grad_parP(Pe, CELL_YLOW) / B0 / Ne0;
      phi_sh -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
    }
    if (thermal_force) {
      phi_sh += 0.71 * Psipara1 * Grad_parP(Te, CELL_YLOW) / B0;
      phi_sh -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
    }
    */

        if (evolve_psi) {
            phi_sh = -eta * Jpar / B0;

            if (eHall) {
                phi_sh += Psipara1 * Grad_parP(Pe, CELL_YLOW) / B0 / Ne0;
                phi_sh -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
            }
            if (thermal_force) {
                phi_sh += 0.71 * Psipara1 * Grad_parP(Te, CELL_YLOW) / B0;
                phi_sh -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
            }
        } else {
            phi_sh = -eta * Jpar;

            if (eHall) {
                phi_sh += Psipara1 * Grad_parP(Pe, CELL_YLOW) / Ne0;
                phi_sh -= Psipara1 * bracket(Apar, Pe0, bm_mag) / Ne0;
            }
            if (thermal_force) {
                phi_sh += 0.71 * Psipara1 * Grad_parP(Te, CELL_YLOW);
                phi_sh -= 0.71 * Psipara1 * bracket(Apar, Te0, bm_mag);
            }
        }
        mesh->communicate(phi_sh);
        phi_sh.applyBoundary();
        SBC_Gradpar(phi, phi_sh, PF_limit, PF_limit_range);

        if (nonlinear) {
            Jpar_sh = Ne_tmp * Nbar * density * ee;
            Jpar_sh *= c_set - vth_et / (2.0 * sqrt(PI)) * exp(-ee * ((phi + phi0) * Va * Lbar * Bbar) / (KB * Te_tmp * Tebar * eV_K));
            // Jpar_sh *= c_set -  vth_et/(2.0*sqrt(PI)) * (1. - ee*((phi+phi0)*Va*Lbar*Bbar)/(KB*Te_tmp*Tebar*eV_K) );
            Jpar_sh -= Jpar_sh0;
            Jpar_sh *= MU0 * Lbar / Bbar;
            SBC_Dirichlet(Jpar, Jpar_sh, PF_limit, PF_limit_range);

            /* phi_sh = -Te_tmp*Tebar;
         if (impurity_prof)
           phi_sh *= log( 2.0*sqrt(PI)*(c_set-(J0+Jpar)*B0*Bbar/(MU0*Lbar)/(Ne_tmp*Nbar*density*ee))/vth_et );
         else
           phi_sh *= log( 2.0*sqrt(PI)*(c_set-(J0+Jpar)*B0*Bbar/(MU0*Lbar)/(N_tmp*Nbar*density*ee*Zi))/vth_et );
         phi_sh -= phi_sh0;
         phi_sh /= Bbar*Va*Lbar;*/
        } else {
            Jpar_sh = Ne0 * Nbar * density * ee;
            // Jpar_sh *= c_set -  vth_et/(2.0*sqrt(PI)) * exp( - ee*((phi+phi0)*Va*Lbar*Bbar)/(KB*Te_tmp*Tebar*eV_K) );
            Jpar_sh *= vth_et / (2.0 * sqrt(PI)) * (ee * ((phi)*Va * Lbar * Bbar) / (KB * Te0 * Tebar * eV_K));
            // Jpar_sh -= Jpar_sh0;
            Jpar_sh *= MU0 * Lbar / (Bbar);
            SBC_Dirichlet(Jpar, Jpar_sh, PF_limit, PF_limit_range);
        }

        if (diffusion_par > 0.0 || kappa_par_i_const > 0.0) {
            q_si = -2. / 3. * gamma_i_BC * Pi * c_se / kappa_par_i; // * (Nbar*density*KB);
        } else {
            q_si = 0.;
        }
        if (diffusion_par > 0.0 || kappa_par_e_const > 0.0) {
            q_se = -2. / 3. * gamma_e_BC * Pe * c_se / kappa_par_e; // * (Nbar*density*KB);
        } else {
            q_se = 0.;
        }

        if (compress0) {
            if (full_sbc)
                SBC_Dirichlet(Vipar, c_set/Va, PF_limit, PF_limit_range);
            else
                SBC_Dirichlet(Vipar, c_se, PF_limit, PF_limit_range);
        }
        // SBC_Gradpar(U, 0.0, PF_limit, PF_limit_range);
        SBC_Gradpar(Ni, 0.0, PF_limit, PF_limit_range);
        if (!Landau) {
            //SBC_Gradpar(Ni, 0.0, PF_limit, PF_limit_range);
            SBC_Gradpar(Ti, q_si, PF_limit, PF_limit_range);
            SBC_Gradpar(Te, q_se, PF_limit, PF_limit_range);
        } else {
            SBC_Gradpar(Ti, 0., PF_limit, PF_limit_range);
            SBC_Gradpar(Te, 0., PF_limit, PF_limit_range);
            //SBC_Gradpar(Vipar, 0., PF_limit, PF_limit_range);
        }
    }

    // update Landau parallel heat flux
    if (diffusion_par && Landau) {
        if (full_sbc) {
            SBC_value_i = gamma_i_BC * (Pi + Pi0) * c_set / (1.6 * N0  * vth_i);
            SBC_value_e = gamma_e_BC * (Pe + Pe0) * c_set / (1.6 * Ne0 * vth_e);
        } else {
            SBC_value_i = gamma_i_BC * Pi * c_se * Va / (1.6 * N0  * vth_i);
            SBC_value_e = gamma_e_BC * Pe * c_se * Va / (1.6 * Ne0 * vth_e);
        }
        if (Landau_coll) {
            q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va * iSign_kpar_wcoll(Ti, kappa_i, 2, nLorentzian, SBC_value_i);
            q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar_wcoll(Te, kappa_e, 2, nLorentzian, SBC_value_e);
        } else {
            q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va * iSign_kpar(Ti, kappa_0, 2, nLorentzian, SBC_value_i);
            q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar(Te, kappa_0, 2, nLorentzian, SBC_value_e);
        }
        mesh->communicate(q_par_i);
        mesh->communicate(q_par_e);
        q_par_i.applyBoundary();
        q_par_e.applyBoundary();

/*    SBC_value_i = gamma_i_BC * (Pi + Pi0) * c_set / Va;
    SBC_value_e = gamma_e_BC * (Pe + Pe0) * c_set / Va;
    SBC_Dirichlet(q_par_i, SBC_value_i, PF_limit, PF_limit_range);
    SBC_Dirichlet(q_par_e, SBC_value_e, PF_limit, PF_limit_range);
*/
        SBC_Gradpar(q_par_i, 0., PF_limit, PF_limit_range);
        SBC_Gradpar(q_par_e, 0., PF_limit, PF_limit_range);
    }

    if (neutral) {
        N_tmp = log(Rcyc_Nn*(N0+Ni));
        SBC_Dirichlet(lNn, N_tmp, PF_limit, PF_limit_range);
        SBC_Dirichlet(Vn, -Rcyc_Vn*c_set/Va, PF_limit, PF_limit_range);
    }

    if (jpar_bndry_width > 0) {
        // Zero j in boundary regions. Prevents vorticity drive
        // at the boundary
        for (jx = 0; jx < jpar_bndry_width; jx++)
            for (jy = 0; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz - 1; jz++) {
                    if (mesh->firstX())
                        Jpar[jx][jy][jz] = 0.0;
                    if (mesh->lastX())
                        Jpar[mesh->ngx - 1 - jx][jy][jz] = 0.0;
                }
    }
    // Smooth j in x
    if (smooth_j_x && (gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
        Jpar = smooth_x(Jpar);
        /*Jpar = smooth_x(Jpar);
    Jpar = smooth_x(Jpar);
    Jpar = smooth_x(Jpar);
    Jpar = smooth_x(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);
    Jpar = smooth_y(Jpar);*/
        mesh->communicate(Jpar);
        Jpar.applyBoundary();
    }

    if (mask_j_x) {
        Jpar *= mask_jx1d;
    }

    if (compresse) {
        if (compress0) {
            if (nonlinear)
                // Ni*Zi + Nimp*Zimp = Ne
                // if (impurity_prof)
                // Vepar = Vipar*N_tmp*Zi/Ne_tmp + Vipar*N_imp0*Z_imp/Ne_tmp - Jpar / Ne_tmp * Vepara;
                // else
                Vepar = Vipar - Jpar / Ne_tmp * Vepara;
            else
                // if (impurity_prof)
                // Vepar = Vipar*N0*Zi/Ne0 + Vipar*N_imp0*Z_imp/Ne0 - Jpar / Ne0 * Vepara;
                // else
                Vepar = Vipar - Jpar / Ne0 * Vepara;
            Vepar.applyBoundary();
            mesh->communicate(Vepar);
        } else {
            if (nonlinear)
                Vepar = -Jpar / Ne_tmp * Vepara;
            else
                Vepar = -Jpar / Ne0 * Vepara;
        }
    }

    // xqx begin
    // Get Delp2(J) from J
    Jpar2 = -Delp2(Jpar);

    Jpar2.applyBoundary();
    mesh->communicate(Jpar2);
    if (jpar_bndry_width > 0) {
        // Zero jpar2 in boundary regions. Prevents vorticity drive
        // at the boundary
        for (jx = 0; jx < jpar_bndry_width; jx++)
            for (jy = 0; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz - 1; jz++) {
                    if (mesh->firstX())
                        Jpar2[jx][jy][jz] = 0.0;
                    if (mesh->lastX())
                        Jpar2[mesh->ngx - 1 - jx][jy][jz] = 0.0;
                }
    }
#if DEBUG_6F>0
    output.write("I see you 3!\n");//xia
#endif

    ////////////////////////////////////////////////////
    // Parallel electric field

    if (!emass) {
        if (evolve_psi) {
            ddt(Psi) = 0.0;
            ddt(Psi) = -Grad_parP(phi, CELL_CENTRE) / B0 - eta * Jpar / B0;

            if (diamag && diamag_phi0) {
                if (diamag_er)
                    ddt(Psi) -= V_dot_Grad(Ve0, Psi);
                else
                    ddt(Psi) -= bracket(phi0, Psi, bm_exb); // Equilibrium flow
            }

            /*
      if(experiment_Er){
        ddt(Psi) -= V_dot_Grad(Ve0_net, Psi);
      }
      */

            if (thermal_force) {
                // grad_par(T_e) correction
                ddt(Psi) += 0.71 * Psipara1 * Grad_parP(Te, CELL_YLOW) / B0;
                ddt(Psi) -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
            }

            if (eHall) {
                ddt(Psi) += Psipara1 * Grad_parP(Pe, CELL_YLOW) / B0 / Ne0;
                ddt(Psi) -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
            }

            if (hyperdiff_par_apar4 > 0.0) {
                tmpA2 = Grad2_par2new(Psi);
                mesh->communicate(tmpA2);
                tmpA2.applyBoundary();
                ddt(Psi) -= hyperdiff_par_apar4 * Grad2_par2new(tmpA2);
            }

            // Hyper-resistivity
            if (hyperresist > 0.0) {
                ddt(Psi) += hyperresist * Delp2(Jpar / B0);
            }
            if (sink_Psir > 0.0) {
                ddt(Psi) -= sink_Psir * sink_tanhxr(P0, Psi, spsi_widthr, spsi_lengthr, false);
            }

        } else { // evolve_psi
            ddt(Apar) = 0.0;
            ddt(Apar) = -Grad_parP(phi, CELL_CENTRE) - eta * Jpar;

            if (diamag && diamag_phi0) {
                if (diamag_er)
                    ddt(Apar) -= V_dot_Grad(Ve0, Apar);
                else
                    ddt(Apar) -= bracket(phi0, Apar, bm_exb);
            }

            if (thermal_force) {
                // grad_par(T_e) correction
                ddt(Apar) += 0.71 * Psipara1 * Grad_parP(Te, CELL_YLOW);
                ddt(Apar) -= 0.71 * Psipara1 * bracket(Apar, Te0, bm_mag);
            }

            if (eHall) {
                ddt(Apar) += Psipara1 * Grad_parP(Pe, CELL_YLOW) / Ne0;
                ddt(Apar) -= Psipara1 * bracket(Apar, Pe0, bm_mag) / Ne0;
            }

            if (hyperdiff_par_apar4 > 0.0) {
                tmpA2 = Grad2_par2new(Apar);
                mesh->communicate(tmpA2);
                tmpA2.applyBoundary();
                ddt(Apar) -= hyperdiff_par_apar4 * Grad2_par2new(tmpA2);
            }

            // Hyper-resistivity
            if (hyperresist > 0.0) {
                //ddt(Psi) += hyperresist * Delp2(Jpar / B0);
                ddt(Apar) += hyperresist * Delp2(Jpar);
            }

            if (sink_Psir > 0.0) {
                ddt(Apar) -= sink_Psir * sink_tanhxr(P0, Apar, spsi_widthr, spsi_lengthr, false);
            }
        }
    } else { // emass
        ddt(Ajpar) = 0.0;
        ddt(Ajpar) = -Grad_parP(phi, CELL_CENTRE) / B0 - eta * Jpar / B0;

        if (diamag && diamag_phi0) {
            if (diamag_er)
                ddt(Ajpar) -= V_dot_Grad(Ve0, Psi);
            else
                ddt(Ajpar) -= bracket(phi0, Psi, bm_exb); // Equilibrium flow
        }

        if (thermal_force) {
            // grad_par(T_e) correction
            ddt(Ajpar) += 0.71 * Psipara1 * Grad_parP(Te, CELL_YLOW) / B0;
            ddt(Ajpar) -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
        }

        if (eHall) {
            ddt(Ajpar) += Psipara1 * Grad_parP(Pe, CELL_YLOW) / B0 / Ne0;
            ddt(Ajpar) -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
        }

        // Hyper-resistivity
        if (hyperresist > 0.0) {
            ddt(Ajpar) += hyperresist * Delp2(Jpar / B0);
        }
    }

    if (output_ohm) {
        ohm_phi = -Grad_parP(phi, CELL_CENTRE) / B0 - bracket(phi0, Psi, bm_exb);
        mesh->communicate(ohm_phi);
        ohm_phi.applyBoundary();
        ohm_hall = Psipara1 * (Grad_parP(Pe, CELL_YLOW) / (B0 * Ne0) - bracket(Psi, Pe0, bm_mag) / Ne0);
        mesh->communicate(ohm_hall);
        ohm_hall.applyBoundary();
        ohm_thermal = 0.71 * Psipara1 * (Grad_parP(Te, CELL_YLOW) / B0 - bracket(Psi, Te0, bm_mag));
        mesh->communicate(ohm_thermal);
        ohm_thermal.applyBoundary();
    }
#if DEBUG_6F>0
    output.write("I see you 4!\n");//xia
#endif

    ////////////////////////////////////////////////////
    // Vorticity equation

    ddt(U) = 0.0;

    if (BScurrent) {
        if (evolve_psi)
            ddt(U) = -(B0 ^ 2) * bracket(Psi, Jpar_BS0 / B0, bm_mag) * B0;
        else
            ddt(U) = -(B0 ^ 2) * bracket(Apar, Jpar_BS0 / B0, bm_mag);
    } else {
        if (evolve_psi)
            ddt(U) = -(B0 ^ 2) * bracket(Psi, J0 / B0, bm_mag) * B0; // Grad j term
        else
            ddt(U) = -(B0 ^ 2) * bracket(Apar, J0 / B0, bm_mag); // Grad j term
    }

    ddt(U) += 2.0 * Upara1 * b0xcv * Grad(P); // curvature term

    ddt(U) += (B0 ^ 2) * Grad_parP(Jpar / B0, CELL_CENTRE); // b dot grad j

    if (diamag && diamag_phi0) {
        if (diamag_er)
            ddt(U) -= V_dot_Grad(Ve0, U);
        else
            ddt(U) -= bracket(phi0, U, bm_exb); // Equilibrium flow
    }

    if (experiment_Er && KH_term)
        ddt(U) -= bracket(phi, U0_net, bm_exb);

    if (compress0 && include_vipar)
        ddt(U) -= Vipar * Grad_parP(U0_net, CELL_CENTRE);

    if (nonlinear) {
        ddt(U) -= bracket(phi, U, bm_exb); // Advection
        /*if (compress0)
      // ddt(U) -= Vipar*Grad_par(U);
      ddt(U) -= Vpar_Grad_par(Vipar, U);*/
        if (compress0 && include_vipar)
            ddt(U) -= Vpar_Grad_par(Vipar, U);
    }

    // xqx: parallel hyper-viscous diffusion for vector potential
    if (hyperdiff_par_u4 > 0.0) {
        tmpU2 = Grad2_par2new(U);
        mesh->communicate(tmpU2);
        tmpU2.applyBoundary();
        ddt(U) -= hyperdiff_par_u4 * Grad2_par2new(tmpU2);
    }

    if (hyperdiff_perp_u4 > 0.0) {
        tmpU2 = Delp2(U);
        mesh->communicate(tmpU2);
        tmpU2.applyBoundary();
        ddt(U) -= hyperdiff_perp_u4 * Delp2(tmpU2);
    }

    if (parallel_viscous && compress0) {
        //ddt(U) += 0.333333 * Upara1 * b0xcv * Grad(pi_ci);
        ddt(U) -= 0.666667 * Upara1 * b0xcv * Grad(pi_ci);
    }

    if (gyroviscous) {
        if (diamag_er && diamag_phi0)
            // Dperp2Phi0 = -Div( B0vec^Ve0 ) - Div(B0vec)*(B0vec*Er0*B0);
            Dperp2Phi0 = Div(Ve0 ^ B0vec);
        else
            Dperp2Phi0 = Field3D(Delp2(phi0));
        Dperp2Phi0.applyBoundary();
        mesh->communicate(Dperp2Phi0);

        Dperp2Phi = Delp2(phi);
        Dperp2Phi.applyBoundary();
        mesh->communicate(Dperp2Phi);

        /*
    GradPhi02 = Field3D( Grad(B0*phi0)*Grad(B0*phi0) / (B0*B0) );
    GradPhi02.applyBoundary();
    mesh->communicate(GradPhi02);

    GradparPhi02 = Field3D( Grad_par(B0*phi0)*Grad_par(B0*phi0) / (B0*B0) );
    GradparPhi02.applyBoundary();
    mesh->communicate(GradparPhi02);

    GradcPhi = Grad(B0*phi0)*Grad(B0*phi) / (B0*B0);
    GradcPhi.applyBoundary();
    mesh->communicate(GradcPhi);

    GradcparPhi = Grad_par(B0*phi0)*Grad_par(B0*phi) / (B0*B0);
    GradcparPhi.applyBoundary();
    mesh->communicate(GradcparPhi);*/

        if (diamag_er && diamag_phi0)
            // GradPhi02 = Er0*Er0; // test
            GradPhi02 = (Ve0 ^ B0vec) * (Ve0 ^ B0vec) / (B0 * B0);
        else
            GradPhi02 = Field3D(Grad_perp(phi0) * Grad_perp(phi0) / (B0 * B0));
        GradPhi02.applyBoundary();
        mesh->communicate(GradPhi02);

        if (diamag_er && diamag_phi0)
            // GradcPhi = -Er0*Grad(phi); //test
            GradcPhi = (Ve0 ^ B0vec) * Grad_perp(phi) / (B0 * B0);
        else
            GradcPhi = Grad_perp(phi0) * Grad_perp(phi) / (B0 * B0);
        GradcPhi.applyBoundary();
        mesh->communicate(GradcPhi);

        Dperp2Pi0 = Field3D(Delp2(Pi0));
        Dperp2Pi0.applyBoundary();
        mesh->communicate(Dperp2Pi0);

        Dperp2Pi = Delp2(Pi);
        Dperp2Pi.applyBoundary();
        mesh->communicate(Dperp2Pi);

        if (diamag_er && diamag_phi0)
            bracketPhi0P = V_dot_Grad(Ve0, Pi);
        else
            bracketPhi0P = bracket(phi0, Pi, bm_exb);
        bracketPhi0P.applyBoundary();
        mesh->communicate(bracketPhi0P);

        bracketPhiP0 = bracket(phi, Pi0, bm_exb);
        bracketPhiP0.applyBoundary();
        mesh->communicate(bracketPhiP0);

        ddt(U) -= 0.5 * Upara2 * bracket(Pi, Dperp2Phi0, bm_exb) / B0;
        ddt(U) -= 0.5 * Upara2 * bracket(Pi0, Dperp2Phi, bm_exb) / B0;
        ddt(U) += Upara3 * B0 * bracket(N0, GradcPhi, bm_exb);
        ddt(U) += 0.5 * Upara3 * B0 * bracket(Ni, GradPhi02, bm_exb);
        ddt(U) += 0.5 * Upara2 * bracket(phi, Dperp2Pi0, bm_exb) / B0;
        if (diamag_er && diamag_phi0)
            ddt(U) += 0.5 * Upara2 * V_dot_Grad(Ve0, Dperp2Pi) / B0;
        else
            ddt(U) += 0.5 * Upara2 * bracket(phi0, Dperp2Pi, bm_exb) / B0;
        ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhi0P) / B0;
        ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhiP0) / B0;

        if (impurity_prof && impurity_gyro) {
            Dperp2Pimp0 = Field3D(Delp2(P_imp0));
            Dperp2Pimp0.applyBoundary();
            mesh->communicate(Dperp2Pimp0);

            bracketPhiPimp0 = bracket(phi, P_imp0, bm_exb);
            bracketPhiPimp0.applyBoundary();
            mesh->communicate(bracketPhiPimp0);

            ddt(U) -= 0.5 * Upara_imp * Upara2 * bracket(P_imp0, Dperp2Phi, bm_exb) / B0;
            ddt(U) += Upara_imp * Upara3 * B0 * bracket(N_imp0, GradcPhi, bm_exb);
            ddt(U) += 0.5 * Upara_imp * Upara2 * bracket(phi, Dperp2Pimp0, bm_exb) / B0;
            ddt(U) -= 0.5 * Upara_imp * Upara2 * Delp2(bracketPhiPimp0) / B0;
        }

        if (nonlinear) {
            /*GradPhi2 = Grad(B0*phi)*Grad(B0*phi) / (B0*B0);
      GradPhi2.applyBoundary();
      mesh->communicate(GradPhi2);

      GradparPhi2 = Grad_par(B0*phi)*Grad_par(B0*phi) / (B0*B0);
      GradparPhi2.applyBoundary();
      mesh->communicate(GradparPhi2);*/

            GradPhi2 = Grad_perp(phi) * Grad_perp(phi) / (B0 * B0);
            GradPhi2.applyBoundary();
            mesh->communicate(GradPhi2);

            bracketPhiP = bracket(phi, Pi, bm_exb);
            bracketPhiP.applyBoundary();
            mesh->communicate(bracketPhiP);

            ddt(U) -= 0.5 * Upara2 * bracket(Pi, Dperp2Phi, bm_exb) / B0;
            ddt(U) += 0.5 * Upara3 * B0 * bracket(N0, GradPhi2, bm_exb);
            ddt(U) += Upara3 * B0 * bracket(Ni, GradcPhi, bm_exb);
            ddt(U) += 0.5 * Upara2 * bracket(phi, Dperp2Pi, bm_exb) / B0;
            ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhiP) / B0;

            if (impurity_prof && impurity_gyro) {
                ddt(U) += 0.5 * Upara_imp * Upara3 * B0 * bracket(N_imp0, GradPhi2, bm_exb);
            }
        }
    }

    // TODO: check
    if (output_transfer) {
        T_R = -bracket(phi, U - Upara0 / B0 * Dperp2Pi, bm_exb);
        mesh->communicate(T_R);
        T_R.applyBoundary();
        T_M = (B0 ^ 2) * Grad_parP(Jpar / B0, CELL_CENTRE);
        mesh->communicate(T_M);
        T_M.applyBoundary();
        T_ID = -bracket(phi, Upara0 / B0 * Dperp2Pi, bm_exb);
        mesh->communicate(T_ID);
        T_ID.applyBoundary();
        T_C = 2.0 * Upara1 * b0xcv * Grad(P);
        mesh->communicate(T_C);
        T_C.applyBoundary();
        T_G = -0.5 * Upara2 * bracket(Pi, Dperp2Phi, bm_exb) / B0;
        T_G += 0.5 * Upara3 * B0 * bracket(N0, GradPhi2, bm_exb);
        T_G += Upara3 * B0 * bracket(Ni, GradcPhi, bm_exb);
        T_G += 0.5 * Upara2 * bracket(phi, Dperp2Pi, bm_exb) / B0;
        T_G -= 0.5 * Upara2 * Delp2(bracketPhiP) / B0;
        T_G.applyBoundary();
        mesh->communicate(T_G);
    }

    // Viscosity terms
    if (viscos_par > 0.0)
        ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity

    if (hyperviscos > 0.0) {
        // Calculate coefficient.

        hyper_mu_x = hyperviscos * mesh->g_11 * SQ(mesh->dx) * abs(mesh->g11 * D2DX2(U)) / (abs(U) + 1e-3);
        hyper_mu_x.applyBoundary("dirichlet"); // Set to zero on all boundaries

        ddt(U) += hyper_mu_x * mesh->g11 * D2DX2(U);

        if (first_run) { // Print out maximum values of viscosity used on this processor
            output.write("   Hyper-viscosity values:\n");
            output.write("      Max mu_x = %e, Max_DC mu_x = %e\n", max(hyper_mu_x), max(hyper_mu_x.DC()));
        }
    }

    // left edge sink terms
    if (sink_Ul > 0.0) {
        ddt(U) -= sink_Ul * sink_tanhxl(P0, U, su_widthl, su_lengthl, false); // core sink
    }

    // right edge sink terms
    if (sink_Ur > 0.0) {
        ddt(U) -= sink_Ur * sink_tanhxr(P0, U, su_widthr, su_lengthr, false); // sol sink
    }
#if DEBUG_6F>0
    output.write("I see you 5!\n");//xia
#endif

    ///////////////////////////////////////////////
    // number density equation

    ddt(Ni) = 0.0;

    if (NiAmp>0) ddt(Ni) += NiSource;

    ddt(Ni) -= bracket(phi, N0, bm_exb);

    if (continuity) {
        ddt(Ni) -= 2.0 * N0 / B0 * b0xcv * Grad(phi);
        // ddt(Ni) -= 2.0 * Ni/B0 *  b0xcv*Grad(phi0);
        // ddt(Ni) -= 2.0 * Nipara1 * b0xcv*Grad(Pi0) / N0/B0;    // balanced by above
        ddt(Ni) -= 2.0 * Ni / B0 * b0xcv * (Ve0_net ^ B0vec); // net flow
        if (diamag) {
            ddt(Ni) -= 2.0 * Nipara1 * b0xcv * Grad(Pi) / B0;
        }
        if (nonlinear) {
            ddt(Ni) -= 2.0 * Ni / B0 * b0xcv * Grad(phi);
        }
    }

    if (diamag && diamag_phi0) {
        if (diamag_er)
            ddt(Ni) -= V_dot_Grad(Ve0, Ni);
        else
            ddt(Ni) -= bracket(phi0, Ni, bm_exb); // Equilibrium flow
    }

    if (nonlinear) {
        ddt(Ni) -= bracket(phi, Ni, bm_exb); // Advection
    }

    if (compress0) {
        // ddt(Ni) -= Vipar * Grad_parP(N0, CELL_YLOW);
        // ddt(Ni) -= Vpar_Grad_par(Vipar, N0);
        if (include_vipar)
            ddt(Ni) -= Vipar * Grad_parP(N0, CELL_CENTRE);

        if (continuity) {
            ddt(Ni) -= N0 * B0 * Grad_parP(Vipar / B0, CELL_CENTRE);
        }
        if (nonlinear) {
            // ddt(Ni) -= Vipar * Grad_par(Ni, CELL_YLOW);

            // ddt(Ni) -= Vpar_Grad_par(Vipar, Ni);
            // ddt(Ni) += Vipar * bracket(Psi, N0, bm_mag)*B0;

            if (include_vipar)
                ddt(Ni) -= Vpar_Grad_par(Vipar, Ni);

            if (continuity)
                ddt(Ni) -= Ni * B0 * Grad_par(Vipar / B0, CELL_CENTRE);
        }
    }

    if (radial_diffusion)
        ddt(Ni) += diff_radial * Delp2(Ni);

    /*  if (neoclassic_i)
      {
        tmpddx2 = D2DX2(Ni.DC());
        mesh->communicate(tmpddx2);
        tmpddx2.applyBoundary();
        ddt(Ni) += mesh->g11*(Dri_neo * tmpddx2).DC();

        partf_neo_i = -Dri_neo * sqrt(mesh->g11)*DDX(Ni);
        mesh->communicate(partf_neo_i);
        partf_neo_i.applyBoundary();
      }
  */
    // M: 4th order Parallel diffusion terms
    if (hyperdiff_par_n4 > 0.0) {
        tmpN2 = Grad2_par2new(Ni);
        mesh->communicate(tmpN2);
        tmpN2.applyBoundary();
        ddt(Ni) -= hyperdiff_par_n4 * Grad2_par2new(tmpN2);
    }

    if (hyperdiff_perp_n4 > 0.0) {
        tmpN2 = Delp2(Ni);
        mesh->communicate(tmpN2);
        tmpN2.applyBoundary();
        ddt(Ni) -= hyperdiff_perp_n4 * Delp2(tmpN2);
    }
#if DEBUG_6F>0
    output.write("I see you 6!\n");//xia
#endif

    ///////////////////////////////////////////////                                // ion temperature equation

    ddt(Ti) = 0.0;

    if (TiAmp>0) ddt(Ti) += TiSource;

    ddt(Ti) -= bracket(phi, Ti0, bm_exb);

    if (continuity) {
        ddt(Ti) -= 4.0 / 3.0 * Ti0 / B0 * b0xcv * Grad(phi);
        // ddt(Ti) -= 4.0/3.0 * Ti * b0xcv*Grad(phi0*B0) / B0;
        // ddt(Ti) -= 4.0/3.0 * Tipara2 * Ti/N0 * b0xcv*Grad(Pi0) / (N0*B0);   // balanced by above
        ddt(Ti) -= 4.0 / 3.0 * Ti / B0 * b0xcv * (Ve0_net ^ B0vec); // net flow
        if (diamag)
            ddt(Ti) -= 4.0 / 3.0 * Tipara2 * Ti0 / N0 * b0xcv * Grad(Pi) / B0;
        if (nonlinear) {
            ddt(Ti) -= 4.0 / 3.0 * Ti / B0 * b0xcv * Grad(phi);
            if (diamag)
                ddt(Ti) -= 4.0 / 3.0 * Tipara2 * Ti / N0 * b0xcv * Grad(Pi) / B0;
        }
    }

    if (energy_flux) {
        // ddt(Ti) -= 10.0/3.0 * Tipara2 * Ti0/B0 * b0xcv*Grad(Ti);
        ddt(Ti) -= 10.0 / 3.0 * Tipara2 / B0 * V_dot_Grad(Ti0 * b0xcv, Ti);
        ddt(Ti) -= 10.0 / 3.0 * Tipara2 * Ti / B0 * b0xcv * Grad(Ti0);
        if (nonlinear)
            // ddt(Ti) -= 10.0/3.0 * Tipara2 * Ti/B0 * b0xcv*Grad(Ti);
            ddt(Ti) -= 10.0 / 3.0 * Tipara2 / B0 * V_dot_Grad(Ti * b0xcv, Ti);
    }

    if (diamag && diamag_phi0) {
        if (diamag_er)
            ddt(Ti) -= V_dot_Grad(Ve0, Ti);
        else
            ddt(Ti) -= bracket(phi0, Ti, bm_exb); // Equilibrium flow
    }

    if (nonlinear) {
        ddt(Ti) -= bracket(phi, Ti, bm_exb); // Advection
    }

    if (compress0) {
        // ddt(Ti) -= Vipar * Grad_parP(Ti0, CELL_YLOW);
        // ddt(Ti) -= Vpar_Grad_par(Vipar, Ti0);
        if (include_vipar)
            ddt(Ti) -= Vipar * Grad_parP(Ti0, CELL_CENTRE);

        if (continuity) {
            ddt(Ti) -= 2.0 / 3.0 * Ti0 * B0 * Grad_parP(Vipar / B0, CELL_CENTRE);
        }
        if (nonlinear) {
            // ddt(Ti) -= Vipar * Grad_par(Ti, CELL_YLOW);

            // ddt(Ti) -= Vpar_Grad_par(Vipar, Ti);
            // ddt(Ti) += Vipar * bracket(Psi, Ti0, bm_mag)*B0;

            if (include_vipar)
                ddt(Ti) -= Vpar_Grad_par(Vipar, Ti);

            if (continuity)
                ddt(Ti) -= 2.0 / 3.0 * Ti * B0 * Grad_par(Vipar / B0, CELL_CENTRE);
        }
    }

    if (energy_exch) {
        ddt(Ti) += 2.0 * Zi * Tbar * nu_e / (ratio_pe * AA) * (Te - Ti);
    }

    if (diffusion_par > 0.0 || kappa_par_i_const > 0.0) {
        if (Landau && diffusion_par > 0.0) {
            if (diff_par_flutter)
                ddt(Ti) -= Grad_parP(q_par_i, CELL_CENTRE) / N0;
            else
                ddt(Ti) -= Grad_par(q_par_i, CELL_CENTRE) / N0;
        } else {
            ddt(Ti) += kappa_par_i * Grad2_par2(Ti) / N0; // Parallel diffusion
            ddt(Ti) += Grad_par(kappa_par_i, CELL_CENTRE) * Grad_par(Ti, CELL_YLOW) / N0;

            if (diff_par_flutter) {
                if (nonlinear) {
                    if (evolve_psi)
                        bracket1i = -bracket(Psi, Ti + Ti0, bm_mag) * B0;
                    else
                        bracket1i = -bracket(Apar, Ti + Ti0, bm_mag);
                } else {
                    if (evolve_psi)
                        bracket1i = -bracket(Psi, Ti0, bm_mag) * B0;
                    else
                        bracket1i = -bracket(Apar, Ti0, bm_mag);
                }
                mesh->communicate(bracket1i);
                bracket1i.applyBoundary();
                gradpar_ti = Grad_par(Ti, CELL_YLOW);
                mesh->communicate(gradpar_ti);
                gradpar_ti.applyBoundary();

                ddt(Ti) += Grad_par(kappa_par_i, CELL_CENTRE) * bracket1i / N0;
                ddt(Ti) += kappa_par_i * Grad_par(bracket1i, CELL_YLOW) / N0;
                if (nonlinear) {
                    if (evolve_psi) {
                        ddt(Ti) -= bracket(Psi, kappa_par_i, bm_mag) * B0 * gradpar_ti / N0;
                        ddt(Ti) -= kappa_par_i * bracket(Psi, gradpar_ti, bm_mag) * B0 / N0;
                        ddt(Ti) -= bracket(Psi, kappa_par_i, bm_mag) * B0 * bracket1i / N0;
                        ddt(Ti) -= kappa_par_i * bracket(Psi, bracket1i, bm_mag) * B0 / N0;
                    } else {
                        ddt(Ti) -= bracket(Apar, kappa_par_i, bm_mag) * gradpar_ti / N0;
                        ddt(Ti) -= kappa_par_i * bracket(Apar, gradpar_ti, bm_mag) / N0;
                        ddt(Ti) -= bracket(Apar, kappa_par_i, bm_mag) * bracket1i / N0;
                        ddt(Ti) -= kappa_par_i * bracket(Apar, bracket1i, bm_mag) / N0;
                    }
                }

                if (output_flux_par) {
                    heatf_par_flutter_i = -kappa_par_i * bracket1i;
                }
            }
        }
    }

    if (diffusion_perp > 0.0) {
        ddt(Ti) += kappa_perp_i * Delp2(Ti) / N0; // Parallel diffusion
        ddt(Ti) += Grad_perp(kappa_perp_i) * Grad_perp(Ti) / N0;
    }

    if (gyroviscous && compress0 && nonlinear) {
        ddt(Ti) -= 1.333333 * Tipara3 * (Ti0 + Ti) * Vipar * b0xcv * Grad(Vipar) / B0;
    }

    if (parallel_viscous) {
        if (compress0)
            ddt(Ti) -= 0.444444 * pi_ci / Tau_ie / N0 / sqrt(B0) * Grad_parP(Vipar * sqrt(B0), CELL_CENTRE);
        if (nonlinear)
            ddt(Ti) += 0.222222 * pi_ci / Tau_ie / N0 / N0 / B0 * bracket(N0+Ni, Ti0+Ti, bm_exb);
    }

    if (neoclassic_i) {
        tmpddx2 = D2DX2(Ti.DC());
        mesh->communicate(tmpddx2);
        tmpddx2.applyBoundary();
        ddt(Ti) += mesh->g11 * (xii_neo * tmpddx2.DC()).DC();

        heatf_neo_i = -xii_neo * sqrt(mesh->g11) * DDX(Ti);
        mesh->communicate(heatf_neo_i);
        heatf_neo_i.applyBoundary();
    }

    // M: 4th order Parallel diffusion terms
    if (hyperdiff_par_ti4 > 0.0) {
        tmpTi2 = Grad2_par2new(Ti);
        mesh->communicate(tmpTi2);
        tmpTi2.applyBoundary();
        ddt(Ti) -= hyperdiff_par_ti4 * Grad2_par2new(tmpTi2);
    }

    if (hyperdiff_perp_ti4 > 0.0) {
        tmpTi2 = Delp2(Ti);
        mesh->communicate(tmpTi2);
        tmpTi2.applyBoundary();
        ddt(Ti) -= hyperdiff_perp_ti4 * Delp2(tmpTi2);
    }
#if DEBUG_6F>0
    output.write("I see you 7!\n");//xia
#endif

    ///////////////////////////////////////////////                                // electron temperature equation
    ddt(Te) = 0.0;

    if (TeAmp>0) ddt(Te) += TeSource;

    ddt(Te) -= bracket(phi, Te0, bm_exb);

    if (continuity) {
        ddt(Te) -= 4.0 / 3.0 * Te0 / B0 * b0xcv * Grad(phi);
        ddt(Te) -= 4.0 / 3.0 * Te / B0 * b0xcv * (Ve0 ^ B0vec);
        ddt(Te) += 4.0 / 3.0 * Tepara2 * Te / Ne0 * b0xcv * Grad(Pe0) / B0;
        if (diamag)
            ddt(Te) += 4.0 / 3.0 * Tepara2 * Te0 / Ne0 * b0xcv * Grad(Pe) / B0;
        if (nonlinear) {
            ddt(Te) -= 4.0 / 3.0 * Te / B0 * b0xcv * Grad(phi);
            if (diamag)
                ddt(Te) += 4.0 / 3.0 * Tepara2 * Te / Ne0 * b0xcv * Grad(Pe) / B0;
        }
    }

    if (energy_flux) {
        // ddt(Te) += 10.0/3.0 * Tepara2 * Te0/B0 * b0xcv*Grad(Te);
        ddt(Te) -= 10.0 / 3.0 * Tepara2 / B0 * V_dot_Grad(-Te0 * b0xcv, Te);
        ddt(Te) += 10.0 / 3.0 * Tepara2 * Te / B0 * b0xcv * Grad(Te0);
        if (thermal_force) {
            ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * Grad_parP(Jpar / B0, CELL_CENTRE);
            if (BScurrent) {
                if (evolve_psi)
                    ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Psi, Jpar_BS0 / B0, bm_mag) * B0;
                else
                    ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Apar, Jpar_BS0 / B0, bm_mag);
                ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te * B0 / Ne0 * Grad_parP(Jpar_BS0 / B0, CELL_CENTRE);
            } else {
                if (evolve_psi)
                    ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Psi, J0 / B0, bm_mag) * B0;
                else
                    ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Apar, J0 / B0, bm_mag);
                ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te * B0 / Ne0 * Grad_parP(J0 / B0, CELL_CENTRE);
            }
        }
        if (nonlinear) {
            // ddt(Te) += 10.0/3.0 * Tepara2 * Te/B0 * b0xcv*Grad(Te);
            ddt(Te) -= 10.0 / 3.0 * Tepara2 / B0 * V_dot_Grad(-Te * b0xcv, Te);
            if (thermal_force)
                ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te * B0 / Ne0 * Grad_par(Jpar / B0, CELL_CENTRE);
        }
    }

    if (diamag && diamag_phi0) {
        if (diamag_er)
            ddt(Te) -= V_dot_Grad(Ve0, Te);
        else
            ddt(Te) -= bracket(phi0, Te, bm_exb); // Equilibrium flow
    }

    if (nonlinear) {
        ddt(Te) -= bracket(phi, Te, bm_exb); // Advection
    }

    if (compresse) {
        // ddt(Te) -= Vepar * Grad_parP(Te0, CELL_YLOW);
        // ddt(Te) -= Vpar_Grad_par(Vepar, Te0);

        if (include_vipar)
            ddt(Te) -= Vepar * Grad_parP(Te0, CELL_CENTRE);

        if (continuity) {
            ddt(Te) -= 2.0 / 3.0 * Te0 * B0 * Grad_parP(Vepar / B0, CELL_CENTRE);
        }
        if (nonlinear) {
            // ddt(Te) -= Vepar * Grad_par(Te, CELL_YLOW);

            // ddt(Te) -= Vpar_Grad_par(Vepar, Te);
            // ddt(Te) += Vepar * bracket(Psi, Te0, bm_mag)*B0;

            if (include_vipar)
                ddt(Te) -= Vpar_Grad_par(Vepar, Te);

            if (continuity)
                ddt(Te) -= 2.0 / 3.0 * Te * B0 * Grad_par(Vepar / B0, CELL_CENTRE);
        }
    }

    if (energy_exch) {
        ddt(Te) -= 2.0 * Tbar * nu_e / (ratio_pe * AA) * (Te - Ti);
        if (BScurrent)
            ddt(Te) += 4.0 / 3.0 * Tepara4 * eta * Jpar_BS0 * Jpar / Ne0;
        else
            ddt(Te) += 4.0 / 3.0 * Tepara4 * eta * J0 * Jpar / Ne0;

        if (nonlinear)
            ddt(Te) += 2.0 / 3.0 * Tepara4 * eta * Jpar * Jpar / Ne0;
    }

    if (diffusion_par > 0.0 || kappa_par_e_const > 0.0) {
        if (Landau && diffusion_par > 0.0) {
            if (diff_par_flutter)
                ddt(Te) -= Grad_parP(q_par_e, CELL_CENTRE) / Ne0;
            else
                ddt(Te) -= Grad_par(q_par_e, CELL_CENTRE) / Ne0;

            if (output_qparcompare) {
                q_par_fl = kappa_par_e * Grad_par(Te);
                q_par_landau = q_par_e;
            }
        } else {
            ddt(Te) += kappa_par_e * Grad2_par2(Te) / Ne0; // Parallel diffusion
            ddt(Te) += Grad_par(kappa_par_e, CELL_CENTRE) * Grad_par(Te, CELL_YLOW) / Ne0;

            if (diff_par_flutter) {
                if (nonlinear) {
                    if (evolve_psi)
                        bracket1e = -bracket(Psi, Te + Te0, bm_mag) * B0;
                    else
                        bracket1e = -bracket(Apar, Te + Te0, bm_mag);
                } else {
                    if (evolve_psi)
                        bracket1e = -bracket(Psi, Te0, bm_mag) * B0;
                    else
                        bracket1e = -bracket(Apar, Te0, bm_mag);
                }
                mesh->communicate(bracket1e);
                bracket1e.applyBoundary();
                gradpar_te = Grad_par(Te, CELL_YLOW);
                mesh->communicate(gradpar_te);
                gradpar_te.applyBoundary();

                ddt(Te) += Grad_par(kappa_par_e, CELL_CENTRE) * bracket1e / Ne0;
                ddt(Te) += kappa_par_e * Grad_par(bracket1e, CELL_YLOW) / Ne0;

                if (nonlinear) {
                    if (evolve_psi) {
                        ddt(Te) -= bracket(Psi, kappa_par_e, bm_mag) * B0 * gradpar_te / Ne0;
                        ddt(Te) -= kappa_par_e * bracket(Psi, gradpar_te, bm_mag) * B0 / Ne0;
                        ddt(Te) -= bracket(Psi, kappa_par_e, bm_mag) * B0 * bracket1e / Ne0;
                        ddt(Te) -= kappa_par_e * bracket(Psi, bracket1e, bm_mag) * B0 / Ne0;
                    } else {
                        ddt(Te) -= bracket(Apar, kappa_par_e, bm_mag) * gradpar_te / Ne0;
                        ddt(Te) -= kappa_par_e * bracket(Apar, gradpar_te, bm_mag) / Ne0;
                        ddt(Te) -= bracket(Apar, kappa_par_e, bm_mag) * bracket1e / Ne0;
                        ddt(Te) -= kappa_par_e * bracket(Apar, bracket1e, bm_mag) / Ne0;
                    }
                }

                if (output_flux_par) {
                    heatf_par_flutter_e = -kappa_par_e * bracket1e;
                }
            }
        }
    }

    if (diffusion_perp > 0.0) {
        ddt(Te) += kappa_perp_e * Delp2(Te) / Ne0; // Parallel diffusion
        ddt(Te) += Grad_perp(kappa_perp_e) * Grad_perp(Te) / Ne0;
    }

    if (neoclassic_e) {
        tmpddx2 = D2DX2(Te.DC());
        mesh->communicate(tmpddx2);
        tmpddx2.applyBoundary();
        ddt(Te) += mesh->g11 * (xie_neo * tmpddx2).DC();

        heatf_neo_e = -xie_neo * sqrt(mesh->g11) * DDX(Te);
        mesh->communicate(heatf_neo_e);
        heatf_neo_e.applyBoundary();
    }

    if (hyperdiff_par_te4 > 0.0) {
        tmpTe2 = Grad2_par2new(Te);
        mesh->communicate(tmpTe2);
        tmpTe2.applyBoundary();
        ddt(Te) -= hyperdiff_par_te4 * Grad2_par2new(tmpTe2);
    }

    if (hyperdiff_perp_te4 > 0.0) {
        tmpTe2 = Delp2(Te);
        mesh->communicate(tmpTe2);
        tmpTe2.applyBoundary();
        ddt(Te) -= hyperdiff_perp_te4 * Delp2(tmpTe2);
    }

    // left edge sink terms
    if (sink_Tel > 0.0) {
        ddt(Te) -= sink_Tel * sink_tanhxl(Te0, Te, ste_widthl, ste_lengthl, false);  //core sink
    }
    // right edge sink terms
    if (sink_Ter > 0.0) {
        ddt(Te) -= sink_Ter * sink_tanhxr(Te0, Te, ste_widthr, ste_lengthr, false); // sol sink
    }

    if (output_flux_par) {
        // gamma_par_i = (N0 + Ni) * Vipar;
        heatf_par_i = -kappa_par_i * Grad_par(Ti, CELL_YLOW);
        mesh->communicate(heatf_par_i);
        heatf_par_i.applyBoundary();
        heatf_par_e = -kappa_par_e * Grad_par(Te, CELL_YLOW);
        mesh->communicate(heatf_par_e);
        heatf_par_e.applyBoundary();
    }

    if (output_vradial) {
        if (diamag_er)
            Vexb = Ve0 + (B0vec ^ Grad(phi)) / (B0 * B0);
        else
            Vexb = (B0vec ^ Grad(phi + phi0)) / (B0 * B0);
        mesh->communicate(Vexb);
        Vexb.applyBoundary();
        if (evolve_psi)
            Vbtilde = -B0vec ^ Grad(Psi);
        else
            Vbtilde = -B0vec ^ Grad(Apar) / B0;
        mesh->communicate(Vbtilde);
        Vbtilde.applyBoundary();
        // Vbti_par = Vipar*Vbtilde.x;
        // Vbte_par = Vepar*Vbtilde.x;
        // mesh->communicate(Vbt_par);
    }
#if DEBUG_6F>0
    output.write("I see you 8!\n");//xia
#endif

    //////////////////////////////////////////////////////////////////////
    if (compress0) // parallel velocity equation
    {
        ddt(Vipar) = 0.0;

        ddt(Vipar) -= Vipara * Grad_parP(P, CELL_YLOW) / N0;
        if (evolve_psi)
            ddt(Vipar) += Vipara * bracket(Psi, P0, bm_mag) * B0 / N0;
        else
            ddt(Vipar) += Vipara * bracket(Apar, P0, bm_mag) / N0;

        if (diamag && diamag_phi0) {
            if (diamag_er)
                ddt(Vipar) -= V_dot_Grad(Ve0, Vipar);
            else
                ddt(Vipar) -= bracket(phi0, Vipar, bm_exb);
        }

        if (nonlinear) {
            ddt(Vipar) -= bracket(phi, Vipar, bm_exb);

            // ddt(Vipar) -= Vipar * Grad_par(Vipar, CELL_CENTRE);
            // ddt(Vipar) -= Vpar_Grad_par(Vipar, Vipar);

            if (include_vipar)
                ddt(Vipar) -= Vpar_Grad_par(Vipar, Vipar);
        }

        // output.write("test point 1.\n");
        if (parallel_viscous && compress0) {
            // Field3D temp_pi;
            // temp_pi = pi_ci / (B0*sqrt(B0));
            ddt(Vipar) -= 0.666667 * Vipara * (B0 * sqrt(B0)) * Grad_parP(pi_ci / (B0 * sqrt(B0)), CELL_YLOW) / N0;
        }

        if (gyroviscous) {
            ddt(Vipar) -= Upara0 * bracket(Pi0, Vipar, bm_exb) / N0;
            if (nonlinear)
                ddt(Vipar) -= Upara0 * bracket(Pi, Vipar, bm_exb) /N0;
        }

        // xqx: parallel hyper-viscous diffusion for vector potential
        if (hyperdiff_par_v4 > 0.0) {
            tmpVp2 = Grad2_par2new(Vipar);
            mesh->communicate(tmpVp2);
            tmpVp2.applyBoundary();
            ddt(Vipar) -= hyperdiff_par_v4 * Grad2_par2new(tmpVp2);
        }

        if (hyperdiff_perp_v4 > 0.0) {
            tmpVp2 = Delp2(Vipar);
            mesh->communicate(tmpVp2);
            tmpVp2.applyBoundary();
            ddt(Vipar) -= hyperdiff_perp_v4 * Delp2(tmpVp2);
        }
        if (sink_vp > 0.0) {
            /// Field2D V0tmp = 0.;
            ddt(Vipar) -= sink_vp * sink_tanhxl(F2D_tmp, Vipar, sp_width, sp_length, false); // sink
        }
    }

    //////////////////////////////////////////////////////////////////////
    if (neutral) // neutral model
    {
        Ti_tmp = Ti0 + Ti;
        N_tmp = N0 + Ni;
        //output.write("Max and min of normalized lNn = %e, %e.\n", max(lNn), min(lNn));
        Nn = exp(lNn);
        //output.write("Max and min of normalized Nn = %e, %e.\n", max(Nn), min(Nn));
        Pn = Nn * Ti_tmp;
//    nu_iz = iz_rate(Ti_tmp*Tibar, AA) * Nbar * density * Lbar / Va;
//    nu_cx = cx_rate(Ti_tmp*Tibar, AA) * Nbar * density * Lbar / Va;
//    nu_rc = rc_rate(Ti_tmp*Tibar, N_tmp*Nbar*density, AA) * Nbar * density * Lbar / Va;
//    sigma_cx = cx_sect(Ti_tmp*Tibar, AA) * Nbar * density * Lbar;
        //output.write("Max and min of normalized nu_iz = %e, %e.\n", max(nu_iz), min(nu_iz));
        //output.write("Max and min of normalized nu_cx = %e, %e.\n", max(nu_cx), min(nu_cx));
        //output.write("Max and min of normalized nu_rc = %e, %e.\n", max(nu_rc), min(nu_rc));
        //output.write("Max and min of normalized si_cx = %e, %e.\n", max(sigma_cx), min(sigma_cx));
        // diffusion and viscosity coefficients
        // zhu: vth_i is caculated in parallel heat flux section
        //vth_i = 9.79e3 * sqrt((Ti0+Ti)*Tibar / AA); // vth_i in m/S.
        Dn = vth_i * vth_i / (Va * Va * N_tmp * nu_cx);
        etan = Nn / (N_tmp + Nn) * vth_i / (Va * sigma_cx);
        // neutral density and parallel velcocity equations
        ddt(lNn) = 0.;
        ddt(Vn) = 0.;
        ddt(lNn) = -Grad_par(Vn) - Vn * Grad_par(lNn);
        ddt(lNn) += Grad_perp(Dn/Ti_tmp) * Grad_perp(Pn) / Nn + Dn / Pn * Delp2(Pn);
        // source/sink terms
        Sn = nu_rc * N_tmp * N_tmp - nu_iz * Nn * N_tmp;
        //output.write("Max and min of normalized Sn = %e, %e.\n", max(Sn), min(Sn));
/*    term1 = -Grad_par(Vn) - Vn * Grad_par(lNn);
    output.write("Max and min of normalized term1 = %e, %e.\n", max(term1), min(term1));
    term2 = Grad_perp(Dn/Ti_tmp) * Grad_perp(Pn) / Nn + Dn / Pn * Delp2(Pn);
    output.write("Max and min of normalized term2 = %e, %e.\n", max(term2), min(term2));
    term3 = Sn/Nn;
    output.write("Max and min of normalized term3 = %e, %e.\n", max(term3), min(term3));
*/
        ddt(lNn) += Sn / Nn + Sn_ext / Nn;
        //ddt(lNn) = 0.; //debug
        ddt(Ni) -= Sn;

        ddt(Vn) = -Vn * Grad_par(Vn) + Dn / Pn * Grad_perp(Pn) * Grad_perp(Vn);
        //ddt(Vn) += Grad_par(etan*Grad_par(Vn)) / Nn;
        //ddt(Vn) += (Grad_par(etan)*Grad_par(Vn)+etan*Grad2_par2new(Vn))/Nn;
        ddt(Vn) -= Upara1 * Tau_ie * Grad_par(Pn);
        Sv = N_tmp * (N_tmp / Nn * nu_rc + nu_cx) * (Vipar - Vn);
/*    term1 = -Vn * Grad_par(Vn);
    term2 = Dn / Pn * Grad_perp(Pn) * Grad_perp(Vn);
    term3 = (Grad_par(etan)*Grad_par(Vn)+etan*Grad2_par2new(Vn))/Nn;
    term4 = -Upara1 * Tau_ie * Grad_par(Pn);
    output.write("Vn Max and min of normalized term1 = %e, %e.\n", max(term1), min(term1));
    output.write("Vn Max and min of normalized term2 = %e, %e.\n", max(term2), min(term2));
    output.write("Vn Max and min of normalized term3 = %e, %e.\n", max(term3), min(term3));
    output.write("Vn Max and min of normalized term4 = %e, %e.\n", max(term4), min(term4));
    output.write("Vn Max and min of normalized term5 = %e, %e.\n", max(Sv), min(Sv));
*/
        //output.write("Max and min of normalized Vn = %e, %e.\n", max(Vn), min(Vn));
        ddt(Vn) += Sv;
        //ddt(Vn) =0.; // debug

        ddt(Vipar) -= Nn * (nu_iz + nu_cx) * (Vipar - Vn);
        ddt(Te) -= nu_iz * Nn * (Te0 + Te + 9.07/Tebar); // neglect recombanation effects
        // strictly speaking, one also need to modify vorticity equs.
    }

    ///////////////////////////////////////////////////////////////////////

    if (PF_limit) {
        // Vipar = PF_filter(Vipar, PF_limit_range);
        // Jpar = PF_filter(Jpar, PF_limit_range);
        if (evolve_psi)
            Psi = PF_filter(Psi, PF_limit_range);
        else
            Apar = PF_filter(Apar, PF_limit_range);
        Ni = PF_filter(Ni, PF_limit_range);
        // Ti = PF_filter(Ti, PF_limit_range);
        // Te = PF_filter(Te, PF_limit_range);
    }

    if (filter_z) {
        // Filter out all except filter_z_mode

        if (!emass) {
            if (evolve_psi)
                ddt(Psi) = filter(ddt(Psi), filter_z_mode);
            else
                ddt(Apar) = filter(ddt(Apar), filter_z_mode);
        } else
            ddt(Ajpar) = filter(ddt(Ajpar), filter_z_mode);

        ddt(U) = filter(ddt(U), filter_z_mode);

        ddt(Ni) = filter(ddt(Ni), filter_z_mode);

        ddt(Ti) = filter(ddt(Ti), filter_z_mode);

        ddt(Te) = filter(ddt(Te), filter_z_mode);

        if (compress0) {
            ddt(Vipar) = filter(ddt(Vipar), filter_z_mode);
        }
    } else if (filter_z_nonlinear) {
        if (!emass) {
            if (evolve_psi)
                ddt(Psi) = filter(ddt(Psi), 0, filter_z_mode);
            else
                ddt(Apar) = filter(ddt(Apar), 0, filter_z_mode);
        } else
            ddt(Ajpar) = filter(ddt(Ajpar), 0, filter_z_mode);

        ddt(U) = filter(ddt(U), 0, filter_z_mode);

        ddt(Ni) = filter(ddt(Ni), 0, filter_z_mode);

        ddt(Ti) = filter(ddt(Ti), 0, filter_z_mode);

        ddt(Te) = filter(ddt(Te), 0, filter_z_mode);

        if (compress0) {
            ddt(Vipar) = filter(ddt(Vipar), 0, filter_z_mode);
        }
    }

#if DEBUG_6F>0
    output.write("I see you 9!\n");//xia
#endif

    if (PF_sink > 0.) {
        if (evolve_psi)
            ddt(Psi) -= PF_sink * sink_PF(F2D_tmp, Psi, PFs_width, PFs_length);
        else
            ddt(Apar) -= PF_sink * sink_PF(F2D_tmp, Apar, PFs_width, PFs_length);
        ddt(Ni) -= PF_sink * sink_PF(N0, Ni, PFs_width, PFs_length);
        ddt(Ti) -= PF_sink * sink_PF(Ti0, Ti, PFs_width, PFs_length);
        ddt(Te) -= PF_sink * sink_PF(Te0, Te, PFs_width, PFs_length);
        ddt(U) -= PF_sink * sink_PF(F2D_tmp, U, PFs_width, PFs_length);
        //	sink_PFtmp = sink_PF(Te0,Te,PFs_width,PFs_length);
    }

    ///////////////////////////////////////////////////////////////////////

    if (low_pass_z > 0) {
        // Low-pass filter, keeping n up to low_pass_z
        if (!emass) {
            if (evolve_psi)
                ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);
            else
                ddt(Apar) = lowPass(ddt(Apar), low_pass_z, zonal_field);
        } else
            ddt(Ajpar) = lowPass(ddt(Ajpar), low_pass_z, zonal_field);

        ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);

        ddt(Ti) = lowPass(ddt(Ti), low_pass_z, zonal_bkgd);
        ddt(Te) = lowPass(ddt(Te), low_pass_z, zonal_bkgd);
        ddt(Ni) = lowPass(ddt(Ni), low_pass_z, zonal_bkgd);

        if (compress0) {
            ddt(Vipar) = lowPass(ddt(Vipar), low_pass_z, zonal_bkgd);
        }

        if (pos_filter) {
            Ti = lowPass_pos2(Ti, Ti);
            Te = lowPass_pos2(Te, Te);
            Ni = lowPass_pos2(Ni, Ni);
        }

        if (pos_filter2) {
            ddt(Ti) = lowPass_pos(ddt(Ti), filter_position_ti);
            ddt(Te) = lowPass_pos(ddt(Te), filter_position_te);
            ddt(Ni) = lowPass_pos(ddt(Ni), filter_position_ni);
        }

        if (pos_filter_zf) {
            ddt(Ti) = sink_zonal_core(ddt(Ti), filter_position_ti);
            ddt(Te) = sink_zonal_core(ddt(Te), filter_position_te);
            ddt(Ni) = sink_zonal_core(ddt(Ni), filter_position_ni);
        }
        if (pos_sink_zf > 0.) {
            ddt(Ti) -= pos_sink_zf * (sink_tanhxl(Ti0, Ti, filter_position_ti * pos_filter_width / int(Grid_NX), filter_position_ti *pos_filter_length / int(Grid_NX))).DC();
            ddt(Te) -= pos_sink_zf * (sink_tanhxl(Te0, Te, filter_position_te * pos_filter_width / int(Grid_NX), filter_position_te *pos_filter_length / int(Grid_NX))).DC();
            ddt(Ni) -= pos_sink_zf * (sink_tanhxl(N0, Ni, filter_position_ni * pos_filter_width / int(Grid_NX), filter_position_ni *pos_filter_length / int(Grid_NX))).DC();

            /* ddt(Ti) -= pos_sink_zf*pos_sink_ti;
         ddt(Te) -= pos_sink_zf*pos_sink_te;
         ddt(Ni) -= pos_sink_zf*pos_sink_ni; */
        }

        if (nonlinear && pos_filter) {
            Pi = lowPass_pos2(Pi, Pi);
            Pe = lowPass_pos2(Pe, Pe);
            P = lowPass_pos2(P, P);
        }
        if (nonlinear && pos_filter2) {
            Pi = lowPass_pos(Pi, position_tmpi);
            Pe = lowPass_pos(Pe, position_tmpe);
            P = lowPass_pos(P, position_tmp);
        }

        if (nonlinear && (pos_filter2 || pos_filter_zf)) {
            Pi = sink_zonal_core(Pi, position_tmpi);
            Pe = sink_zonal_core(Pe, position_tmpe);
            P = sink_zonal_core(P, position_tmp);
        }
    }

    if (damp_width > 0) {
        for (jx = 0; jx < damp_width; jx++) {
            for (jy = 0; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    if (mesh->firstX())
                        ddt(U)[jx][jy][jz] -= U[jx][jy][jz] / damp_t_const;
                    if (mesh->lastX())
                        ddt(U)[mesh->ngx - 1 - jx][jy][jz] -= U[mesh->ngx - 1 - jx][jy][jz] / damp_t_const;
                }
        }
    }

    if (filter_nl > 0) {
        ddt(Ni) = nl_filter(ddt(Ni), filter_nl);
    }
#if DEBUG_6F>0
    output.write("I see you 10!\n");//xia
#endif

    if (fakerun) {
        // output.write("I see you 11!\n");//xia
        ddt(U) = 0.;
        ddt(Ni) = 0.;
        if (evolve_psi)
            ddt(Psi) = 0.;
        else
            ddt(Apar) = 0.;
        ddt(Ti) = 0.;
        ddt(Te) = 0.;
        if (compress0)
            ddt(Vipar) = 0.;
        // output.write("I see you 12!\n");//xia
    }

    first_run = false;

    return 0;
}

/****************BOUNDARY FUNCTIONS*****************************/
// Sheath Boundary Conditions on Phi
// Linearized
void SBC_Dirichlet(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) // let the boundary equall to the value next to the boundary
{
    SBC_yup_eq(var, value, PF_limit, PF_limit_range);
    SBC_ydown_eq(var, -value, PF_limit, PF_limit_range);
}

void SBC_Gradpar(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    SBC_yup_Grad_par(var, value, PF_limit, PF_limit_range);
    SBC_ydown_Grad_par(var, -value, PF_limit, PF_limit_range);
}

// Boundary to specified Field3D object
void SBC_yup_eq(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {

    RangeIterator xrup = mesh->iterateBndryUpperY();

    // for(xrup->first(); !xrup->isDone(); xrup->next())
    for (; !xrup.isDone(); xrup++) {
        xind = xrup.ind;
        indx = mesh->XGLOBAL(xind);
        if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
            for (jy = mesh->yend + 1 - Sheath_width; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = value[xind][jy][jz];
                }
        } else if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
            for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->ngy; jy++)
                for (int jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = 0.;
                }
        } else if (!PF_limit) {
            for (jy = mesh->yend + 1 - Sheath_width; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = value[xind][jy][jz];
                }
        }
    }
    // mesh->communicate(var);
}

void SBC_ydown_eq(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {

    RangeIterator xrdn = mesh->iterateBndryLowerY();

    // for(xrdn->first(); !xrdn->isDone(); xrdn->next())
    for (; !xrdn.isDone(); xrdn++) {
        xind = xrdn.ind;
        indx = mesh->XGLOBAL(xind);
        // output.write("Boundary index: %i   Boundary limit: %e   ix_sep: %e\n", indx, ixsep*PF_limit_range, ixsep);
        if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
            for (jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = value[xind][jy][jz];
                }
        } else if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
            for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
                for (int jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = 0.;
                }
        } else if (!PF_limit) {
            for (jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = value[xind][jy][jz];
                }
        }
    }
    // mesh->communicate(var);
}

// Boundary gradient to specified Field3D object
void SBC_yup_Grad_par(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {

    RangeIterator xrup = mesh->iterateBndryUpperY();

    // for(xrup->first(); !xrup->isDone(); xrup->next())
    for (; !xrup.isDone(); xrup++) {
        xind = xrup.ind;
        indx = mesh->XGLOBAL(xind);
        if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
            for (jy = mesh->yend + 1 - Sheath_width; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = var[xind][jy - 1][jz] + mesh->dy[xind][jy] * sqrt(mesh->g_22[xind][jy]) * value[xind][jy][jz];
                }
        } else if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
            for (jy = mesh->yend + 1 - Sheath_width; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = var[xind][jy - 1][jz];
                }
        } else if (!PF_limit) {
            for (jy = mesh->yend + 1 - Sheath_width; jy < mesh->ngy; jy++)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = var[xind][jy - 1][jz] + mesh->dy[xind][jy] * sqrt(mesh->g_22[xind][jy]) * value[xind][jy][jz];
                }
        }
    }
    // mesh->communicate(var);
}

void SBC_ydown_Grad_par(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {

    RangeIterator xrdn = mesh->iterateBndryLowerY();

    // for(xrdn->first(); !xrdn->isDone(); xrdn->next())
    for (; !xrdn.isDone(); xrdn++) {
        xind = xrdn.ind;
        indx = mesh->XGLOBAL(xind);
        if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
            for (jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = var[xind][jy + 1][jz] - mesh->dy[xind][jy] * sqrt(mesh->g_22[xind][jy]) * value[xind][jy][jz];
                }
        }
        if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
            for (jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = var[xind][jy + 1][jz];
                }
        } else if (!PF_limit) {
            for (jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
                for (jz = 0; jz < mesh->ngz; jz++) {
                    var[xind][jy][jz] = var[xind][jy + 1][jz] - mesh->dy[xind][jy] * sqrt(mesh->g_22[xind][jy]) * value[xind][jy][jz];
                }
        }
    }
    // mesh->communicate(var);
}

const Field3D BS_ft(const int index) {
    Field3D result;
    if (!ft_simple) {
        Field3D result1;
        result.allocate();
        result1.allocate();
        result1 = 0.;

        BoutReal xlam, dxlam;
        dxlam = 1. / max(B0, true) / index;
        xlam = 0.;

        for (int i = 0; i < index; i++) {
            result1 += xlam * dxlam / sqrt(1. - xlam * B0);
            xlam += dxlam;
        }
        result = 1. - 0.75 * B0 * B0 * result1;
    } else {
        result = 1 - (1 - epsilon) * (1 - epsilon) / sqrt(1 - epsilon * epsilon) / (1 + 1.46 * sqrt(epsilon));
    }

    return result;
}

const Field3D F31(const Field3D input) {
    Field3D result;
    result.allocate();

    result = (1 + 1.4 / (Zi + 1.)) * input;
    result -= 1.9 / (Zi + 1.) * input * input;
    result += 0.3 / (Zi + 1.) * input * input * input;
    result += 0.2 / (Zi + 1.) * input * input * input * input;

    return result;
}

const Field3D F32ee(const Field3D input) {
    Field3D result;
    result.allocate();

    result = (0.05 + 0.62 * Zi) / (Zi * (1 + 0.44 * Zi)) * (input - input * input * input * input);
    result += 1. / (1. + 0.22 * Zi) * (input * input - input * input * input * input - 1.2 * (input * input * input - input * input * input * input));
    result += 1.2 / (1. + 0.5 * Zi) * input * input * input * input;

    return result;
}

const Field3D F32ei(const Field3D input) {
    Field3D result;
    result.allocate();

    result = -(0.56 + 1.93 * Zi) / (Zi * (1 + 0.44 * Zi)) * (input - input * input * input * input);
    result += 4.95 / (1. + 2.48 * Zi) * (input * input - input * input * input * input - 0.55 * (input * input * input - input * input * input * input));
    result -= 1.2 / (1. + 0.5 * Zi) * input * input * input * input;

    return result;
}

const Field3D F33(const Field3D input) {
    Field3D result;
    result.allocate();

    result = 1. - (1. + 0.36 / Zeff) * input;
    result += 0.59 / Zeff * input * input;
    result -= 0.23 / Zeff * input * input * input;

    return result;
}

int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
    // output << "Starting precon" << endl;

    ni_tmp = ddt(Ni);
    ti_tmp = ddt(Ti);
    te_tmp = ddt(Te);
    vi_tmp = ddt(Vipar);
    psi_tmp = ddt(Psi);
    u_tmp = ddt(U);
    ni_tmp.applyBoundary("neumann");
    ti_tmp.applyBoundary("neumann");
    te_tmp.applyBoundary("neumann");
    vi_tmp.applyBoundary("neumann");
    psi_tmp.applyBoundary("neumann");
    u_tmp.applyBoundary("neumann");

    // first metrix

    mesh->communicate(ni_tmp, ti_tmp, te_tmp, psi_tmp);
    jpar1 = -B0 * Delp2(psi_tmp);
    mesh->communicate(vi_tmp, u_tmp, jpar1);
    jpar1.applyBoundary();
    if (smooth_j_x) {
        jpar1 = smooth_x(jpar1);
        jpar1 = smooth_x(jpar1);
        jpar1 = smooth_x(jpar1);
        jpar1 = smooth_x(jpar1);
        jpar1 = smooth_x(jpar1);
        jpar1 = smooth_y(jpar1);
        jpar1 = smooth_y(jpar1);
        jpar1 = smooth_y(jpar1);
        jpar1 = smooth_y(jpar1);
        mesh->communicate(jpar1);
        jpar1.applyBoundary();
    }
    // output << "jpar1: " << min(jpar1) << " : " << max(jpar1) << endl;

    p_tmp = N0 * (Tau_ie * ti_tmp + te_tmp) + ni_tmp * (Tau_ie * Ti0 + Te0);
    u_tmp1 = u_tmp + gamma * ((B0 ^ 2) * Grad_par(jpar1 / B0, CELL_CENTRE) + 2.0 * Upara1 * b0xcv * Grad(p_tmp)); // curvature term
    mesh->communicate(u_tmp1);
    u_tmp1.applyBoundary();
    Vipar = vi_tmp;
    mesh->communicate(Vipar);
    Vipar.applyBoundary();
    // output << "u_tmp1: " << min(u_tmp1) << " : " << max(u_tmp1) << endl;

    // second metrix
    if (!diffusion_par) {
        kappa_par_i_lin = 0.;
        kappa_par_e_lin = 0.;
    }

    static InvertPar *invi = 0;
    if (!invi) {
        invi = InvertPar::Create(); // Create parallel solver
        invi->setCoefA(1.0);
    }
    invi->setCoefB(-gamma * kappa_par_i_lin);
    ti_tmp2 = invi->solve(ti_tmp); // Solve Pshur
    mesh->communicate(ti_tmp2);
    ti_tmp2.applyBoundary();
    // output << "ti_tmp2: " << min(ti_tmp2) << " : " << max(ti_tmp2) << endl;
    static InvertPar *inve = 0;
    if (!inve) {
        inve = InvertPar::Create(); // Create parallel solver
        inve->setCoefA(1.0);
    }
    inve->setCoefB(-gamma * kappa_par_e_lin);
    te_tmp2 = inve->solve(te_tmp); // Solve Pshur
    mesh->communicate(te_tmp2);
    te_tmp2.applyBoundary();
    // output << "te_tmp2: " << min(te_tmp2) << " : " << max(te_tmp2) << endl;
    static InvertPar *invu = 0;
    if (!invu) {
        invu = InvertPar::Create(); // Create parallel solver
        Field2D rb_tmp = Grad_par(Rxy * Bpxy);
        mesh->communicate(rb_tmp);
        rb_tmp.applyBoundary("dirichlet");
        invu->setCoefA(1.0 + 2. * gamma * gamma * Grad_par(rb_tmp) / (Rxy * Bpxy * SQ(B0)));
    }
    invu->setCoefB(-SQ(gamma * B0));
    U = invu->solve(u_tmp1); // Solve Pshur
    mesh->communicate(U);
    U.applyBoundary();
    // output << "U: " << min(U) << " : " << max(U) << endl;

    // third metrix
    BoutReal Ntemp = max(N0, true);
    // TODO: diamag in U?
    phi_tmp = invert_laplace(U * B0 / Ntemp, phi_flags, NULL);
    mesh->communicate(phi_tmp);
    phi_tmp.applyBoundary();
    Ni = ni_tmp - gamma * bracket(phi_tmp, N0, bm_exb);
    Ti = ti_tmp2 - gamma * bracket(phi_tmp, Ti0, bm_exb);
    Te = te_tmp2 - gamma * bracket(phi_tmp, Te0, bm_exb);
    Psi = psi_tmp - gamma * Grad_par(phi_tmp) / B0;
    mesh->communicate(Ni, Ti, Te, Psi);
    Ni.applyBoundary();
    Ti.applyBoundary();
    Te.applyBoundary();
    Psi.applyBoundary();
    // output << "Ni: " << min(Ni) << " : " << max(Ni) << endl;
    // output << "Ti: " << min(Ti) << " : " << max(Ti) << endl;
    // output << "Te: " << min(Te) << " : " << max(Te) << endl;
    // output << "Psi: " << min(Psi) << " : " << max(Psi) << endl;

    return 0;
}
