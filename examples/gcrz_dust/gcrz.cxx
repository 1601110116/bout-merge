/*******************************************************************************
*guiding-center module for divertor configuration(with X-points) in RZ coordinate
 *******************************************************************************/
#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <iostream>
#include <fstream>
#include <invert_parderiv.hxx>
#include <sourcex.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>
#include <bout/globalfield.hxx>
//#include "interp2d.cxx"
//#include <bout/globalfield.hxx>
#include "cylinder.cxx"
#include <bout/globalfield.hxx>
#include <string>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_integration.h>

//constant
const BoutReal Mp = 1.6726e-27; // proton mass
const BoutReal E_charge = 1.6022e-19; // elementary charge
const BoutReal PI=3.141592653589793238462643383279502884197;
const BoutReal TWOPI=2.0*PI;
const int neq=6;
const int kr=4,kz=4,kp=4,kpsi=4,ktheta=4;
const BoutReal epsilon_map = 1.0e-6;
const int niterv = 10;
/**** test ****/
//BoutReal gsl_test=gsl_sf_lambert_W0(1.0);
BoutReal gammaa,gammae,kusaia,kusaie,kusai,vd,deltath,deltasec,meanre;
BoutReal md,ii,ie,see,aaa,bbb,ccc,submd,submd1,submd11;
BoutReal const me=9.109382616e-31,mi=1.6726231e-27,kb=1.0,kb1=1.601554e-19,kbther=0.861734e-4,qe=1.6e-19,qi=1.6e-19,e0=8.8542e-12,boltz=1.3806488e-23,daa=1.20173e6,segma=5.6703e-8;
BoutReal const speh =132.02,lath=6.94e16,wf=4.55,row=1.93e4,tmelt=6203;

BoutReal delta,inte1,inte11,inte2,inte22,inte10,inte110,ener,etf,re,jth,je,Vp1,pottl;
//initial parameters
BoutReal const td0=300.0,rd0=1.e-10;
BoutReal rd1=rd0;
/*************/


//parameters
BoutReal Ti,v_th,bmag,rmag_rho,rmag,epsilon_l,epsilon_g,v_para,v_perp,currenttime=0.,vmag,tmag,R_limiter, Z_limiter,t;
BoutReal yperiod;
BoutReal BB,para;
BoutReal Fdrag;
//interpolation schemeGlobalField2D
int interpn;
//Electric field inclued or not
int flagEle,flagphi;
//particle's mass, charge, magnetic moment. unchanged
BoutReal AA,ZZ,mu,energy,y[6],Pzeta;
BoutReal r1,z1,phi1,vpara1,mu1,energy1,Pzeta1;
BoutReal r2,z2,phi2,vpara2,mu2,energy2,Pzeta2;
BoutReal td,rd,time1,time2,te1,ti1,ni1,ne1,te2,ti2,ne2,ni2,vi_para;
//Field2D Ti1,Ni1,Te1,Ne1;
//Field2D Ti0,Ni0,Te0,Ne0;
//Field2D speh;
BoutReal test1,test2,test3,test4,test5,test6,test7,test8,test9;
BoutReal test11,test12,test13,test14,test15,test16,test17,test18,test19;
BoutReal test21,test22,test23,test24,test25,test26,test27,test28,test29;

int region;

//particle's rho_para, position vector. Which will be evolved.
BoutReal rho_para,rr,rz,rp,orbitR,orbitZ,dt,B_init;
int nresult;

//field quantities
int nx,ny,ngx,ngy,ngz,zperiod,archive,pex,pey,nxsub,nysub,xstart,ystart,nxpe,nype;
int jyseps1_1,jyseps2_2,ixseps1,MYPE,tsize;
BoutReal zlength;
Field2D Rxy,Zxy;
BoutReal psia,psib;

Field3D fake;

	mesh1D meshpsi;
	scal1D fpol,pres,workk1,workk2,qpsi;
	meshRZ gmesh,solcoremesh,pf1mesh,pf2mesh,gmesh_fine;
	BoutReal rmin_fine,rmax_fine,zmin_fine,zmax_fine,dr_fine,dz_fine;
	int nr_fine,nz_fine;
        scalRZ psip,bphi,scal1,scal2,scal3,scal4,scal5,scal6,rsolcore,zsolcore,rpf1,zpf1,rpf2,zpf2,phi0rz,TTN1,TTN2,TTN3;
	vecRZ vec_b,vec_gradb,vec_curlb,vec_bxgradb,vec_gradphi,vec_bxgradphi;
	//int npsi_core,ntheta_core,npsi_sol,ntheta_sol,npsi_pf,ntheta_pf;
	int npsi_solcore,ntheta_solcore,npsi_pf1,ntheta_pf1,npsi_pf2,ntheta_pf2;
	int nlcfs,nwall,ndomain,ncore,npf,npf1,npf2,nsolcore;
	int merr_solcore,merr_pf1,merr_pf2;
	int lost_pf1,lost_pf2,lost_sol1,lost_sol2,lost_sol3,lost_sol4,lost_total,tmp_core,tmp_blank,lost_sol5,lost_blank,lost_other;
	int lost_pf1_rec,lost_pf2_rec,lost_sol1_rec,lost_sol2_rec,lost_sol3_rec,lost_sol4_rec,lost_total_rec,tmp_core_rec,tmp_blank_rec;
	int lost_sol5_rec,lost_blank_rec,lost_other_rec;
	BoutReal err_psi;
	
	//BoutReal *gridpsi_core,*gridtheta_core,*gridpsi_sol,*gridtheta_sol,*gridpsi_pf,*gridtheta_pf,dtheta,dpsi;
	//BoutReal *r_lcfs,*z_lcfs,*r_wall,*z_wall,*r_domain,*z_domain,*r_core,*z_core,*r_pf,*z_pf;
	BoutReal *gridpsi_solcore,*gridtheta_solcore,*gridpsi_pf1,*gridtheta_pf1,*gridpsi_pf2,*gridtheta_pf2,dpsi,dtheta;
	BoutReal *r_lcfs,*z_lcfs,*r_core,*z_core;
	BoutReal *r_solcore,*z_solcore,*r_wall,*z_wall,*r_domain,*z_domain,*r_pf1,*z_pf1,*r_pf2,*z_pf2;

int npic;
int indexrun=0;	
particle *impurity;
particle rznode[132][64];
BoutReal **density,**density_rec;
int flag_tb,nstep_tb;
BoutReal D_tb,length_tb,time_tb,dl_tb,dt_tb;


BoutReal interp3d ( BoutReal *xarray, int nx, BoutReal *yarray, int ny, BoutReal *zarray, int nz, Field3D &fxyz, BoutReal x, BoutReal y, BoutReal z);
void gcf(BoutReal t, BoutReal *y, int rkneq, BoutReal *dydt);
void rk4(BoutReal *y, int rkneq, BoutReal t, BoutReal dt, BoutReal *yout);
BoutReal gaussian();
int pnpoly(int nvert, BoutReal *vertx, BoutReal *verty, BoutReal testx, BoutReal testy);
void checkpsi(GlobalField2D &gpsixy);
void domain(GlobalField2D &gRxy, GlobalField2D &gZxy, int nv, BoutReal *vr, BoutReal *vz);

void region_y(BoutReal x, BoutReal &y, int indexregion);
void region_locate(BoutReal x, BoutReal y, int &indexregion);
void map_core(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek);
void map_sol(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek);
void map_pf(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek);
void maprz2pt(particle &impurity);
void randomwalk(BoutReal dr, BoutReal dz, BoutReal &dr_tb, BoutReal &dz_tb);

BoutReal dot (BoutReal v1r, BoutReal v1p, BoutReal v1z, BoutReal v2r,BoutReal v2p,BoutReal v2z) {
	return (v1r*v2r+v1p*v2p+v1z*v2z);
}

BoutReal dot (BoutReal v1r, BoutReal v1p, BoutReal v1z) {
	return (v1r*v1r+v1p*v1p+v1z*v1z);
}

cylinder cross (BoutReal v1r, BoutReal v1p, BoutReal v1z, BoutReal v2r,BoutReal v2p,BoutReal v2z) {
	cylinder res;
	res.r = v1p*v2z-v1z*v2p;
	res.p = v1z*v2r-v1r*v2z;
	res.z = v1r*v2p-v1p*v2r;
	return res;
}

BoutReal *gridr,*gridz,*gridpsi,*gridr_fine,*gridz_fine;
int nxefit,nyefit,maxstep;
GlobalField2D *gRxy, *gZxy, *gpsixy,*gphish0,*gti0,*gni0,*gvi_para0;
BoutReal simagx,sibdry,va,normphi;
Field2D ti0,ni0,vi_para0;
int physics_init ( bool restarting ) {

	Field2D Bxy,Bpxy,Btxy,hthe,I,psixy,phish0;

	// TEST:
	//output.write("gsl: %f\n", gsl_test);
	// Load metrics
	mesh->get(nx, "nx");
	mesh->get(ny, "ny");
	ngx = mesh->ngx;
	ngy = mesh->ngy;
	ngz = mesh->ngz;
	mesh->get(Rxy, "Rxy");
	mesh->get(Zxy, "Zxy");
	mesh->get(Bxy, "Bxy");
	mesh->get(Bpxy, "Bpxy");
	mesh->get(Btxy, "Btxy");
	mesh->get(hthe, "hthe");
	mesh->get(I,    "sinty");
	mesh->get(psixy, "psixy");
	mesh->get(jyseps1_1, "jyseps1_1");
	mesh->get(jyseps2_2, "jyseps2_2");
	mesh->get(ixseps1, "ixseps1");
	//these values are from nc file. Not mapped. Used for obtaining mapping errors. 
        mesh->get(phish0,"phi");
        mesh->get(ni0,"Niexp");
        mesh->get(ti0,"Tiexp");
        mesh->get(vi_para0,"Vi_par");
     
	// Load normalisation values
        mesh->get(bmag, "bmag");
	mesh->get(rmag, "rmag");

	/*************** READ OPTIONS *************************/
	Options *globalOptions = Options::getRoot();
	Options *options = globalOptions->getSection("gc");
	OPTION(options, AA, 1.0);
	OPTION(options, ZZ, -1.0);
       //OPTION(options, td, 1.0);
       //OPTION(options, rd, 1.0);
        OPTION(options, Ti, 100);
	OPTION(options, v_para,0.8);
	OPTION(options, rr, 1.45);//1.487278);//2.24395);
	OPTION(options, rz, 0.0);//-0.8);//-0.0440347);
	OPTION(options, rp, 0.0);
	OPTION(options, dt, 0.01);
	OPTION(options, maxstep, 1000);
	OPTION(options, interpn, 4);
	OPTION(options, flagEle, 0);
	OPTION(options, flagphi, 0);
	OPTION(options, flag_tb, 0);
	OPTION(options, D_tb, 1.0);
	OPTION(options, npic, 1);
	globalOptions->get("zperiod", zperiod, 1);
	globalOptions->get("archive", archive, 20);

//	v_th = sqrt(2.0*E_charge*Ti/Mp/AA);
	//v_th = sqrt(2.0*E_charge*Ti/Mp);
	v_th = 23.0;
	rmag_rho = Mp*v_th/E_charge/bmag;
	epsilon_l=rmag_rho/rmag;
	epsilon_g=Mp*9.8*rmag/(E_charge*Ti);
        v_perp = sqrt(1.0-v_para*v_para);
	//va=1.523144e07;	//va may vary!!!
//	normphi=va*rmag*bmag/Ti;

	//normalise
	Rxy = Rxy / rmag;
	Zxy = Zxy / rmag;
	Bxy = Bxy / bmag;
	mesh->Bxy = Bxy;
	Bpxy = Bpxy / bmag;
	Btxy = Btxy / bmag;
	hthe = hthe / rmag;
	I = I*bmag*rmag*rmag;
	psixy = psixy/bmag/rmag/rmag;
	mesh->dx = mesh->dx /(rmag*rmag*bmag);
	mesh->get(psia,"psi_axis");
	psia = psia /bmag/rmag/rmag;
	mesh->get(psib,"psi_bndry");
	psib = psib /bmag/rmag/rmag;
	//psixy,psia,psib is normalised to (0,1) futher.
	//psixy = (psixy-psia)/(psib-psia);
	
	/**************** CALCULATE METRICS ******************/
	mesh->g11 = (Rxy*Bpxy)^2;
	mesh->g22 = 1.0 / (hthe^2);
	mesh->g33 = (I^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
	mesh->g12 = 0.0;
	mesh->g13 = -I*mesh->g11;
	mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
	mesh->J = hthe / Bpxy;
	mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
	mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
	mesh->g_33 = Rxy*Rxy;
	mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
	mesh->g_13 = I*Rxy*Rxy;
	mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
	mesh->geometry();
        /*end of reading nc file*/

	gRxy = new GlobalField2D(mesh);
	gZxy = new GlobalField2D(mesh);
	gpsixy = new GlobalField2D(mesh);
	gphish0 = new GlobalField2D(mesh);
        gvi_para0 = new GlobalField2D(mesh);
        gti0 = new GlobalField2D(mesh);
        gni0 = new GlobalField2D(mesh);
	density = new BoutReal*[nx];
	density[0] = new BoutReal[nx*ny];
	for (int i=1;i<nx;i++) {
		density[i] = density[i-1] + ny;
	}
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++) {
			density[i][j]=0;
	}
		density_rec = new BoutReal*[nx];
		density_rec[0] = new BoutReal[nx*ny];
		for (int i=1;i<nx;i++) {
			density_rec[i] = density_rec[i-1] + ny;
		}
		for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++) {
				density_rec[i][j]=0;
		}
	gRxy->gather(Rxy);
	gZxy->gather(Zxy);
	gpsixy->gather(psixy);
	gphish0->gather(phish0);
        gvi_para0->gather(vi_para0);
	gti0->gather(ti0);
        gni0->gather(ni0);
        dpsi = (psib-psia)/(nx-1);
	dtheta = TWOPI/ny;

	//3 sub meshes. solcore, pf1, pf2.
	npsi_solcore=nx;
	ntheta_solcore=jyseps2_2 - jyseps1_1;
	gridpsi_solcore = new BoutReal[npsi_solcore];
	gridtheta_solcore = new BoutReal[ntheta_solcore];
	npsi_pf1=nx;
	ntheta_pf1=jyseps1_1+1;
	gridpsi_pf1 = new BoutReal[npsi_pf1];
	gridtheta_pf1 = new BoutReal[ntheta_pf1];
	npsi_pf2=nx;
	ntheta_pf2=ny-jyseps2_2-1;
	gridpsi_pf2 = new BoutReal[npsi_pf1];
	gridtheta_pf2 = new BoutReal[ntheta_pf1];	
	
	MPI_Comm_rank(BoutComm::get(), &MYPE);
	if (MYPE!=0) {
		//Initialise solcore, pf1, pf2 meshes except root core, in order to receive them from root core.
		solcoremesh.init(nx,ntheta_solcore,kpsi,ktheta);
		rsolcore.init(solcoremesh);
		zsolcore.init(solcoremesh);
		pf1mesh.init(nx,ntheta_pf1,kpsi,ktheta);
		rpf1.init(pf1mesh);
		zpf1.init(pf1mesh);
		pf2mesh.init(nx,ntheta_pf2,kpsi,ktheta);
		rpf2.init(pf2mesh);
		zpf2.init(pf2mesh);
	}

	if (MYPE==0) {
		//construct solcore mesh	
		for(int i=0;i<nx;i++) {
			gridpsi_solcore[i] = (*gpsixy)(i,jyseps1_1+1);
		}
		for(int j=0;j<ntheta_solcore;j++) {
			gridtheta_solcore[j] = dtheta*(j+jyseps1_1+1);
		}
		solcoremesh.init(nx,gridpsi_solcore,ntheta_solcore,gridtheta_solcore,kpsi,ktheta);
		rsolcore.init(solcoremesh);
		zsolcore.init(solcoremesh);		
		for(int i=0;i<nx;i++) {
			for(int j=0;j<ntheta_solcore;j++) {
				rsolcore.f[i][j] = (*gRxy)(i,jyseps1_1+j+1);
				zsolcore.f[i][j] = (*gZxy)(i,jyseps1_1+j+1);
			}
		}
		rsolcore.coef(solcoremesh);
		zsolcore.coef(solcoremesh);
		//construct PF1 mesh
		for(int i=0;i<nx;i++) {
			gridpsi_pf1[i] = (*gpsixy)(i,0);
		}
		for(int j=0;j<ntheta_pf1;j++) {
			gridtheta_pf1[j] = dtheta*j;
		}
		pf1mesh.init(nx,gridpsi_pf1,ntheta_pf1,gridtheta_pf1,kpsi,ktheta);
		rpf1.init(pf1mesh);
		zpf1.init(pf1mesh);		
		for(int i=0;i<nx;i++) {
			for(int j=0;j<ntheta_pf1;j++) {
				rpf1.f[i][j] = (*gRxy)(i,j);
				zpf1.f[i][j] = (*gZxy)(i,j);
			}
		}
		rpf1.coef(pf1mesh);
		zpf1.coef(pf1mesh);
		//construct PF2 mesh
		for(int i=0;i<nx;i++) {
			gridpsi_pf2[i] = (*gpsixy)(i,ny-1);
		}
		for(int j=0;j<ntheta_pf2;j++) {
			gridtheta_pf2[j] = dtheta*(jyseps2_2+1+j);
		}
		pf2mesh.init(nx,gridpsi_pf2,ntheta_pf2,gridtheta_pf2,kpsi,ktheta);
		rpf2.init(pf2mesh);
		zpf2.init(pf2mesh);		
		for(int i=0;i<nx;i++) {
			for(int j=0;j<ntheta_pf2;j++) {
				rpf2.f[i][j] = (*gRxy)(i,jyseps2_2+1+j);
				zpf2.f[i][j] = (*gZxy)(i,jyseps2_2+1+j);
			}
		}
		rpf2.coef(pf2mesh);
		zpf2.coef(pf2mesh);
	}

	output.write("debug. befor broadcast 3 meshes.\n");
	//Broadcast 3 sub meshes: solcore, pf1 and pf2
	MPI_Bcast(&solcoremesh.gridr[0],solcoremesh.nr,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&solcoremesh.gridz[0],solcoremesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&solcoremesh.knotr[0],solcoremesh.nr+solcoremesh.kr,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&solcoremesh.knotz[0],solcoremesh.nz+solcoremesh.kz,MPI_DOUBLE,0,BoutComm::get());
	output.write("debug. after broadcast solcore mesh.\n");

	MPI_Bcast(&pf1mesh.gridr[0],pf1mesh.nr,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&pf1mesh.gridz[0],pf1mesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&pf1mesh.knotr[0],pf1mesh.nr+pf1mesh.kr,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&pf1mesh.knotz[0],pf1mesh.nz+pf1mesh.kz,MPI_DOUBLE,0,BoutComm::get());
	output.write("debug. after broadcast pf1 mesh.\n");

	MPI_Bcast(&pf2mesh.gridr[0],pf2mesh.nr,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&pf2mesh.gridz[0],pf2mesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&pf2mesh.knotr[0],pf2mesh.nr+pf2mesh.kr,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&pf2mesh.knotz[0],pf2mesh.nz+pf2mesh.kz,MPI_DOUBLE,0,BoutComm::get());
	output.write("debug. after broadcast pf2 mesh.\n");

	//Broadcast R,Z bscoef in the 3 sub meshes.
	MPI_Bcast(&rsolcore.bscoef[0][0],solcoremesh.nr*solcoremesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&zsolcore.bscoef[0][0],solcoremesh.nr*solcoremesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&rpf1.bscoef[0][0],pf1mesh.nr*pf1mesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&zpf1.bscoef[0][0],pf1mesh.nr*pf1mesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&rpf2.bscoef[0][0],pf2mesh.nr*pf2mesh.nz,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&zpf2.bscoef[0][0],pf2mesh.nr*pf2mesh.nz,MPI_DOUBLE,0,BoutComm::get());
	output.write("debug. after broadcast 3 meshes.\n");

/*	//solcore
	for(int i=0;i<solcoremesh.nr;i++) {
		output.write("%d,%e\n",i,solcoremesh.gridr[i]);
	}
	for(int i=0;i<solcoremesh.nz;i++) {
		output.write("%d,%e\n",i,solcoremesh.gridz[i]);
	}
	for(int i=0;i<solcoremesh.nr+solcoremesh.kr;i++) {
		output.write("%d,%e\n",i,solcoremesh.knotr[i]);
	}
	for(int i=0;i<solcoremesh.nz+solcoremesh.kz;i++) {
		output.write("%d,%e\n",i,solcoremesh.knotz[i]);
	}
	for(int i=0;i<solcoremesh.nr;i++)
		for(int j=0;j<solcoremesh.nz;j++) {
			output.write("%d,%d,%e,%e\n",i,j,rsolcore.bscoef[i][j],zsolcore.bscoef[i][j]);
	}
	//pf1
	for(int i=0;i<pf1mesh.nr;i++) {
		output.write("%d,%e\n",i,pf1mesh.gridr[i]);
	}
	for(int i=0;i<pf1mesh.nz;i++) {
		output.write("%d,%e\n",i,pf1mesh.gridz[i]);
	}
	for(int i=0;i<pf1mesh.nr+pf1mesh.kr;i++) {
		output.write("%d,%e\n",i,pf1mesh.knotr[i]);
	}
	for(int i=0;i<pf1mesh.nz+pf1mesh.kz;i++) {
		output.write("%d,%e\n",i,pf1mesh.knotz[i]);
	}
	for(int i=0;i<pf1mesh.nr;i++)
		for(int j=0;j<pf1mesh.nz;j++) {
			output.write("%d,%d,%e,%e\n",i,j,rpf1.bscoef[i][j],zpf1.bscoef[i][j]);
	}
	//pf2
	for(int i=0;i<pf2mesh.nr;i++) {
		output.write("%d,%e\n",i,pf2mesh.gridr[i]);
	}
	for(int i=0;i<pf2mesh.nz;i++) {
		output.write("%d,%e\n",i,pf2mesh.gridz[i]);
	}
	for(int i=0;i<pf2mesh.nr+pf2mesh.kr;i++) {
		output.write("%d,%e\n",i,pf2mesh.knotr[i]);
	}
	for(int i=0;i<pf2mesh.nz+pf2mesh.kz;i++) {
		output.write("%d,%e\n",i,pf2mesh.knotz[i]);
	}
	for(int i=0;i<pf2mesh.nr;i++)
		for(int j=0;j<pf2mesh.nz;j++) {
			output.write("%d,%d,%e,%e\n",i,j,rpf2.bscoef[i][j],zpf2.bscoef[i][j]);
	}
*/


	//BoutReal tmpb,energy,Pzeta,tmpbt,tmpphi,tmpphi0,tmpphi1;
	ifstream inFile;

	//----------------------------------------------read gfile----------------------------------------------------------
	int pos,pos2;
	BoutReal xdim,zdim,rcentr,rgrid1,zmid,rmagx,zmagx,bcentr,cpasma,xdum;
 	istringstream stream1, stream2;

	string line;

	ifstream gfile ("data/g056129.05550_1");

	if(!gfile) {
		cout << "Error opening the file.\n";
		return 1;
	}

	//The last two elements of fi line are read into nxefit, nyefit
	getline(gfile,line);
	pos = line.find_last_of(" ");
  	stream1.str(line.substr(pos+1));
	stream1 >> nyefit;
	pos2 = line.find_last_of(" ",pos-1);
	stream2.str(line.substr(pos2,pos-pos2));
	stream2 >> nxefit;
	//cout << nxefit << "\t" << nyefit << "\n" ;

	gfile >> xdim >> zdim >> rcentr >> rgrid1 >> zmid;	//line 2
	gfile >> rmagx >> zmagx >> simagx >> sibdry >> bcentr;	//line 3
	gfile >> cpasma >> simagx >> xdum >> rmagx >> xdum;	//line 4
	gfile >> zmagx >> xdum >> sibdry >> xdum >>xdum;	//line 5

	/* The means of each variables
	xdim,zdim	; Size of the domain in meters
	rcentr,bcentr	; Reference vacuum toroidal field (m, T)
	rgrid1		; R of left side of domain
	zmid		; Z at the middle of the domain
	rmagx,zmagx	; Location of magnetic axis
	simagx		; Poloidal flux at the axis (Weber / rad)
	sibdry		; Poloidal flux at plasma boundary (Weber / rad)
	*/

	gridr = new BoutReal [nxefit];
	gridz = new BoutReal [nyefit];
	gridpsi = new BoutReal [nxefit];
	//rmag,bmag is read from nc file. They are keep the same in PIC part.
	//generate the grid points on R and Z direction, and normolised by rmag,bmag which read from nc file.
	for(int i=0;i<nxefit;i++) {gridr[i] = (rgrid1 + xdim*i/(nxefit-1))/rmag;}
	for(int j=0;j<nyefit;j++) {gridz[j] = ((zmid-0.5*zdim)+zdim*j/(nyefit-1))/rmag;}
	//generate the uniform flux grid points
	for(int i=0;i<nxefit;i++) {gridpsi[i] = (simagx+(sibdry-simagx)*i/(nxefit-1))/bmag/rmag/rmag;}
	
        //In order to desribe electric potential accurately, high resolution fine mesh is generated besides gmesh.
	rmin_fine=1.3/rmag;
	rmax_fine=2.3/rmag;
	zmin_fine=-1.2/rmag;
	zmax_fine=0.8/rmag;
	
        nr_fine=601;
	nz_fine=1201;
	dr_fine=(rmax_fine-rmin_fine)/(nr_fine-1);
	dz_fine=(zmax_fine-zmin_fine)/(nz_fine-1);
	gridr_fine = new BoutReal [nr_fine];
	gridz_fine = new BoutReal [nz_fine];
	for(int i=0;i<nr_fine;i++) {gridr_fine[i] = rmin_fine + dr_fine*i;}
	for(int j=0;j<nz_fine;j++) {gridz_fine[j] = zmin_fine + dz_fine*j;}

	//initial the mesh, scal2D and vec2D variables.
	meshpsi.init(nxefit,gridpsi,kpsi);
	fpol.init(meshpsi);
	pres.init(meshpsi);
	workk1.init(meshpsi);
	workk2.init(meshpsi);
	qpsi.init(meshpsi);
	gmesh.init(nxefit,gridr,nyefit,gridz,kr,kz);
	psip.init(gmesh);
	bphi.init(gmesh);
	scal1.init(gmesh);
	scal2.init(gmesh);
	scal3.init(gmesh);
	scal4.init(gmesh);
	vec_b.init(gmesh);
	vec_gradb.init(gmesh);
	vec_curlb.init(gmesh);
	vec_bxgradb.init(gmesh);

        gmesh_fine.init(nr_fine,gridr_fine,nz_fine,gridz_fine,kr,kz);
	scal5.init(gmesh_fine);
	scal6.init(gmesh_fine);
	vec_gradphi.init(gmesh_fine);
	vec_bxgradphi.init(gmesh_fine);


	for (int i=0; i<nxefit; i++) {
		gfile >> fpol.f[i];				//line 6 to 31 Poloidal current function on uniform flux grid
	}
	for (int i=0; i<nxefit; i++) {
		gfile >> pres.f[i];				//line 32 to 57 Plasma pressure in nt/m^2 on uniform flux grid
	}
	for (int i=0; i<nxefit; i++) {
		gfile >> workk1.f[i];				//line 58 to 83
	}
	for (int i=0; i<nxefit; i++) {
		gfile >> workk2.f[i];				//line 84 to 109
	}
	for (int j=0; j<nyefit; j++) {
		for (int i=0; i<nxefit; i++) {
			gfile >> psip.f[i][j];			//line 110 to 3438 Poloidal flux in Weber/rad on grid points
		}
	}
	for (int i=0; i<nxefit; i++) {
		gfile >> qpsi.f[i];				//line 3439 to 3464 q values on uniform flux grid
	}

	gfile >> nlcfs >> nwall;				//number of plasma boundary points and wall boundary points

	if(nlcfs>0) {
		r_lcfs = new BoutReal [nlcfs];
		z_lcfs = new BoutReal [nlcfs];
		for(int i=0;i<nlcfs;i++) {
			gfile >> r_lcfs[i] >> z_lcfs[i];	//line 3466 to 3502 Plasma boundary
		}
	}

	if(nwall>0) {
		r_wall = new BoutReal [nwall];
		z_wall = new BoutReal [nwall];
		for(int i=0;i<nwall;i++) {
			gfile >> r_wall[i] >> z_wall[i];	//line 3503 to 3523 Wall boundary
		}
	}
	gfile.close();

	//--------------------------------end of read gfile--------------------------------------------------------------------
	output.write("debug1.\n");
	phi0rz.init(gmesh_fine);
        TTN1.init(gmesh_fine);
        TTN2.init(gmesh_fine);
        TTN3.init(gmesh_fine);
        
 	ifstream pfile1 ("data/ttn1.dat");

	if(!pfile1) {
		cout << "Error opening ttn1.dat.\n";
		return 1;
	}


	for (int i=0; i<nr_fine; i++) {	
		for (int j=0; j<nz_fine; j++) {
//			pfile >> phi0rz.f[i][j];				//phirz.dat is phirz[1:nz][1:nr], generated by matlab,
			//phi0rz.f[i][j]=phi0rz.f[i][j]*normphi*0.1;
//			phi0rz.f[i][j]=phi0rz.f[i][j]*normphi*4.0;
                       // phi0rz.f[i][j]=0.0;
                        pfile1 >> TTN1.f[i][j];
         //phi0rz.f[i][j]= 0.1*gridr[i]*gridr[i]+0.2*gridz[j];
			//phi0rz.f[i][j]=phi0rz.f[i][j]*0.1;	//decrease by 0.1
		}
	}

        ifstream pfile2 ("data/ttn2.dat");

        if(!pfile2) {
                cout << "Error opening ttn2.dat.\n";
                return 1;
        }


        for (int i=0; i<nr_fine; i++) {
                for (int j=0; j<nz_fine; j++) {
                        pfile2 >> TTN2.f[i][j];
                }
        }

       ifstream pfile3 ("data/ttn3.dat");

        if(!pfile3) {
                cout << "Error opening ttn3.dat.\n";
                return 1;
        }


        for (int i=0; i<nr_fine; i++) {
                for (int j=0; j<nz_fine; j++) {
                        pfile3 >> TTN3.f[i][j];
                  }
        }

      ifstream pfile4 ("data/ttn4.dat");

         if(!pfile4) {
                cout << "Error opening ttn4.dat.\n";
                return 1;
         
         }

        for (int i=0; i<nr_fine; i++) {
                for (int j=0; j<nz_fine; j++) {
                        pfile4 >> phi0rz.f[i][j];
                  }
        }

	output.write("debug2.\n");
	phi0rz.coef(gmesh_fine);
        TTN1.coef(gmesh_fine);
        TTN2.coef(gmesh_fine);
        TTN3.coef(gmesh_fine);

	output.write("debug3.\n");

	if (MYPE==0) {
		BoutReal tmp,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8;
		int i,j;
		ofstream debugfile ("data/errortest.dat");
              // ofstream debugfile ("data/xphi_090106SBC_Li_analytical3_fine.dat");
		//modify gphish0 in core region as genphi.pro doing.
		//for (i=0;i<100;i++) {
		//	for (j=8;j<56;j++) {
		//		(*gphish0)(i,j) = (*gphish0)(100,j);
		//	}
		//}
		for (j=0;j<ny;j++) {
			for (i=0;i<nx;i++) {
				//change original phish0 to analitical one.
				//(*gphish0)(i,j) =0.1*(*gRxy)(i,j)*(*gRxy)(i,j)+0.2*(*gZxy)(i,j);
				//tmp=phi0rz.val((*gRxy)(i,j),(*gZxy)(i,j),0,0,gmesh)/normphi*10.0; //normphi;
				//tmp=phi0rz.val((*gRxy)(i,j),(*gZxy)(i,j),0,0,gmesh)/normphi*0.25;
				tmp = phi0rz.val((*gRxy)(i,j),(*gZxy)(i,j),0,0,gmesh_fine);
                                tmp2 = ((*gphish0)(i,j)-tmp) / (*gphish0)(i,j);
                                debugfile << (*gRxy)(i,j) <<" "<< (*gZxy)(i,j)  <<" "<< (*gphish0)(i,j)  <<" "<< tmp  <<" "<< tmp2  <<"\n";
			}
		}
	}
	output.write("debug4.\n");
	//--------------------------------end of read phirz--------------------------------------------------------------------

	simagx = simagx/rmag/rmag/bmag;
	sibdry = sibdry/rmag/rmag/bmag;
	for(int i=0;i<nxefit;i++) {
		fpol.f[i] = fpol.f[i] / bmag /rmag;
	}
	for (int i=0; i<nxefit; i++) {
		for (int j=0; j<nyefit; j++) {
			psip.f[i][j] = psip.f[i][j] / bmag /rmag /rmag;
		}
	}
	//r_lcfs,z_lcfs: R,Z coordinates of polygon describing separatrix.
	for(int i=0;i<nlcfs;i++) {
		r_lcfs[i] = r_lcfs[i]/rmag;
		z_lcfs[i] = z_lcfs[i]/rmag;
	}
	//r_wall,z_wall: R,Z coordinates of polygon describing wall.
	for(int i=0;i<nwall;i++) {
		r_wall[i] = r_wall[i]/rmag;
		z_wall[i] = z_wall[i]/rmag;
	}
	//computation domain of MHD part is described by ndomain vertex, which's R,Z coodinates are r_domain,z_domain.
	ndomain = ny+ny-jyseps2_2+jyseps1_1;
	r_domain = new BoutReal[ndomain];
	z_domain = new BoutReal[ndomain];

	//normalisation of advancing orbits module
	vmag = sqrt(2.0*E_charge*Ti/Mp);
	rmag_rho = Mp*vmag/E_charge/bmag;
	epsilon_l=rmag_rho/rmag;
	tmag = rmag/vmag;
        epsilon_g=Mp*9.8*rmag/(E_charge*Ti);
	

	R_limiter=R_limiter/rmag;
	Z_limiter=Z_limiter/rmag;

	//flag_tb = 1;
	//length_tb = 1.4e-3;
	//time_tb = 6.0e-7;
	time_tb = dt*tmag;
	length_tb = sqrt(D_tb*time_tb);
	dl_tb = length_tb/rmag;
	dt_tb = time_tb/tmag;
	if(flag_tb ==1) {
		if (dt_tb > 0.1) { 
			output.write("The turbulence frequency is not much higher than particle transit time.\n");
			output.write("time_tb,tmag,dt_tb are: %20.10e%20.10e%20.10e\n",time_tb,tmag,dt_tb);
		}
		else if (dt > dt_tb ) {
			output.write("The time step isn't smaller than turbulence characteristic time.");
			output.write("time step dt,  and dt_tb are(normalized to tmag): %20.10e%20.10e\n",dt,dt_tb);
		}
		nstep_tb = int(dt_tb/dt+0.5);
		dt_tb = float(nstep_tb)*dt;
		time_tb = dt_tb*tmag;
	}
	
	output.write("length_tb(m) and time_tb(sec) are: %20.10e%20.10e\n",length_tb,time_tb);
	output.write("length_tb/rmag_rho and time_tb/tmag/dt are: %20.10e%20.10e\n",length_tb/rmag_rho,dt_tb/dt);
	output.write("%20.10e%5i\n",dt_tb/dt,int(dt_tb/dt+0.5));
	/*----------------------------------------------------------------------------*/

	dump.add(rmag,"rmag",0);
	dump.add(bmag,"bmag",0);

	fpol.coef(meshpsi);
	psip.coef(gmesh);

	BoutReal tmpfpol;
	for(int i=0;i<nxefit;i++) {
		for(int j=0;j<nyefit;j++) {
			vec_b.vr.f[i][j] =  psip.val(gridr[i],gridz[j],0,1,gmesh)/gridr[i];
			vec_b.vz.f[i][j] = -psip.val(gridr[i],gridz[j],1,0,gmesh)/gridr[i];
			//tmpfpol = inrcsval(gridpsi,cscoeffpol,fold(i,j),0)
			tmpfpol = fpol.val(psip.f[i][j],0,meshpsi);
			if (psip.f[i][j] >= gridpsi[nxefit-1] or psip.f[i][j] <= gridpsi[0]) {
				tmpfpol = fpol.f[nxefit-1];
			}
			vec_b.vp.f[i][j] =  tmpfpol / gridr[i];
			scal1.f[i][j] = sqrt(dot(vec_b.vr.f[i][j],vec_b.vp.f[i][j],vec_b.vz.f[i][j]));
		}
	}

	//calculate coef of bphi on (R,Z)
	for(int i=0;i<nxefit;i++) {
		for(int j=0;j<nyefit;j++) {
			bphi.f[i][j] = vec_b.vp.f[i][j];
		}
	}
	bphi.coef(gmesh);
	scal1.coef(gmesh);

	cylinder tmpv;
	for(int i=0;i<nxefit;i++) {
		for(int j=0;j<nyefit;j++) {
			vec_gradb.vr.f[i][j] = scal1.val(gridr[i],gridz[j],1,0,gmesh);
			vec_gradb.vz.f[i][j] = scal1.val(gridr[i],gridz[j],0,1,gmesh);
			vec_gradb.vp.f[i][j] = 0.0;
			//vec_bxgradb[i][j] = cross( vec_b[i][j], vec_gradb[i][j]); !!!!
			tmpv = cross(vec_b.vr.f[i][j],vec_b.vp.f[i][j],vec_b.vz.f[i][j],vec_gradb.vr.f[i][j],vec_gradb.vp.f[i][j],vec_gradb.vz.f[i][j]);
			vec_bxgradb.vr.f[i][j]=tmpv.r;
			vec_bxgradb.vp.f[i][j]=tmpv.p;
			vec_bxgradb.vz.f[i][j]=tmpv.z;

			vec_gradphi.vr.f[i][j] = phi0rz.val(gridr[i],gridz[j],1,0,gmesh);
			vec_gradphi.vz.f[i][j] = phi0rz.val(gridr[i],gridz[j],0,1,gmesh);
			vec_gradphi.vp.f[i][j] = 0.0;
			tmpv = cross(vec_b.vr.f[i][j],vec_b.vp.f[i][j],vec_b.vz.f[i][j],vec_gradphi.vr.f[i][j],vec_gradphi.vp.f[i][j],vec_gradphi.vz.f[i][j]);
			vec_bxgradphi.vr.f[i][j]=tmpv.r;
			vec_bxgradphi.vp.f[i][j]=tmpv.p;
			vec_bxgradphi.vz.f[i][j]=tmpv.z;
		}
	}


	BoutReal tmp1,tmp2;
	for(int i=0;i<nxefit;i++) {
		for(int j=0;j<nyefit;j++) {
			//vec_curlb(i,j)%r = -cdbbsval2d(knotr,knotz,kr,kz,bscoefbphi,gridr(i),gridz(j),0,1)
			vec_curlb.vr.f[i][j] = -bphi.val(gridr[i],gridz[j],0,1,gmesh);
			//tmp1 = cdbbsval2d(knotr,knotz,kr,kz,bscoepsip,gridr(i),gridz(j),2,0)
			tmp1 = psip.val(gridr[i],gridz[j],2,0,gmesh);
			//tmp2 = cdbbsval2d(knotr,knotz,kr,kz,bscoepsip,gridr(i),gridz(j),0,2)
			tmp2 = psip.val (gridr[i],gridz[j],0,2,gmesh);
                        vec_curlb.vp.f[i][j] = (tmp1+tmp2+vec_b.vz.f[i][j])/gridr[i];
			//vec_curlb(i,j)%z = vec_b(i,j)%p/gridr(i) +  cdbbsval2d(knotr,knotz,kr,kz,bscoefbphi,gridr(i),gridz(j),1,0)
			vec_curlb.vz.f[i][j] = (vec_b.vp.f[i][j]) / gridr[i] +  bphi.val(gridr[i],gridz[j],1,0,gmesh);
		}
	}


	for(int i=0;i<nxefit;i++) {
		for(int j=0;j<nyefit;j++) {
			//(vecB .dot. CurlB / B^2)
			scal2.f[i][j] = dot(vec_b.vr.f[i][j],vec_b.vp.f[i][j],vec_b.vz.f[i][j],vec_curlb.vr.f[i][j],vec_curlb.vp.f[i][j],vec_curlb.vz.f[i][j]) / (scal1.f[i][j] * scal1.f[i][j]);
			//(GradB .dot. vecB)
			scal3.f[i][j] = dot(vec_gradb.vr.f[i][j],vec_gradb.vp.f[i][j],vec_gradb.vz.f[i][j],vec_b.vr.f[i][j],vec_b.vp.f[i][j],vec_b.vz.f[i][j]);

			//(GradB .dot. CurlB)
			scal4.f[i][j] = dot (vec_gradb.vr.f[i][j],vec_gradb.vp.f[i][j],vec_gradb.vz.f[i][j],vec_curlb.vr.f[i][j],vec_curlb.vp.f[i][j],vec_curlb.vz.f[i][j]);
			//GradPhi .dot. vecB)
			}
	}

	//scal1.coef(gmesh);
	scal2.coef(gmesh);
	scal3.coef(gmesh);
	scal4.coef(gmesh);

	vec_b.coef(gmesh);		//vec1
	vec_curlb.coef(gmesh);		//vec2
	vec_bxgradb.coef(gmesh);	//vec3
        //generate the variables which envolve phi(fine mesh).
	cylinder tmpvec_b,tmpvec_curlb;
	for(int i=0;i<nr_fine;i++) {
		for(int j=0;j<nz_fine;j++) {
			vec_gradphi.vr.f[i][j] = phi0rz.val(gridr_fine[i],gridz_fine[j],1,0,gmesh_fine);
			vec_gradphi.vz.f[i][j] = phi0rz.val(gridr_fine[i],gridz_fine[j],0,1,gmesh_fine);
			vec_gradphi.vp.f[i][j] = 0.0;
			tmpvec_b = vec_b.val(gridr_fine[i],gridz_fine[j],0,0,gmesh);
			tmpv = cross(tmpvec_b.r,tmpvec_b.p,tmpvec_b.z,vec_gradphi.vr.f[i][j],vec_gradphi.vp.f[i][j],vec_gradphi.vz.f[i][j]);
			vec_bxgradphi.vr.f[i][j]=tmpv.r;
			vec_bxgradphi.vp.f[i][j]=tmpv.p;
			vec_bxgradphi.vz.f[i][j]=tmpv.z;
                        //GradPhi .dot. vecB)
                         scal5.f[i][j] = dot (vec_gradphi.vr.f[i][j],vec_gradphi.vp.f[i][j],vec_gradphi.vz.f[i][j],tmpvec_b.r,tmpvec_b.p,tmpvec_b.z);
                        //(GradPhi .dot. CurlB)
                        tmpvec_curlb = vec_curlb.val(gridr_fine[i],gridz_fine[j],0,0,gmesh);
			scal6.f[i][j] = dot (vec_gradphi.vr.f[i][j],vec_gradphi.vp.f[i][j],vec_gradphi.vz.f[i][j],tmpvec_curlb.r,tmpvec_curlb.p,tmpvec_curlb.z);
                }
        }

	scal5.coef(gmesh_fine);
	scal6.coef(gmesh_fine);
	vec_bxgradphi.coef(gmesh_fine);
	
	scal1.clearf();
	scal2.clearf();
	scal3.clearf();
	scal4.clearf();
	scal5.clearf();
	scal6.clearf();
	vec_b.clearf();
	vec_curlb.clearf();
	vec_bxgradb.clearf();
	vec_bxgradphi.clearf();
	

	//initial the distribution function of impurity.
	int xindex,yindex;
	//rr=(*gRxy)(120,39)*rmag; gRxy only works in root core!!!
	//rz=(*gZxy)(120,39)*rmag; gZxy only works in root core!!!
	output.write("rr,rz is reset to be: %f%f.\n",rr,rz);
	B_init = scal1.val(rr/rmag,rz/rmag,0,0,gmesh);
	srand (time(NULL));
	output.write("debug 0\n");

	//construct 5 ploygons.
	//polygon: whole computation domain.
	ndomain = ny+jyseps1_1+ny-jyseps2_2;
	r_domain= new BoutReal[ndomain];
	z_domain= new BoutReal[ndomain];
	//polygon: inner boundary of core region
	ncore = jyseps2_2-jyseps1_1;
	r_core = new BoutReal[ncore];
	z_core = new BoutReal[ncore];
	//polygon, sol-core region
	nsolcore = 2*(jyseps2_2-jyseps1_1)+6;
	r_solcore = new BoutReal[nsolcore];
	z_solcore = new BoutReal[nsolcore];
	//polygon, pf1 region
	npf1=2*(jyseps1_1+1)+3;
	r_pf1 = new BoutReal[npf1];
	z_pf1 = new BoutReal[npf1];
	//polygon, pf2 region
	npf2=2*(ny-jyseps2_2-1)+3;
	r_pf2 = new BoutReal[npf2];
	z_pf2 = new BoutReal[npf2];

	if (MYPE==0) {
		//polygon: whole computation domain in MHD part.
		for(int i=0;i<ny;i++) {
			r_domain[i] = (*gRxy)(nx-1,i);
			z_domain[i] = (*gZxy)(nx-1,i);
		}
		for(int i=0;i<ny-jyseps2_2-1;i++) {
			r_domain[ny+i] = (*gRxy)(0,ny-1-i);
			z_domain[ny+i] = (*gZxy)(0,ny-1-i);
		}
		for(int i=0;i<jyseps1_1+1;i++) {
			r_domain[ny+ny-jyseps2_2-1+i] = (*gRxy)(0,jyseps1_1-i);
			z_domain[ny+ny-jyseps2_2-1+i] = (*gZxy)(0,jyseps1_1-i);
		}
		//polygon: inner boundary of core region
		for(int i=0;i<ncore;i++) {
			r_core[i] = (*gRxy)(0,jyseps1_1+i+1);
			z_core[i] = (*gZxy)(0,jyseps1_1+i+1);
		}
		//polygon, sol-core region
		for(int i=0;i<jyseps2_2-jyseps1_1;i++) {
			r_solcore[i] = (*gRxy)(nx-1,jyseps1_1+i+1);
			z_solcore[i] = (*gZxy)(nx-1,jyseps1_1+i+1);
		}
		r_solcore[jyseps2_2-jyseps1_1] = (*gRxy)( (ixseps1+nx-1)/2, jyseps2_2);
		r_solcore[jyseps2_2-jyseps1_1+1] = (*gRxy)( ixseps1, jyseps2_2);
		r_solcore[jyseps2_2-jyseps1_1+2] = (*gRxy)( ixseps1/2, jyseps2_2);
		z_solcore[jyseps2_2-jyseps1_1] = (*gZxy)( (ixseps1+nx-1)/2, jyseps2_2);
		z_solcore[jyseps2_2-jyseps1_1+1] = (*gZxy)( ixseps1, jyseps2_2);
		z_solcore[jyseps2_2-jyseps1_1+2] = (*gZxy)( ixseps1/2, jyseps2_2);
		for(int i=0;i<jyseps2_2-jyseps1_1;i++) {
			r_solcore[jyseps2_2-jyseps1_1+i+3] = (*gRxy)(0,jyseps2_2-i);
			z_solcore[jyseps2_2-jyseps1_1+i+3] = (*gZxy)(0,jyseps2_2-i);
		}
		r_solcore[2*(jyseps2_2-jyseps1_1)+3] = (*gRxy)( ixseps1/2, jyseps1_1+1);
		r_solcore[2*(jyseps2_2-jyseps1_1)+4] = (*gRxy)( ixseps1,jyseps1_1+1);
		r_solcore[2*(jyseps2_2-jyseps1_1)+5] = (*gRxy)( (ixseps1+nx-1)/2, jyseps1_1+1);
		z_solcore[2*(jyseps2_2-jyseps1_1)+3] = (*gZxy)( ixseps1/2, jyseps1_1+1);
		z_solcore[2*(jyseps2_2-jyseps1_1)+4] = (*gZxy)( ixseps1,jyseps1_1+1);
		z_solcore[2*(jyseps2_2-jyseps1_1)+5] = (*gZxy)( (ixseps1+nx-1)/2, jyseps1_1+1);
		//polygon, pf1 region
		for(int i=0;i<jyseps1_1+1;i++) {
			r_pf1[i] = (*gRxy)(nx-1,i);
			z_pf1[i] = (*gZxy)(nx-1,i);
		}
		r_pf1[jyseps1_1+1] = (*gRxy)( (ixseps1+nx-1)/2, jyseps1_1);
		r_pf1[jyseps1_1+2] = (*gRxy)( ixseps1, jyseps1_1);
		r_pf1[jyseps1_1+3] = (*gRxy)( ixseps1/2, jyseps1_1);
		z_pf1[jyseps1_1+1] = (*gZxy)( (ixseps1+nx-1)/2, jyseps1_1);
		z_pf1[jyseps1_1+2] = (*gZxy)( ixseps1, jyseps1_1);
		z_pf1[jyseps1_1+3] = (*gZxy)( ixseps1/2, jyseps1_1);
		for(int i=0;i<jyseps1_1+1;i++) {
			r_pf1[jyseps1_1+4+i] = (*gRxy)(0,jyseps1_1-i);
			z_pf1[jyseps1_1+4+i] = (*gZxy)(0,jyseps1_1-i);
		}
		//polygon, pf2 region
		for(int i=0;i<ny-jyseps2_2-1;i++) {
			r_pf2[i] = (*gRxy)(nx-1,i+jyseps2_2+1);
			z_pf2[i] = (*gZxy)(nx-1,i+jyseps2_2+1);
		}
		for(int i=0;i<ny-jyseps2_2-1;i++) {
			r_pf2[ny-jyseps2_2-1+i] = (*gRxy)(0,ny-1-i);
			z_pf2[ny-jyseps2_2-1+i] = (*gZxy)(0,ny-1-i);
		}
		r_pf2[2*(ny-jyseps2_2-1)]   = (*gRxy)(ixseps1/2,jyseps2_2+1);
		r_pf2[2*(ny-jyseps2_2-1)+1] = (*gRxy)(ixseps1,jyseps2_2+1);
		r_pf2[2*(ny-jyseps2_2-1)+2] = (*gRxy)((ixseps1+nx-1)/2,jyseps2_2+1);
		z_pf2[2*(ny-jyseps2_2-1)]   = (*gZxy)(ixseps1/2,jyseps2_2+1);
		z_pf2[2*(ny-jyseps2_2-1)+1] = (*gZxy)(ixseps1,jyseps2_2+1);
		z_pf2[2*(ny-jyseps2_2-1)+2] = (*gZxy)((ixseps1+nx-1)/2,jyseps2_2+1);
	}

	output.write("debug. befor broadcast polygons.\n");
	//Broadcast polygons.
	MPI_Bcast(&r_domain[0],ndomain,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&z_domain[0],ndomain,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&r_core[0],ncore,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&z_core[0],ncore,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&r_solcore[0],nsolcore,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&z_solcore[0],nsolcore,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&r_pf1[0],npf1,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&z_pf1[0],npf1,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&r_pf2[0],npf2,MPI_DOUBLE,0,BoutComm::get());
	MPI_Bcast(&z_pf2[0],npf2,MPI_DOUBLE,0,BoutComm::get());
	output.write("debug. after broadcast polygons.\n");

/*	output.write("debug. domain:i,r[i],z[i]\n");
	for(int i=0;i<ndomain;i++) {
		output.write("%d,%e,%e\n",i,r_domain[i],z_domain[i]);
	}
	output.write("debug. core:i,r[i],z[i]\n");
	for(int i=0;i<ncore;i++) {
		output.write("%d,%e,%e\n",i,r_core[i],z_core[i]);
	}
	output.write("debug. solcore:i,r[i],z[i]\n");
	for(int i=0;i<nsolcore;i++) {
		output.write("%d,%e,%e\n",i,r_solcore[i],z_solcore[i]);
	}
	output.write("debug. pf1:i,r[i],z[i]\n");
	for(int i=0;i<npf1;i++) {
		output.write("%d,%e,%e\n",i,r_pf1[i],z_pf1[i]);
	}
	output.write("debug. pf2:i,r[i],z[i]\n");
	for(int i=0;i<npf2;i++) {
		output.write("%d,%e,%e\n",i,r_pf2[i],z_pf2[i]);
	}
*/

/*	particle testmap;
	testmap.y[0]=  rr/rmag;	//7.7384763490e-01;
	testmap.y[1]=  rz/rmag;	//2.5233608594e-01;
	testmap.psi = psip.val(testmap.y[0],testmap.y[1],0,0,gmesh);
	testmap.theta=PI;
	testmap.region=0;
	maprz2pt(testmap);
	output.write("%e,%e,%d,%e,%e\n", rr/rmag,rz/rmag,testmap.region,testmap.psi,testmap.theta);

	output.write("launch point R,Z and its psi,theta:\n");
	output.write("%20.10e%20.10e%20.10e%20.10e\n", rr,rz,(testmap.psi-psia)/(psib-psia),testmap.theta);
	output.write("map from RZ nodes to (psi,theta) mesh.\n");
	output.write("    i    j         R         Z       psi     theta idex:\n");
	output.write("debug 0.1 \n");
	for (int i=0;i<nxefit;i++) {
		for(int j=0;j<nyefit;j++) {
			rznode[i][j].y[0] = gridr[i];
			rznode[i][j].y[1] = gridz[j];
			rznode[i][j].psi = psip.val(rznode[i][j].y[0],rznode[i][j].y[1],0,0,gmesh);
			maprz2pt(rznode[i][j]);
			}
	}	
	output.write("debug after rznode. \n");
*/

	merr_solcore=0;
	err_psi = 0.0;
	merr_pf1=0;
	merr_pf2=0;
	tmp_core=0;
	tmp_blank=0;
	lost_total=0;

lost_pf1=0;
lost_pf2=0;
lost_sol1=0;
lost_sol2=0;
lost_sol3=0;
lost_sol4=0;
lost_sol5=0;
lost_blank=0;
lost_other=0;

	BoutReal tmpvth;
	impurity = new particle[npic];
	particle tmpparticle;
	char locname[50];
	sprintf(locname,"data/impurity%d.dat",MYPE);
	if(restarting){
		ifstream fin;
		output.write("start to open data/impurity.dat. \n");
		fin.open(locname,ios_base::in | ios_base::binary);
		if (!fin.is_open()) {
			output.write("fail to open data/impurity.dat. \n");
		}
		output.write("success to open data/impurity.dat. \n");
		for(int i=0;i<npic;i++) {
			fin.read((char*)&impurity[i],sizeof tmpparticle);
		}
		fin.close();
	}
	else {
		for(int i=0;i<npic;i++) {
			//Zxy(*,39) is near middle plane of low field side in this file: g38300_132X64_0.9-1.02.nc.
			impurity[i].y[0] = rr/rmag;
			impurity[i].y[1] = rz/rmag;
			impurity[i].y[2] = rp;
		       //BB = impurity[i].y[2];
                        tmpvth = gaussian();
		       //v_para = sin(PI*(rand()/(BoutReal)(RAND_MAX)-0.5));
			tmpvth = 1.0;
		        // v_para=0.01*(i-99);
			//tmpvth =1.0;// .4;
			//v_para = 0.8;//0.3;
			v_perp = sqrt(1.0-v_para*v_para);
			v_para = v_para*abs(tmpvth);
			v_perp = v_perp*abs(tmpvth);
			impurity[i].y[3] = AA*v_para/ZZ/B_init;
		        impurity[i].y[4] = td0;//300.0/11605.0;//initial dust temperature 300k. Normalised to eV
                       //BB = impurity[i].y[4];
                        impurity[i].y[5] = rd0;
                        impurity[i].mu = AA*v_perp*v_perp/B_init;
			impurity[i].psi = psip.val(impurity[i].y[0],impurity[i].y[1],0,0,gmesh);
			impurity[i].theta = PI;
			maprz2pt(impurity[i]);
		}
	}

	output.write("debug index.\n");
	xstart=mesh->xstart;
	ystart=mesh->ystart;
	nxpe=mesh->getNXPE();
	nype=mesh->getNYPE();
	output.write("%d,%d\n",mesh->xstart,mesh->xend);
	output.write("%d,%d\n",mesh->ystart,mesh->yend);
	output.write("%d,%d\n",mesh->getNXPE(),mesh->getNYPE());
	output.write("%d,%d,%d\n",ngx,ngy,ngz);
	output.write("end debug index.\n");

/*	
	locdens=MYPE+0.1;
	for(int i=0;i<ngx;i++)
		for(int j=0;j<ngy;j++) {
			output.write("%d,%d,%f\n",i,j,locdens[i][j]);
	}

	MPI_Reduce(&density[0][0],&density[0][0],nx*ny,MPI_DOUBLE,MPI_SUM,0,BoutComm::get());
	MPI_Bcast(&density[0][0],nx*ny,MPI_DOUBLE,0,BoutComm::get());
	// Get the X and Y indices
	pex = MYPE % nxpe;
	pey = MYPE / nxpe;
	nxsub=nx/nxpe;
	nysub=ny/nype;
	//scatter date from density[nx][ny] to Field2D locdens;
	output.write("%d,%d,%d,%d\n",nx,ny,ngx,ngy);
	output.write("%d,%d\n",xstart,ystart);
	output.write("%d,%d,%d,%d\n",pex,pey,nxsub,nysub);
	for(int i=0;i<nxsub;i++)
		for(int j=0;j<nysub;j++) {
			//output.write("%d,%d,%d,%d\n",i+xstart,j+ystart,pex*nxsub+i,pey*nysub+j);
			//seems no guard cell in x direction? ngx=nx and ngy=ny+2*2
			//locdens[i+xstart][j+ystart]=density[pex*nxsub+i][pey*nysub+j];
	}
*/
	//density is gathered to root core.

	MPI_Reduce(&density[0][0],&density_rec[0][0],nx*ny,MPI_DOUBLE,MPI_SUM,0,BoutComm::get());
	if (MYPE==0) {
		ofstream fout("data/density.dat",ios_base::out | ios_base::binary);
		for(int j=0;j<ny;j++)
			for(int i=0;i<nx;i++) {
				fout.write((char*)&density_rec[i][j], sizeof (double));
		}
		fout.close();

		ofstream fout2("data/samporbit.dat",ios_base::out | ios_base::binary);
		for(int i=0;i<npic;i++) {
			fout2.write((char*)&impurity[i].y[0], sizeof (double));
			fout2.write((char*)&impurity[i].y[1], sizeof (double));
		}
		fout2.close();
	}

	output.write("r,z is: %10.5f%10.5f",rr/rmag,rz/rmag);
	output.write("r,z is: %10.5f%10.5f",impurity[npic-1].y[0],impurity[npic-1].y[1]);

	y[0] = impurity[npic-1].y[0];
	y[1] = impurity[npic-1].y[1];
        time1 = 0;
 
        //AA=1.0;
        //ZZ=1.0;
        //dump.add(locdens,"locdens",1);
	SOLVE_FOR(fake);
	dump.add(merr_solcore,"merr_solcore",1);	
	dump.add(err_psi,"err_psi",1);
	dump.add(merr_pf1,"merr_pf1",1);
	dump.add(merr_pf2,"merr_pf2",1);

	dump.add(tmp_core_rec,"tmp_core",1);
	dump.add(tmp_blank_rec,"tmp_blank",1);
	dump.add(lost_total_rec,"lost_total",1);
	dump.add(lost_pf1_rec,"lost_pf1",1);
	dump.add(lost_pf2_rec,"lost_pf2",1);
	dump.add(lost_sol1_rec,"lost_sol1",1);
	dump.add(lost_sol2_rec,"lost_sol2",1);
	dump.add(lost_sol3_rec,"lost_sol3",1);
	dump.add(lost_sol4_rec,"lost_sol4",1);
	dump.add(lost_sol5_rec,"lost_sol5",1);
	dump.add(lost_blank_rec,"lost_blank",1);
	dump.add(lost_other_rec,"lost_other",1);

        //dump.add(v_para,"para",1);
        dump.add(y[0],"r",1);
        dump.add(y[1],"z",1);
        dump.add(y[2],"phi",1);
        dump.add(y[3],"vpara",1);
        dump.add(y[4],"td",1);
        dump.add(y[5],"rd",1);
        dump.add(rd1,"rd1",1);
        dump.add(time1,"time1",1);
        dump.add(time2,"time2",1);
        dump.add(AA,"aa",1);
        dump.add(ZZ,"zz",1);
        dump.add(energy,"energy",1);
        dump.add(Pzeta,"pzeta",1);
        dump.add(region,"region",1);

        dump.add(r1,"r1",1);
        dump.add(z1,"z1",1);
        dump.add(phi1,"phi1",1);
        dump.add(vpara1,"vpara1",1);
        dump.add(mu1,"mu1",1);
        dump.add(energy1,"energy1",1);
        dump.add(Pzeta1,"pzeta1",1);
        dump.add(r2,"r2",1);
        dump.add(z2,"z2",1);
        dump.add(phi2,"phi2",1);
        dump.add(vpara2,"vpara2",1);
        dump.add(mu2,"mu2",1);
        dump.add(energy2,"energy2",1);
        dump.add(Pzeta2,"pzeta2",1);

        dump.add(ti1,"ti1",1); 
        dump.add(te1,"te1",1);
        dump.add(ni1,"ni1",1);
        dump.add(ne1,"ne1",1);
        
        dump.add(Vp1,"vp",1);
        dump.add(pottl,"pottl",1);
        
        dump.add(test1,"test1",1);
        dump.add(test11,"test11",1);
        dump.add(test21,"test21",1);
        dump.add(test2,"test2",1);
        dump.add(test12,"test12",1);
        dump.add(test22,"test22",1);
        dump.add(test3,"test3",1);
        dump.add(test13,"test13",1);
        dump.add(test23,"test23",1);
        dump.add(test4,"test4",1);
        dump.add(test14,"test14",1);
        dump.add(test24,"test24",1);
        dump.add(test5,"test5",1);
        dump.add(test15,"test15",1);
        dump.add(test25,"test25",1);
        dump.add(test6,"test6",1);
        dump.add(test16,"test16",1);
        dump.add(test26,"test26",1);
        dump.add(test7,"test7",1);
        dump.add(test17,"test17",1);
        dump.add(test27,"test27",1);
        dump.add(test8,"test8",1);
        dump.add(test9,"test9",1);
        dump.add(test18,"test18",1);
        dump.add(test19,"test19",1);
        dump.add(test28,"test28",1);
        dump.add(test29,"test29",1);
         
        dump.add(ti0,"ti0",1);
        dump.add(ni0,"ni0",1);
        dump.add(kusai,"kusai",1);
        dump.add(kusaia,"kusaia",1);
        dump.add(kusaie,"kusaie",1);
        dump.add(gammaa,"gammaa",1);
        dump.add(gammae,"gammae",1);
        dump.add(vd,"vd",1);
        dump.add(deltath,"deltath",1);
        dump.add(deltasec,"deltasec",1);
        dump.add(meanre,"meanre",1);




return 0;

}

int physics_run(BoutReal t) {
	int tmpindex,xindex,yindex;
	BoutReal x=0.0,tmpb,tmpbt,tmppsi,tmpfpol,yout[6],BB;
	BoutReal dr,dz,dr_tb,dz_tb,tmpb0,tmpb1,tmpfac,tmpv,tmpphi;
	particle tmpparticle;

	MPI_Comm_rank(BoutComm::get(), &MYPE);

	char locname[50];
	sprintf(locname,"data/impurity%d.dat",MYPE);
	if( t>0.0 and fmod(t,1.0) ==0.0) {
	        time1=time1+dt;
                time2=time1*tmag;
                 //if(time2>=tpt){
                 //output.write("Run Finished!") 
                 //exit(EXIT_SUCCESS);}

                     for(int i=0;i<nx;i++)
			for(int j=0;j<ny;j++) {
				density[i][j]=0;
		}
		tmp_blank=0;
		tmp_core=0;
		for(int i=0;i<npic;i++) {
			mu = impurity[i].mu;
			rk4(impurity[i].y, neq, t, dt, yout);
		      
                	if(flag_tb == 1 and fmod(t,BoutReal(nstep_tb))==0.0){
				if(i==npic-1) {
					output.write("The time is:%20.10e\n",t);
				}
		      
                 	dr = yout[0]-impurity[i].y[0];
			dz = yout[1]-impurity[i].y[1];
			randomwalk(dr,dz,dr_tb,dz_tb);
			y[0] = yout[0] + dr_tb;
			y[1] = yout[1] + dz_tb;

			tmpb0 = scal1.val(yout[0],yout[1],0,0,gmesh);
			tmpb1 = scal1.val(y[0],y[1],0,0,gmesh);
			tmpfac = tmpb0/tmpb1;
			impurity[i].mu = impurity[i].mu*tmpfac;
			mu = impurity[i].mu;
			yout[3] = yout[3]*tmpfac;
			yout[0]=y[0];
			yout[1]=y[1];
		       // yout[4]=y[4];
                       //yout[5]=y[5];
                                               
 }
		//y[4]=impurity[i].y[4];
                        BB = impurity[i].y[4];
			impurity[i].y[0] = yout[0];
			impurity[i].y[1] = yout[1];
			impurity[i].y[2] = yout[2];
		        impurity[i].y[3] = yout[3];
		        impurity[i].y[4] = yout[4];
                        impurity[i].y[5] = yout[5];
                       
                         
                        y[0]= yout[0];
			y[1]= yout[1];
		        y[2]= yout[2];
			y[3]= yout[3];
                        y[4]= yout[4];
                        y[5]= yout[5];
                         
                         // BB = impurity[i].y[2];
			tmpb = scal1.val(y[0],y[1],0,0,gmesh);
			tmpbt = bphi.val(y[0],y[1],0,0,gmesh);
			tmppsi =  psip.val(y[0],y[1],0,0,gmesh);
			tmpphi = phi0rz.val(y[0],y[1],0,0,gmesh_fine);
			tmpphi=0.0;
			energy = mu*tmpb + (ZZ*y[3]*tmpb)*(ZZ*y[3]*tmpb)/AA+ZZ*tmpphi;

			tmpfpol = fpol.val(tmppsi, 0, meshpsi);
			if (tmppsi >= gridpsi[nxefit-1] or tmppsi <= gridpsi[0]) {
				tmpfpol = fpol.f[nxefit-1];
			}
			Pzeta = tmppsi-epsilon_l*y[3]*tmpfpol;

			if(i==0) {
				r1 = y[0];
				z1 = y[1];
				phi1 = y[2];
				vpara1 = y[3];
				mu1 = mu;
				energy1 = energy;
				Pzeta1 = Pzeta;
			}

			if(i==158) {
				r2 = y[0];
				z2 = y[1];
				phi2 = y[2];
				vpara2 = y[3];
				mu2 = mu;
				energy2 = energy;
				Pzeta2 = Pzeta;
			}

			maprz2pt(impurity[i]);
		}
                
                
		MPI_Reduce(&density[0][0],&density_rec[0][0],nx*ny,MPI_DOUBLE,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_pf1,&lost_pf1_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_pf2,&lost_pf2_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_sol1,&lost_sol1_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_sol2,&lost_sol2_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_sol3,&lost_sol3_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_sol4,&lost_sol4_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_sol5,&lost_sol5_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_other,&lost_other_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_blank,&lost_blank_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&lost_total,&lost_total_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&tmp_blank,&tmp_blank_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		MPI_Reduce(&tmp_core,&tmp_core_rec,1,MPI_INT,MPI_SUM,0,BoutComm::get());
		if (MYPE==0) {
			ofstream fout("data/density.dat",ios_base::out | ios_base::app | ios_base::binary);
			for(int j=0;j<ny;j++)
				for(int i=0;i<nx;i++) {
						fout.write((char*)&density_rec[i][j], sizeof (double));
			}
			fout.close();

			ofstream fout2("data/samporbit.dat",ios_base::out | ios_base::app | ios_base::binary);
			for(int i=0;i<npic;i++) {
				fout2.write((char*)&impurity[i].y[0], sizeof (double));
				fout2.write((char*)&impurity[i].y[1], sizeof (double));
			}
			fout2.close();
		}
	}

	y[0] = yout[0];
        y[1] = yout[1];
       //y[4]=yout[4];
       // y[5]=yout[5];
	/*if( t>0.0 and fmod(t,1.0) ==0.0) {
	ndes = 0.0;
	gndes->gather(ndes);

	nerr_solcore=0;
	err_psi = 0.0;
	nerr_pf1=0;
	nerr_pf2=0;
	nerr_core=0;
	nerr_blank=0;
	MPI_Comm_rank(BoutComm::get(), &MYPE);
	if (MYPE==0) {
		for(int i=0;i<npic;i++) {
			if(impurity[i].region == 9) continue;
			mu = impurity[i].mu;
			rk4(impurity[i].y, neq, t, dt, yout);
			if(flag_tb == 1 and fmod(t,BoutReal(nstep_tb))==0.0){
				if(i==npic-1) {
					output.write("The time is:%20.10e\n",t);
				}
			dr = yout[0]-impurity[i].y[0];
			dz = yout[1]-impurity[i].y[1];
			randomwalk(dr,dz,dr_tb,dz_tb);
			y[0] = yout[0] + dr_tb;
			y[1] = yout[1] + dz_tb;

			tmpb0 = scal1.val(yout[0],yout[1],0,0,gmesh);
			tmpb1 = scal1.val(y[0],y[1],0,0,gmesh);
		tmpfac = tmpb0/tmpb1;
			impurity[i].mu = impurity[i].mu*tmpfac;
			mu = impurity[i].mu;
			yout[3] = yout[3]*tmpfac;
			yout[0]=y[0];
			yout[1]=y[1];
			}

			y[0] = yout[0];
			y[1] = yout[1];
			y[2] = yout[2];
			y[3] = yout[3];
			impurity[i].y[0] = yout[0];
			impurity[i].y[1] = yout[1];
			impurity[i].y[2] = yout[2];
			impurity[i].y[3] = yout[3];

			tmpb = scal1.val(y[0],y[1],0,0,gmesh);
			tmpbt = bphi.val(y[0],y[1],0,0,gmesh);
			tmppsi =  psip.val(y[0],y[1],0,0,gmesh);
			energy = mu*tmpb + (ZZ*y[3]*tmpb)*(ZZ*y[3]*tmpb)/AA;
			tmpfpol = fpol.val(tmppsi, 0, meshpsi);
			if (tmppsi >= gridpsi[nxefit-1] or tmppsi <= gridpsi[0]) {
				tmpfpol = fpol.f[nxefit-1];
			}
			Pzeta = tmppsi-epsilon_l*y[3]*tmpfpol;
			maprz2pt(impurity[i]);

			region = impurity[i].region;	

		}
	}
	ndes = gndes->scatter();
	indexrun++;
	}*/

	if( t>0.0 and fmod(t,archive) ==0.0) {
	      
                ofstream fout(locname,ios_base::out | ios_base::binary);
		for(int i=0;i<npic;i++) {
			fout.write((char*)&impurity[i], sizeof tmpparticle);
		}
		fout.close();
	}

	ddt(fake) = 0.;
	return 0;
}

void gcf(BoutReal x, BoutReal *y, int neq, BoutReal *dy) {
	BoutReal D,tmpscal1,tmpscal2,tmpscal3,tmpscal4,tmpscal5,tmpscal6;
       	cylinder tmpvec1,tmpvec2,tmpvec3,tmpvec4,tmpvec5;
       
        tmpscal1 = scal1.val(y[0],y[1],0,0,gmesh);      //b
	tmpscal2 = scal2.val(y[0],y[1],0,0,gmesh);      //b_dot_curlb
	tmpscal3 = scal3.val(y[0],y[1],0,0,gmesh);      //gradb_dot_b
	tmpscal4 = scal4.val(y[0],y[1],0,0,gmesh);      //gradb_dot_curlb
	tmpscal5 = scal5.val(y[0],y[1],0,0,gmesh_fine);	//b_dot_gradphi
	tmpscal6 = scal6.val(y[0],y[1],0,0,gmesh_fine);	//curlb_dot_gradphi

	tmpscal5 = 0.0;
	tmpscal6 = 0.0;
        
	tmpvec1 = vec_b.val(y[0],y[1],0,0,gmesh);         //b
	tmpvec2 = vec_curlb.val(y[0],y[1],0,0,gmesh);     //curl_b
	tmpvec3 = vec_bxgradb.val(y[0],y[1],0,0,gmesh);   //b_cross_gradb
	tmpvec4 = vec_bxgradphi.val(y[0],y[1],0,0,gmesh_fine); //b_cross_gradphi

	tmpvec4.r =0.0;
	tmpvec4.z =0.0;
	tmpvec4.p =0.0;

	D = 1.0 + epsilon_l*y[3]*tmpscal2;
	dy[0] = ZZ*y[3]/AA/D*(tmpvec1.r + epsilon_l*y[3]*tmpvec2.r) + epsilon_l*0.5/D/tmpscal1/tmpscal1*((2.0*ZZ*y[3]*y[3]/AA*tmpscal1+mu/ZZ)*tmpvec3.r+tmpvec4.r);
	dy[1] = ZZ*y[3]/AA/D*(tmpvec1.z + epsilon_l*y[3]*tmpvec2.z) + epsilon_l*0.5/D/tmpscal1/tmpscal1*((2.0*ZZ*y[3]*y[3]/AA*tmpscal1+mu/ZZ)*tmpvec3.z+tmpvec4.z);
        dy[2] = ZZ*y[3]/AA/D*(tmpvec1.p + epsilon_l*y[3]*tmpvec2.p) + epsilon_l*0.5/D/tmpscal1/tmpscal1*((2.0*ZZ*y[3]*y[3]/AA*tmpscal1+mu/ZZ)*tmpvec3.p+tmpvec4.p+epsilon_g*AA*tmpvec1.r/ZZ);
	dy[3] = -0.5/D/tmpscal1/tmpscal1*((2.0*ZZ*y[3]*y[3]/AA*tmpscal1+mu/ZZ)*(tmpscal3+epsilon_l*y[3]*tmpscal4)+tmpscal5+epsilon_l*y[3]*tmpscal6+epsilon_g*AA/ZZ*(tmpvec1.z+epsilon_l*y[3]*tmpvec2.z))+Fdrag/fabs(ZZ)/tmpscal1;

        Fdrag=PI*rd0*rd0*ni1*1e20*(Vp1/v_th-y[3])*fabs(Vp1/v_th-y[3])*bmag;      
        dy[4]=3.0*kusai*tmag*dt/(rd1*speh*row);
        dy[5] = -kusai*tmag*dt/(row*lath);
           
        rd1=rd1+dy[5];

      
}
void rk4(BoutReal *y, int rkneq, BoutReal t, BoutReal dt, BoutReal *yout) {
	BoutReal dydx[rkneq],yt[rkneq],dyt[rkneq],dym[rkneq],dth,dt6,th;

	dth=dt*0.5;
	dt6=dt/6.0;
	th=t+dth;
        ti1 = TTN1.val(y[0],y[1],0,0,gmesh_fine);
        ni1 = TTN2.val(y[0],y[1],0,0,gmesh_fine);
        Vp1 = TTN3.val(y[0],y[1],0,0,gmesh_fine);
        pottl = phi0rz.val(y[0],y[1],0,0,gmesh_fine);
        
        te1=ti1;
        ne1=ni1;
	md=row*4.0*PI*rd1*rd1*rd1/3.0;
        kusai=1.0e9;
        vd=-2.5;
        AA=md/Mp;
        ZZ=vd*4.0*PI*e0*rd1*kb1*ti1/(qe*qe);


	gcf(t,y,rkneq,dydx);
	for(int i=0;i < rkneq;i++)
	{yt[i] = y[i]+dth*dydx[i];}
       
        if(yt[4]>tmelt){
        submd=yt[4]-tmelt;
        submd1=submd1+submd;;
        yt[4]=tmelt;
        }
        

	gcf(th,yt,rkneq,dyt);
	for(int i=0;i < rkneq;i++)
	{yt[i] = y[i]+dth*dyt[i];}
        
        if(yt[4]>tmelt){
        submd=yt[4]-tmelt;
        submd1=submd1+submd;;
        yt[4]=tmelt;
        }

	gcf(th,yt,rkneq,dym);
	for(int i=0;i < rkneq;i++)
	{yt[i]=y[i]+dt*dym[i];dym[i]=dyt[i]+dym[i];}
        
        if(yt[4]>tmelt){
        submd=yt[4]-tmelt;
        submd1=submd1+submd;;
        yt[4]=tmelt;
        }
     

	gcf(t+dt,yt,rkneq,dyt);
	for(int i=0;i < rkneq;i++)
	{yout[i]=y[i]+dt6*(dydx[i]+dyt[i]+2.0*dym[i]);}

        if(yout[4]>tmelt){
        submd=yout[4]-tmelt;
        submd1=submd1+submd;;
        yout[4]=tmelt;
        }

        submd11=lath/speh;
        test6=submd1;
        test7=submd11;
        if(submd1>submd11){
        output.write("Dust totally ablated. The remaining radius is ", y[5]);
        exit(0);}


}

/*BoutReal gaussian() {
	int static iset=0;
	BoutReal static res2;
	BoutReal rsq,v1,v2,fac;

	if(iset ==0) {
		rsq=1.1;
		while (rsq >= 1.0 or rsq==0.0) {
			v1 = rand()/(BoutReal)(RAND_MAX);
			v1=2.0*v1-1.0;
			v2 = rand()/(BoutReal)(RAND_MAX);
			v2=2.0*v2-1.0;
			rsq=v1*v1+v2*v2;
		}
		fac = sqrt(-2.0*log(rsq)/rsq);
		iset = 1;
		res2 = v2*fac;
		return v1*fac;
	}
	else {
		iset=0;
		return res2;
	}
}
*/



void checkpsi( GlobalField2D &gpsixy) {
	output.write("%3i%3i%5i\n",nx,ny,ixseps1);
	for(int i=0;i<ny;i++) {
		output.write("%3i%20.10e%20.10e\n",i,gpsixy(ixseps1-1,i),gpsixy(ixseps1,i));
	}
}

//generate the polygon, which describing the computation domain of MHD part.
void domain(GlobalField2D &gRxy, GlobalField2D &gZxy, int nv, BoutReal *vr, BoutReal *vz) {
	int i,npf1,npf2;
	npf1 = jyseps1_1+1;
	npf2 = ny-jyseps2_2-1;
	nv = ny+npf1+npf2;	//for single null case
	for(i=0;i<ny;i++) {
		vr[i] = gRxy(nx-1,i);
		vz[i] = gZxy(nx-1,i);
	}
	for(i=0;i<npf2;i++) {
		vr[ny+i] = gRxy(0,ny-1-i);
		vz[ny+i] = gZxy(0,ny-1-i);
	}
	for(i=0;i<npf1;i++) {
		vr[ny+npf2+i] = gRxy(0,jyseps1_1-i);
		vz[ny+npf2+i] = gZxy(0,jyseps1_1-i);
	}
}

void map_solcore(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek) {
	BoutReal rk,zk,det,drdpsi,drdthe,dzdpsi,dzdthe,psi0,the0,dr2,tmppsi,r,z,tmp,thetamin,thetamax,thetaperiod;
	int xindex,yindex;
	//psi0 = psik;
	//the0 = thek;
	//output.write("%20.10e%20.10e\n",rt,zt);
	thetamin = solcoremesh.gridz[0];
	thetamax = solcoremesh.gridz[solcoremesh.nz-1];
	thetaperiod = thetamax-thetamin;
	psi0 = psip.val(rt,zt,0,0,gmesh);
	psi0 = max(min(psi0,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);
	the0 = fmod(thek-thetamin,thetaperiod)+thetamin;
	the0 = max(min(the0,thetamax),thetamin);
	psik = psi0;
	for (int i=0; i<niterv;i++) {
		//make psik,thek falls inside the domian of grid!!!
		//output.write("%5i%20.10e%20.10e\n",i,psik,thek);
		rk = rsolcore.val(psik,thek,0,0,solcoremesh);
		zk = zsolcore.val(psik,thek,0,0,solcoremesh);
		
		dr2 = ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		if (  dr2 < epsilon_map) {
			xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}
     // output.write("%5i%20.10e%20.10e%20.10e%20.10e%20.10e\n",i,psik,thek,rk,zk,(rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		drdpsi = rsolcore.val(psik,thek,1,0,solcoremesh);
		drdthe = rsolcore.val(psik,thek,0,1,solcoremesh);
		dzdpsi = zsolcore.val(psik,thek,1,0,solcoremesh);
		dzdthe = zsolcore.val(psik,thek,0,1,solcoremesh);
		det = drdpsi*dzdthe - drdthe*dzdpsi;
		psik = psik + ( (rt-rk)*dzdthe-(zt-zk)*drdthe )/det;
		thek = thek + (-(rt-rk)*dzdpsi+(zt-zk)*drdpsi )/det;	

		psik = max(min(psik,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);
		thek = fmod(thek-thetamin,thetaperiod)+thetamin;
		thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);	
	}
	if (dr2 < 100.0*epsilon_map) {
		//psik = max(min(psik,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);	
		//thek = fmod(thek-thetamin,thetaperiod)+thetamin;
		//thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);
		xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
		yindex = int(thek/dtheta+0.5);
		//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
		density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
		return;
	}
	//psi0 = psip.val(rt,zt,0,0,gmesh);
	//psi0 = max(min(psi0,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);
	thek = the0;
	for (int i=0; i<niterv;i++) {
		rk = rsolcore.val(psi0,thek,0,0,solcoremesh);
		zk = zsolcore.val(psi0,thek,0,0,solcoremesh);	
		dr2 = ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		if (  dr2 < 100.0*epsilon_map) {
			xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}
		drdthe = rsolcore.val(psi0,thek,0,1,solcoremesh);
		if(abs(drdthe)<epsilon_map) {
			continue;
		}
		r = rsolcore.val(psi0,thek,0,0,solcoremesh);
		tmp = r / drdthe;

		thek = thek - tmp;
		thek = fmod(thek-thetamin,thetaperiod)+thetamin;	
		thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);
	}
	thek = the0;
	for (int i=0; i<niterv;i++) {
		rk = rsolcore.val(psi0,thek,0,0,solcoremesh);
		zk = zsolcore.val(psi0,thek,0,0,solcoremesh);
		
		dr2 = ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		if (  dr2 < 100.0*epsilon_map) {
			xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}
		dzdthe = zsolcore.val(psi0,thek,0,1,solcoremesh);
		if(abs(dzdthe)<epsilon_map) {
			continue;
		}
		z = zsolcore.val(psi0,thek,0,0,solcoremesh);
		tmp = z / dzdthe;
		thek = thek - tmp;
		thek = fmod(thek-thetamin,thetaperiod)+thetamin;
		thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);
	}

	psik = psi0;
	thek = PI;
////////////////////////
	for (int i=0; i<niterv;i++) {
		//make psik,thek falls inside the domian of grid!!!
		//output.write("%5i%20.10e%20.10e\n",i,psik,thek);

		rk = rsolcore.val(psik,thek,0,0,solcoremesh);
		zk = zsolcore.val(psik,thek,0,0,solcoremesh);
		
		dr2 = ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		if (  dr2 < epsilon_map) {
			xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}
     // output.write("%5i%20.10e%20.10e%20.10e%20.10e%20.10e\n",i,psik,thek,rk,zk,(rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		drdpsi = rsolcore.val(psik,thek,1,0,solcoremesh);
		drdthe = rsolcore.val(psik,thek,0,1,solcoremesh);
		dzdpsi = zsolcore.val(psik,thek,1,0,solcoremesh);
		dzdthe = zsolcore.val(psik,thek,0,1,solcoremesh);
		det = drdpsi*dzdthe - drdthe*dzdpsi;
		psik = psik + ( (rt-rk)*dzdthe-(zt-zk)*drdthe )/det;
		thek = thek + (-(rt-rk)*dzdpsi+(zt-zk)*drdpsi )/det;	

		psik = max(min(psik,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);
		thek = fmod(thek-thetamin,thetaperiod)+thetamin;
		thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);	
	}
	if (dr2 < 100.0*epsilon_map) {
		//psik = max(min(psik,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);	
		//thek = fmod(thek-thetamin,thetaperiod)+thetamin;
		//thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);
		xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
		yindex = int(thek/dtheta+0.5);
		//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
		density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
		return;
	}
	//psi0 = psip.val(rt,zt,0,0,gmesh);
	//psi0 = max(min(psi0,solcoremesh.gridr[solcoremesh.nr-1]),solcoremesh.gridr[0]);
	thek = the0;
	for (int i=0; i<niterv;i++) {
		rk = rsolcore.val(psi0,thek,0,0,solcoremesh);
		zk = zsolcore.val(psi0,thek,0,0,solcoremesh);	
		dr2 = ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		if (  dr2 < 100.0*epsilon_map) {
			xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}
		drdthe = rsolcore.val(psi0,thek,0,1,solcoremesh);
		if(abs(drdthe)<epsilon_map) {
			continue;
		}
		r = rsolcore.val(psi0,thek,0,0,solcoremesh);
		tmp = r / drdthe;

		thek = thek - tmp;
		thek = fmod(thek-thetamin,thetaperiod)+thetamin;	
		thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);
	}
	thek = the0;
	for (int i=0; i<niterv;i++) {
		rk = rsolcore.val(psi0,thek,0,0,solcoremesh);
		zk = zsolcore.val(psi0,thek,0,0,solcoremesh);
		
		dr2 = ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt));
		if (  dr2 < 100.0*epsilon_map) {
			xindex = cdblocate(solcoremesh.gridr, solcoremesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}
		dzdthe = zsolcore.val(psi0,thek,0,1,solcoremesh);
		if(abs(dzdthe)<epsilon_map) {
			continue;
		}
		z = zsolcore.val(psi0,thek,0,0,solcoremesh);
		tmp = z / dzdthe;
		thek = thek - tmp;
		thek = fmod(thek-thetamin,thetaperiod)+thetamin;
		thek = max(min(thek,solcoremesh.gridz[solcoremesh.nz-1]),solcoremesh.gridz[0]);
	}

//////////////
	merr_solcore = merr_solcore+1;
	output.write("merr_solcore is: %5i\n",merr_solcore);
	tmppsi = psi0;
	tmppsi = (tmppsi-psia)/(psib-psia);
	err_psi = err_psi+tmppsi;
	//output.write("map_solcore err: %20.10e%20.10e%20.10e%20.10e%20.10e\n",rt,zt,tmppsi,the0,thek);
	output.write("%20.10e%20.10e%5i\n",rt,zt,0);
}

void map_pf1(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek) {
	BoutReal rk,zk,det,drdpsi,drdthe,dzdpsi,dzdthe;
	int xindex,yindex;
	for (int i=0; i<niterv;i++) {
		//make psik,thek falls inside the domian of grid!!!
		psik = max(min(psik,pf1mesh.gridr[pf1mesh.nr-1]),pf1mesh.gridr[0]);	
		thek = max(min(thek,pf1mesh.gridz[pf1mesh.nz-1]),pf1mesh.gridz[0]);
		rk = rpf1.val(psik,thek,0,0,pf1mesh);
		zk = zpf1.val(psik,thek,0,0,pf1mesh);
		if ( ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt)) < epsilon_map) {
			xindex = cdblocate(pf1mesh.gridr, pf1mesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("PF1. xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}

		drdpsi = rpf1.val(psik,thek,1,0,pf1mesh);
		drdthe = rpf1.val(psik,thek,0,1,pf1mesh);
		dzdpsi = zpf1.val(psik,thek,1,0,pf1mesh);
		dzdthe = zpf1.val(psik,thek,0,1,pf1mesh);
		det = drdpsi*dzdthe - drdthe*dzdpsi;
		psik = psik + ( (rt-rk)*dzdthe-(zt-zk)*drdthe )/det;
		thek = thek + (-(rt-rk)*dzdpsi+(zt-zk)*drdpsi )/det;
	}
	merr_pf1 = merr_pf1+1;
	//output.write("map_pf1 err: %20.10e%20.10e%20.10e%20.10e\n",rt,zt,psik,thek);
	output.write("merr_pf1 is: %5i\n",merr_pf1);
	output.write("%20.10e%20.10e%5i\n",rt,zt,1);
}

void map_pf2(BoutReal rt, BoutReal zt, BoutReal &psik, BoutReal &thek) {
	BoutReal rk,zk,det,drdpsi,drdthe,dzdpsi,dzdthe;
	int xindex,yindex;
	thek= 0.5*(pf2mesh.gridz[pf2mesh.nz-1] + pf2mesh.gridz[0]);
	for (int i=0; i<niterv;i++) {
		//make psik,thek falls inside the domian of grid!!!
		psik = max(min(psik,pf2mesh.gridr[pf2mesh.nr-1]),pf2mesh.gridr[0]);	
		thek = max(min(thek,pf2mesh.gridz[pf2mesh.nz-1]),pf2mesh.gridz[0]);
		rk = rpf2.val(psik,thek,0,0,pf2mesh);
		zk = zpf2.val(psik,thek,0,0,pf2mesh);
		if ( ((rk-rt)*(rk-rt)+(zk-zt)*(zk-zt)) < epsilon_map) {
			xindex = cdblocate(pf2mesh.gridr, pf2mesh.nr, psik);
			yindex = int(thek/dtheta+0.5);
			//(*gndes)(xindex,yindex) = (*gndes)(xindex,yindex)+1;
			density[xindex][yindex] = density[xindex][yindex] +1;
			if(xindex>nx-1 or xindex<0) {output.write("pf2.xindex,rt,zt,psik are: %d,%e,%e,%e\n",xindex,rt,zt,psik);}
			if(yindex>ny-1 or yindex<0) {output.write("yindex is: %d\n",yindex);}
			return;
		}

		drdpsi = rpf2.val(psik,thek,1,0,pf2mesh);
		drdthe = rpf2.val(psik,thek,0,1,pf2mesh);
		dzdpsi = zpf2.val(psik,thek,1,0,pf2mesh);
		dzdthe = zpf2.val(psik,thek,0,1,pf2mesh);
		det = drdpsi*dzdthe - drdthe*dzdpsi;
		psik = psik + ( (rt-rk)*dzdthe-(zt-zk)*drdthe )/det;
		thek = thek + (-(rt-rk)*dzdpsi+(zt-zk)*drdpsi )/det;
	}
	merr_pf2 = merr_pf2+1;
	output.write("merr_pf2 is: %5i\n",merr_pf2);
	//output.write("map_pf2 err: %20.10e%20.10e%20.10e%20.10e\n",rt,zt,psik,thek);
	output.write("%20.10e%20.10e%5i\n",rt,zt,2);
}


//map from RZ to (psi,theta) in whole MHD mesh.
//pnpoly=1, inside polygon; pnpoly=0, outside polygon.
//inner most region:-1, sol+core:0, PF1:1, PF2:2, blank region:5, outside whole domain: 9
void maprz2pt(particle &impurity) {
	int xindex,yindex;
	BoutReal theta0=impurity.theta;
	if(impurity.region==9||impurity.region==-1) 	return;
	//if(impurity.region==9) 	return;
	if ( pnpoly(ndomain, r_domain, z_domain, impurity.y[0],impurity.y[1]) == 1){
		if( pnpoly(nsolcore, r_solcore, z_solcore, impurity.y[0],impurity.y[1]) == 1) {
			impurity.region = 0;
			map_solcore(impurity.y[0],impurity.y[1],impurity.psi,impurity.theta);
		}
		else if(pnpoly(npf1, r_pf1, z_pf1, impurity.y[0],impurity.y[1]) == 1){
			impurity.region = 1;
			map_pf1(impurity.y[0],impurity.y[1],impurity.psi,impurity.theta);
		}
		else if(pnpoly(npf2, r_pf2, z_pf2, impurity.y[0],impurity.y[1]) == 1) {
			impurity.region = 2;
			map_pf2(impurity.y[0],impurity.y[1],impurity.psi,impurity.theta);
		}
		else if(pnpoly(ncore, r_core, z_core, impurity.y[0],impurity.y[1]) == 1) {
			impurity.region = -1;
			tmp_core = tmp_core+1;
		}
		else {
			impurity.region = 5;
			tmp_blank = tmp_blank+1;
		}	
	}
	else {
		yindex=int(theta0/dtheta+0.5);
		if(impurity.region==1) {
			lost_pf1=lost_pf1+1;
		}
		else if(impurity.region==2) {
			lost_pf2=lost_pf2+1;
		}
		else if(impurity.region==5) {
			lost_blank=lost_blank+1;
		}
		else if(impurity.region==0) {
			//sol is divided to 4 sections: [4,[22),[30),[39),59] in g38300_132X64_0.9-1.02.nc configuration.
			//sol is divided to 4 sections: [8,[23),[30),[37),55] in g056129.05550_2_x260y64_psi070to106_fitTe_decreasePF_gaussianY.nc configuration.
			if(yindex >= 8 and yindex<23) {
				lost_sol1=lost_sol1+1;
			}
			else if(yindex >= 23 and yindex<30) {
				lost_sol2=lost_sol2+1;
			}
			else if(yindex >= 30 and yindex<37) {
				lost_sol3=lost_sol3+1;
			}
			else if(yindex >= 37 and yindex<=55) {
				lost_sol4=lost_sol4+1;
			}
			else {
				lost_sol5=lost_sol5+1;
				//output.write("lost_sol5 errs, the angle and yindex are: %f,%d\n",theta0,yindex);
			}
		}
		else {
			lost_other=lost_other+1;
		}
		impurity.region = 9;
		lost_total = lost_total+1;
	}
}

BoutReal gaussian() {
	int static iset=0;
	BoutReal static res2;
	BoutReal rsq,v1,v2,fac;

	if(iset ==0) {
		rsq=1.1;
		while (rsq >= 1.0 or rsq==0.0) {
			v1 = rand()/(BoutReal)(RAND_MAX);
			v1=2.0*v1-1.0;
			v2 = rand()/(BoutReal)(RAND_MAX);
			v2=2.0*v2-1.0;
			rsq=v1*v1+v2*v2;
		}
		fac = sqrt(-2.0*log(rsq)/rsq);
		iset = 1;
		res2 = v2*fac;
		return v1*fac;
	}
	else {
		iset=0;
		return res2;
	}
}


void randomwalk(BoutReal dr, BoutReal dz, BoutReal &dr_tb, BoutReal &dz_tb) {
	BoutReal  r,dl,tmp;

	r=sqrt(dr*dr+dz*dz);
	//dl = gaussian()+dl_tb;
	tmp = gaussian();
	dl = tmp*dl_tb;
	

	dr_tb=dz/r*dl;
	dz_tb=-dr/r*dl;
//output.write("random of gaussian(0,1) is: %20.10e%20.10e20.10e%20.10e20.10e%20.10e\n", tmp,dl/r,dr_tb,dz_tb);
}
