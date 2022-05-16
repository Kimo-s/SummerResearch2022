/*---------------------------------------------------------------------------------------
  Simulation of electron-ion recombination in nuclear recoils tracks in LAr at 87 K.
  Variable time step. Space-energy criterion of reaction.
  Mobile cations included.
  Cross sections available up to 100 eV.
  An energy distribution of secondary electrons included.
  External electric field applied.
  MAGNETIC FIELD APPLIED.
                              by Mariusz Wojcik (Lodz Univ. of Technology), Dec. 2, 2019
----------------------------------------------------------------------------------------*/
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ran_gen.h"
#include "LAr_crsec.h"
#include "LAr_fEk_sec_el.h"

#define	GO_ON	0
#define	FIN		1

#define	EL_CHG	1.6022e-19					//e, C
#define	EL_MASS	9.1095e-31					//m, kg
#define	N_AV	6.022e23					//Avogadro number
#define M_AR	39.95						//molar mass of Ar
#define	EPS_0	8.8542e-12					//epsilon0, [SI]
#define	KB		1.3807e-23					//Boltzmann constant, J/K
#define	PI		3.14159265358979323846		//pi
#define	TWO_PI	(2.0*PI)

#define	TEMP	87.0			//T, K
#define	DENSITY	2.11e28			//n, m^(-3)
#define	EPS		1.49			//dielectric constant
#define	R_CAT	2.0				//radius of the cation, A
#define	R_CAT_0	1.9				//radius used by the potential smoothing function, A
#define	R_MIN_CAT	0.1			//minimum initial distance between the cations, A
#define	EK0_EV_MAX	100.		//max. Ek0 [eV]

#define	IN_FNAME	"LAr_dm_B.inp"	//input filename
#define	RES_FNAME	"LAr_dm_B.res"	//results filename
#define	DET_FNAME	"LAr_dm_B.det"	//detailed results filename
#define	OUT_FNAME	"LAr_dm_B.out"	//screen output filename

#define	NP	128	//maximum number of cation-electron pairs
#define	N_EL_CAT	(NP*NP)
#define	N_RE	20000	//max. size of the detailed output file
#define	NTAB	100		//max. no. of recorded Prec values

#define	NDT		33			//number of time step values
#define	DT_INIT	19			//initial time step number
double dts_fs[NDT]={1.e-7,1.8e-7,3.2e-7,5.6e-7,1.e-6,1.8e-6,3.2e-6,5.6e-6,
                    1.e-5,1.8e-5,3.2e-5,5.6e-5,1.e-4,1.8e-4,3.2e-4,5.6e-4,
					1.e-3,1.8e-3,3.2e-3,5.6e-3,0.01,0.018,0.032,0.056,
					0.1,0.18,0.32,0.56,1.,1.8,3.2,5.6,10.};//time steps, fs

int Npair;			//initial number of cation-electron pairs
double rk_A;		//cation-cation distance, A
int rk_dist;		//cation-cation distances: equal (rk_dist=0) or exponentially distributed (rk_dist=1)
double cat_mob_cm2Vs;	//cation mobility, cm2/Vs
double cat_mass_MAr;	//cation mass (multiple of argon atom mass)
double F_Vcm;		//external electric field, V/cm
double Fd_x;		//direction of F - x component
double Fd_y;		//direction of F - y component
double Fd_z;		//direction of F - z component
double B_T;			//magnetic field, T
double Bd_x;		//direction of B - x component
double Bd_y;		//direction of B - y component
double Bd_z;		//direction of B - z component
int Ek0_dist;		//is Ek0 distributed (0-No, 1-Yes)
double r0_A;		//initial electron-cation distance, A
double Ek0_eV;		//initial kinetic energy of electrons, eV (not used if Ek0_dist=1)
double Rcrit_A;		//critical distance for the reaction criterion, A
double Ekcrit_eV;	//critical kinetic energy for the reaction criterion, eV
double Rmax_A;		//critical distance for the escape criterion, A
long Nrun;			//number of runs in the simulation
long Nav;			//number of runs between sub-averaging

// Time step decrease/increase thresholds:
double g_E,h_E;		// energy drift
double g_tau,h_tau;	// dt/tau

double rk,r0,Ek0,Rcrit,Ekcrit,Rmax,Rmax_sq;
double F,Fx,Fy,Fz,eFx,emFx,eMcatFx,eFy,emFy,eMcatFy,eFz,emFz,eMcatFz;
double B,Bx,By,Bz,emBx,emBy,emBz;
double dens_inv;
double Rcat,Rcat0,Rcat1;
double Rmin_cat;
double t,dt,dt_half,dt_sq_half;
double dts[NDT],dts_half[NDT],dts_sq_half[NDT],r_tau_cats[NDT],P_nccs[NDT];
int idt,dt_warn;
int ndt_max_cat;
int Ncat,Nel,Ncatm1,Nelm1;
double xcat[NP],ycat[NP],zcat[NP],xcat_n[NP],ycat_n[NP],zcat_n[NP];
int cat_id[NP];
double vcatx[NP],vcaty[NP],vcatz[NP],vcatx_n[NP],vcaty_n[NP],vcatz_n[NP];
double acatx[NP],acaty[NP],acatz[NP],acatx_n[NP],acaty_n[NP],acatz_n[NP];
double x[NP],y[NP],z[NP],x_n[NP],y_n[NP],z_n[NP];
int el_id[NP];
double vx[NP],vy[NP],vz[NP],vx_n[NP],vy_n[NP],vz_n[NP];
double ax[NP],ay[NP],az[NP],ax_n[NP],ay_n[NP],az_n[NP];
double r_el_cat[N_EL_CAT],r_el_cat_n[N_EL_CAT];
double Ep,Ekelt,Ekcatt,Et,Ekel_av,Ep_n,Ekelt_n,Ekcatt_n,Et_n,Ekel_n_av,EpRcat;
double Ekel_r;
double Ek[NP],Ek_n[NP];
double Ekcat[NP],Ekcat_n[NP];
double Acat,Bcat,Ccat,Amcat;
double E_Coul_r0;
double vMx[NP],vMy[NP],vMz[NP];
double ux[NP],uy[NP],uz[NP],u[NP];
double sigma_E[NP],tau[NP],r_tau[NP];
int coll[NP],coll_E[NP],coll_cat[NP];
double tau_cat;
double r_tau_cat,P_ncoll_cat;
double cat_mob;
double M,m_half,Mcat,Mcat_half;
double k,two_k,km,kMcat;
double kBT,kBTM,kBTMcat;
double msm,Msm;
double mMcat;
int kkmax;
long irun;
long currentSeed;
int Nrec,Nout,Ngem;
int i_rm;
double t_rm[NP],Ek_rm[NP],x_rm[NP],y_rm[NP],z_rm[NP];
char type_rm[NP];
int id_rm[NP];
long Nrec_tot,Nout_tot,Ngem_tot;
int i_re,run_re[N_RE];
double t_re[N_RE],Ek_re[N_RE],x_re[N_RE],y_re[N_RE],z_re[N_RE];
char type_re[N_RE];
int id_re[N_RE];
double Prec,Prec_tab[NTAB];
double Pgem,Pgem_tab[NTAB];
int n_av;
double sum_Prec,sum2_Prec,Prec_av,Prec_stderr;
double sum_Pgem,sum2_Pgem,Pgem_av,Pgem_stderr;
FILE *fp_out;

void traj(void);
void count_up(void);
void init_el(void);
int check_traj(void);
void dt_set(void);
void dt_dec(void);
void dt_inc(void);
void input(void);
void output(int);
void err_exit(char[]);


void main(int argc, char *argv[])
{
	char *eptr;
	if(argc == 2){
		currentSeed = strtol(argv[1],&eptr,10);
	} else {
		srand(time(0));
		currentSeed = (long)rand();
	}
	ran_init(currentSeed); //initialize the random number generator (with the seed of long integer type)
	input();		//read input data
	output(1);		//initialize output procedures

	
	printf("Seed supplied from command line: %ld\n", currentSeed);
	//	ran_init(528645);
	crsec_init();	//initialize the cross section data
	if (Ek0_dist)
		fEk_sec_el_init();	//initialize the energy distribution data

	for (irun=1;irun<=Nrun;irun++) {			//main loop
		printf("run no. %d\t",irun);
		fprintf(fp_out,"run no. %d\t",irun);

		traj();		//calculate a complete trajectory for a cluster of electron-ion pairs

		printf("\n");
		fprintf(fp_out,"\n");
		if (irun%Nav == 0)		//periodic sub-averaging and output
			output(2);
	}
	output(3);		//final output
}


// Initialization of a cluster of cation-electron pairs
void init_cat_el()
{
	int i,j;
	double ctheta,stheta,phi;
	double v0;
	double zz,xav,yav,zav;
	double rkk;
	double xcat_ij,ycat_ij,zcat_ij,rcat_ij;
	int rcat_short;
	double f_Ek,r_accpt;

	Ncat=Npair;
	Nel=Npair;
	Ncatm1=Ncat-1;
	Nelm1=Nel-1;

	if (rk_dist==0) {				//linear cluster, equal cation-cation distances
		xcat[0]=0.;
		ycat[0]=0.;
		zcat[0]=-0.5*(Ncat-1)*rk;
		cat_id[0]=0;
		for (i=1;i<Ncat;i++) {			
			xcat[i]=0.;
			ycat[i]=0.;
			zcat[i]=zcat[i-1]+rk;
			cat_id[i]=i;
		}

	} else if (rk_dist == 1) {		//linear cluster, exponentially distributed cation-cation distances
		zz=zav=0.;
		for (i=0;i<Ncat;i++) {
			xcat[i]=0.;
			ycat[i]=0.;
			zcat[i]=zz;
			cat_id[i]=i;
			zav+=zz;
			do {
				rkk=ran_exp()*rk;
			} while (rkk < Rmin_cat);
			zz+=rkk;
		}
		zav/=Ncat;
		for (i=0;i<Ncat;i++)
			zcat[i]-=zav;

	} else if (rk_dist == 5) {		//random cluster, equal cation-cation distances
		xcat[0]=0.;
		ycat[0]=0.;
		zcat[0]=0.;
		cat_id[0]=0;
		xav=yav=zav=0.;
		for (i=1;i<Ncat;i++) {
			do {
				ctheta=ran()*2.-1.;
				stheta=sqrt(1.-ctheta*ctheta);
				phi=ran()*TWO_PI;
				xcat[i]=xcat[i-1]+rk*stheta*cos(phi);
				ycat[i]=ycat[i-1]+rk*stheta*sin(phi);
				zcat[i]=zcat[i-1]+rk*ctheta;
				rcat_short=0;
				for (j=0;j<i;j++) {
					xcat_ij=xcat[j]-xcat[i];
					ycat_ij=ycat[j]-ycat[i];
					zcat_ij=zcat[j]-zcat[i];
					rcat_ij=sqrt(xcat_ij*xcat_ij+ycat_ij*ycat_ij+zcat_ij*zcat_ij);
					if (rcat_ij < Rmin_cat)
						rcat_short=1;
				}
			} while (rcat_short);
			cat_id[i]=i;
			xav+=xcat[i];
			yav+=ycat[i];
			zav+=zcat[i];
		}
		xav/=Ncat;
		yav/=Ncat;
		zav/=Ncat;
		for (i=0;i<Ncat;i++) {
			xcat[i]-=xav;
			ycat[i]-=yav;
			zcat[i]-=zav;
		}

	} else {						//random cluster, exponentially distributed cation-cation distances
		xcat[0]=0.;
		ycat[0]=0.;
		zcat[0]=0.;
		cat_id[0]=0;
		xav=yav=zav=0.;
		for (i=1;i<Ncat;i++) {
			do {
				ctheta=ran()*2.-1.;
				stheta=sqrt(1.-ctheta*ctheta);
				phi=ran()*TWO_PI;
				rkk=ran_exp()*rk;
				xcat[i]=xcat[i-1]+rkk*stheta*cos(phi);
				ycat[i]=ycat[i-1]+rkk*stheta*sin(phi);
				zcat[i]=zcat[i-1]+rkk*ctheta;
				rcat_short=0;
				for (j=0;j<i;j++) {
					xcat_ij=xcat[j]-xcat[i];
					ycat_ij=ycat[j]-ycat[i];
					zcat_ij=zcat[j]-zcat[i];
					rcat_ij=sqrt(xcat_ij*xcat_ij+ycat_ij*ycat_ij+zcat_ij*zcat_ij);
					if (rcat_ij < Rmin_cat)
						rcat_short=1;
				}
			} while (rcat_short);
			cat_id[i]=i;
			xav+=xcat[i];
			yav+=ycat[i];
			zav+=zcat[i];
		}
		xav/=Ncat;
		yav/=Ncat;
		zav/=Ncat;
		for (i=0;i<Ncat;i++) {
			xcat[i]-=xav;
			ycat[i]-=yav;
			zcat[i]-=zav;
		}
	}
	for (i=0;i<Ncat;i++) {				//initial cation velocities
		vcatx[i]=kBTMcat*ran_gauss();
		vcaty[i]=kBTMcat*ran_gauss();
		vcatz[i]=kBTMcat*ran_gauss();
	}
	
	for (i=0;i<Nel;i++) {				//initial electron positions
		ctheta=ran()*2.-1.;				
		stheta=sqrt(1.-ctheta*ctheta);
		phi=ran()*TWO_PI;
		x[i]=xcat[i]+r0*stheta*cos(phi);
		y[i]=ycat[i]+r0*stheta*sin(phi);
		z[i]=zcat[i]+r0*ctheta;
		el_id[i]=i;
	}

	if (Ek0_dist == 0) {
		v0=sqrt(2.*Ek0/EL_MASS);
		for (i=0;i<Nel;i++) {				//initial electron velocities
			ctheta=ran()*2.-1.;
			stheta=sqrt(1.-ctheta*ctheta);
			phi=ran()*TWO_PI;
			vx[i]=v0*stheta*cos(phi);
			vy[i]=v0*stheta*sin(phi);
			vz[i]=v0*ctheta;
		}
	} else {
		for (i=0;i<Nel;i++) {
			do {
				Ek0_eV=ran()*EK0_EV_MAX;
				f_Ek=fEk(Ek0_eV);
				r_accpt=ran();
			} while (r_accpt > f_Ek);
			Ek0=Ek0_eV*EL_CHG;
			Ek0-=E_Coul_r0;
			v0=sqrt(2.*Ek0/EL_MASS);
			ctheta=ran()*2.-1.;
			stheta=sqrt(1.-ctheta*ctheta);
			phi=ran()*TWO_PI;
			vx[i]=v0*stheta*cos(phi);
			vy[i]=v0*stheta*sin(phi);
			vz[i]=v0*ctheta;
		}
	}
	Nrec=0;
	Nout=0;
	Ngem=0;
	i_rm=0;
}


// Calculates a complete trajectory for a cluster of electron-cation pairs
void traj()
{
	int i,j,kk;
	double xij,yij,zij,rij_sq,rij,axij,ayij,azij;
	double kmr3,Amcatr;
	double kMcatr3,AMccatr;
	int fate;
	int dt_high,dt_increase;
	double r_E;
	int colls,colls_cat;
	double P_ncoll,R_coll;
	double sigma_p;
	double r_p,R_p;
	double ctheta,stheta,phi,nx,ny,nz;
	double Msmu;
	double vxc,vyc,vzc,r_vc_v;

	init_cat_el();		//initialize the cations and electrons
	dt_set();			//set the time step
	t=0.;

// Calculate the accelerations and energies at t
	Ep=0.;
	for (i=0;i<Nel;i++) ax[i]=ay[i]=az[i]=0.;
	for (j=0;j<Ncat;j++) acatx[j]=acaty[j]=acatz[j]=0.;

	kk=0;
	for (i=0;i<Nel;i++) {			//electron-cation contributions
		for (j=0;j<Ncat;j++) {
			xij=x[i]-xcat[j];
			yij=y[i]-ycat[j];
			zij=z[i]-zcat[j];
			rij_sq=xij*xij+yij*yij+zij*zij;
			r_el_cat[kk++]=rij=sqrt(rij_sq);
			if (rij > Rcat1) {
				kmr3=-km/(rij*rij_sq);
				ax[i]+=kmr3*xij;
				ay[i]+=kmr3*yij;
				az[i]+=kmr3*zij;
				kMcatr3=-kmr3*mMcat;
				acatx[j]+=kMcatr3*xij;
				acaty[j]+=kMcatr3*yij;
				acatz[j]+=kMcatr3*zij;
				Ep-=k/rij;
			} else if (rij > Rcat0) {
				Amcatr=Amcat*(Rcat0/rij-1.0);
				ax[i]+=Amcatr*xij;
				ay[i]+=Amcatr*yij;
				az[i]+=Amcatr*zij;
				AMccatr=-Amcatr*mMcat;
				acatx[j]+=AMccatr*xij;
				acaty[j]+=AMccatr*yij;
				acatz[j]+=AMccatr*zij;
				Ep+=Acat*rij_sq+Bcat*rij+Ccat;
			} else
				Ep+=EpRcat;
		}
	}
	kkmax=kk;

	for (i=0;i<Nelm1;i++) {			//electron-electron contributions
		for (j=i+1;j<Nel;j++) {
			xij=x[i]-x[j];
			yij=y[i]-y[j];
			zij=z[i]-z[j];
			rij_sq=xij*xij+yij*yij+zij*zij;
			rij=sqrt(rij_sq);
			kmr3=km/(rij*rij_sq);
			axij=kmr3*xij;
			ayij=kmr3*yij;
			azij=kmr3*zij;
			ax[i]+=axij;
			ay[i]+=ayij;
			az[i]+=azij;
			ax[j]-=axij;
			ay[j]-=ayij;
			az[j]-=azij;
			Ep+=k/rij;
		}
	}

	for (i=0;i<Ncatm1;i++) {			//cation-cation contributions
		for (j=i+1;j<Ncat;j++) {
			xij=xcat[i]-xcat[j];
			yij=ycat[i]-ycat[j];
			zij=zcat[i]-zcat[j];
			rij_sq=xij*xij+yij*yij+zij*zij;
			rij=sqrt(rij_sq);
			kMcatr3=kMcat/(rij*rij_sq);
			axij=kMcatr3*xij;
			ayij=kMcatr3*yij;
			azij=kMcatr3*zij;
			acatx[i]+=axij;
			acaty[i]+=ayij;
			acatz[i]+=azij;
			acatx[j]-=axij;
			acaty[j]-=ayij;
			acatz[j]-=azij;
			Ep+=k/rij;
		}
	}

	for (i=0;i<Nel;i++) {			//external electric field contributions
		ax[i]-=emFx;
		ay[i]-=emFy;
		az[i]-=emFz;
		Ep+=eFx*x[i]+eFy*y[i]+eFz*z[i];
	}
	for (j=0;j<Ncat;j++) {
		acatx[j]+=eMcatFx;
		acaty[j]+=eMcatFy;
		acatz[j]+=eMcatFz;
		Ep-=eFx*xcat[j]+eFy*ycat[j]+eFz*zcat[j];
	}

	for (i=0;i<Nel;i++) {			//magnetic field contributions (only for electrons)
		ax[i]-=emBz*vy[i]-emBy*vz[i];
		ay[i]-=emBx*vz[i]-emBz*vx[i];
		az[i]-=emBy*vx[i]-emBx*vy[i];
	}

	Ekelt=0.;
	for (i=0;i<Nel;i++) {							//kinetic energies for electrons
		Ek[i]=m_half*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
		Ekelt+=Ek[i];
	}
	Ekcatt=0.;
	for (j=0;j<Ncat;j++) {							//kinetic energies for cations
		Ekcat[j]=Mcat_half*(vcatx[j]*vcatx[j]+vcaty[j]*vcaty[j]+vcatz[j]*vcatz[j]);
		Ekcatt+=Ekcat[j];
	}
	Et=Ep+Ekelt+Ekcatt;
	Ekel_av=Ekelt/Nel;

	while ((fate=check_traj())==GO_ON) {	//check for reaction or escape

		for (i=0;i<Nel;i++) {
			vMx[i]=kBTM*ran_gauss();		//generate the velocities of atoms
			vMy[i]=kBTM*ran_gauss();
			vMz[i]=kBTM*ran_gauss();
		}

		do {
			dt_high=0;
			dt_increase=1;

			for (i=0;i<Nel;i++) {
				x_n[i]=x[i]+vx[i]*dt+ax[i]*dt_sq_half;	//electron positions at t+dt
				y_n[i]=y[i]+vy[i]*dt+ay[i]*dt_sq_half;
				z_n[i]=z[i]+vz[i]*dt+az[i]*dt_sq_half;
			}
			for (j=0;j<Ncat;j++) {
				xcat_n[j]=xcat[j]+vcatx[j]*dt+acatx[j]*dt_sq_half;	//cation positions at t+dt
				ycat_n[j]=ycat[j]+vcaty[j]*dt+acaty[j]*dt_sq_half;
				zcat_n[j]=zcat[j]+vcatz[j]*dt+acatz[j]*dt_sq_half;
			}
			for (i=0;i<Nel;i++) {
				vx_n[i]=vx[i]+ax[i]*dt;					//predicted electron velocities at t+dt
				vy_n[i]=vy[i]+ay[i]*dt;
				vz_n[i]=vz[i]+az[i]*dt;
			}
/*
			for (j=0;j<Ncat;j++) {
				vcatx_n[j]=vcatx[j]+acatx[j]*dt;	//predicted cation velocities at t+dt (not needed)
				vcaty_n[j]=vcaty[j]+acaty[j]*dt;
				vcatz_n[j]=vcatz[j]+acatz[j]*dt;
			}
*/			
			Ep_n=0.;		//accelerations and potential energy at t+dt
			for (i=0;i<Nel;i++) ax_n[i]=ay_n[i]=az_n[i]=0.;
			for (j=0;j<Ncat;j++) acatx_n[j]=acaty_n[j]=acatz_n[j]=0.;

			kk=0;
			for (i=0;i<Nel;i++) {			//electron-cation contributions
				for (j=0;j<Ncat;j++) {
					xij=x_n[i]-xcat_n[j];
					yij=y_n[i]-ycat_n[j];
					zij=z_n[i]-zcat_n[j];
					rij_sq=xij*xij+yij*yij+zij*zij;
					r_el_cat_n[kk++]=rij=sqrt(rij_sq);
					if (rij > Rcat1) {
						kmr3=-km/(rij*rij_sq);
						ax_n[i]+=kmr3*xij;
						ay_n[i]+=kmr3*yij;
						az_n[i]+=kmr3*zij;
						kMcatr3=-kmr3*mMcat;
						acatx_n[j]+=kMcatr3*xij;
						acaty_n[j]+=kMcatr3*yij;
						acatz_n[j]+=kMcatr3*zij;
						Ep_n-=k/rij;
					} else if (rij > Rcat0) {
						Amcatr=Amcat*(Rcat0/rij-1.0);
						ax_n[i]+=Amcatr*xij;
						ay_n[i]+=Amcatr*yij;
						az_n[i]+=Amcatr*zij;
						AMccatr=-Amcatr*mMcat;
						acatx_n[j]+=AMccatr*xij;
						acaty_n[j]+=AMccatr*yij;
						acatz_n[j]+=AMccatr*zij;
						Ep_n+=Acat*rij_sq+Bcat*rij+Ccat;
					} else
						Ep_n+=EpRcat;
				}
			}
			kkmax=kk;

			for (i=0;i<Nelm1;i++) {			//electron-electron contributions
				for (j=i+1;j<Nel;j++) {
					xij=x_n[i]-x_n[j];
					yij=y_n[i]-y_n[j];
					zij=z_n[i]-z_n[j];
					rij_sq=xij*xij+yij*yij+zij*zij;
					rij=sqrt(rij_sq);
					kmr3=km/(rij*rij_sq);
					axij=kmr3*xij;
					ayij=kmr3*yij;
					azij=kmr3*zij;
					ax_n[i]+=axij;
					ay_n[i]+=ayij;
					az_n[i]+=azij;
					ax_n[j]-=axij;
					ay_n[j]-=ayij;
					az_n[j]-=azij;
					Ep_n+=k/rij;
				}
			}

			for (i=0;i<Ncatm1;i++) {			//cation-cation contributions
				for (j=i+1;j<Ncat;j++) {
					xij=xcat_n[i]-xcat_n[j];
					yij=ycat_n[i]-ycat_n[j];
					zij=zcat_n[i]-zcat_n[j];
					rij_sq=xij*xij+yij*yij+zij*zij;
					rij=sqrt(rij_sq);
					kMcatr3=kMcat/(rij*rij_sq);
					axij=kMcatr3*xij;
					ayij=kMcatr3*yij;
					azij=kMcatr3*zij;
					acatx_n[i]+=axij;
					acaty_n[i]+=ayij;
					acatz_n[i]+=azij;
					acatx_n[j]-=axij;
					acaty_n[j]-=ayij;
					acatz_n[j]-=azij;
					Ep_n+=k/rij;
				}
			}

			for (i=0;i<Nel;i++) {			//external electric field contributions
				ax_n[i]-=emFx;
				ay_n[i]-=emFy;
				az_n[i]-=emFz;
				Ep_n+=eFx*x_n[i]+eFy*y_n[i]+eFz*z_n[i];
			}
			for (j=0;j<Ncat;j++) {
				acatx_n[j]+=eMcatFx;
				acaty_n[j]+=eMcatFy;
				acatz_n[j]+=eMcatFz;
				Ep_n-=eFx*xcat_n[j]+eFy*ycat_n[j]+eFz*zcat_n[j];
			}

			for (i=0;i<Nel;i++) {			//magnetic field contributions (only for electrons)
				ax_n[i]-=emBz*vy_n[i]-emBy*vz_n[i];
				ay_n[i]-=emBx*vz_n[i]-emBz*vx_n[i];
				az_n[i]-=emBy*vx_n[i]-emBx*vy_n[i];
			}

			for (i=0;i<Nel;i++) {				//electron velocities at t+dt
				vx_n[i]=vx[i]+(ax[i]+ax_n[i])*dt_half;
				vy_n[i]=vy[i]+(ay[i]+ay_n[i])*dt_half;
				vz_n[i]=vz[i]+(az[i]+az_n[i])*dt_half;
			}
			for (j=0;j<Ncat;j++) {				//cation velocities at t+dt
				vcatx_n[j]=vcatx[j]+(acatx[j]+acatx_n[j])*dt_half;
				vcaty_n[j]=vcaty[j]+(acaty[j]+acaty_n[j])*dt_half;
				vcatz_n[j]=vcatz[j]+(acatz[j]+acatz_n[j])*dt_half;
			}

			Ekelt_n=0.;
			for (i=0;i<Nel;i++) {						//kinetic energies at t+dt
				Ek_n[i]=m_half*(vx_n[i]*vx_n[i]+vy_n[i]*vy_n[i]+vz_n[i]*vz_n[i]);
				Ekelt_n+=Ek_n[i];
			}
			Ekcatt_n=0.;
			for (j=0;j<Ncat;j++) {							//kinetic energies for cations
				Ekcat_n[j]=Mcat_half*(vcatx_n[j]*vcatx_n[j]+vcaty_n[j]*vcaty_n[j]+vcatz_n[j]*vcatz_n[j]);
				Ekcatt_n+=Ekcat_n[j];
			}
			Et_n=Ep_n+Ekelt_n+Ekcatt_n;
			Ekel_n_av=Ekelt_n/Nel;

			if (Ekel_av > kBT)
				Ekel_r=Ekel_av;
			else
				Ekel_r=kBT;
			r_E=fabs(Et_n-Et)/Ekel_r;	//check the energy drift
			if (r_E > g_E) {			//if it's too large, decrease dt and restart
				dt_high=1;
				dt_dec();
				if (dt_warn) 
					err_exit("traj(): Unable to reduce the time step");
				continue;
			}
			if (r_E > h_E) dt_increase=0; //if the drift is low, allow to increase dt

			if (r_tau_cat > h_tau) dt_increase=0; //if dt/tau_cat is not too large, allow to increase dt

			for (i=0;i<Nel;i++) {			//relative electron-atom velocities
				ux[i]=vx_n[i]-vMx[i];
				uy[i]=vy_n[i]-vMy[i];
				uz[i]=vz_n[i]-vMz[i];
				u[i]=sqrt(ux[i]*ux[i]+uy[i]*uy[i]+uz[i]*uz[i]);
				sigma_E[i]=crsec_E(u[i]);			//energy cross sections at t+dt
			}

			for (i=0;i<Nel;i++) {
				tau[i]=dens_inv/(u[i]*sigma_E[i]);		//mean free times to collisions
				r_tau[i]=dt/tau[i];			//check if dt is not too large compared to tau[i]
				if (r_tau[i] > g_tau) {		//if it is, decrease dt and restart
					dt_high=1;
					dt_dec();
					if (dt_warn) 
						err_exit("traj(): Unable to reduce the time step");
					break;
				}
				if (r_tau[i] > h_tau) dt_increase=0; //if it isn't, allow to increase dt
			}

		} while (dt_high);

		colls=0;
		for (i=0;i<Nel;i++) {
			P_ncoll=exp(-r_tau[i]);		//decide if an electron collision occurs
			R_coll=ran();
			if (R_coll > P_ncoll) {		//collision
				colls=1;
				coll[i]=1;
				sigma_p=crsec_p(u[i]);	//momentum cross section
				r_p=sigma_p/sigma_E[i];
				R_p=ran();
				if (R_p < r_p)			//decide the type of collission
					coll_E[i]=0;
				else
					coll_E[i]=1;
			} else
				coll[i]=0;
		}
		colls_cat=0;
		for (j=0;j<Ncat;j++) {
			R_coll=ran();
			if (R_coll > P_ncoll_cat) {		//decide if a cation collision occurs
				colls_cat=1;
				coll_cat[j]=1;
			} else
				coll_cat[j]=0;
		}

		t+=dt;				//increase time

		if (colls) {
			Ekelt=0.;
			for (i=0;i<Nel;i++) {
				if (coll[i]) {
					ctheta=ran()*2.-1.;
					stheta=sqrt(1.-ctheta*ctheta);
					phi=ran()*TWO_PI;
					nx=stheta*cos(phi);		//random vector of unit length
					ny=stheta*sin(phi);
					nz=ctheta;
					Msmu=Msm*u[i];
					vxc=Msmu*nx+msm*vx_n[i]+Msm*vMx[i];
					vyc=Msmu*ny+msm*vy_n[i]+Msm*vMy[i];
					vzc=Msmu*nz+msm*vz_n[i]+Msm*vMz[i];
					Ek[i]=m_half*(vxc*vxc+vyc*vyc+vzc*vzc);
					if (coll_E[i]) {
						r_vc_v=sqrt(Ek[i]/Ek_n[i]);		//energy collision
						vxc=r_vc_v*vx_n[i];
						vyc=r_vc_v*vy_n[i];
						vzc=r_vc_v*vz_n[i];
					}
					vx[i]=vxc;
					vy[i]=vyc;
					vz[i]=vzc;
				} else {					//no collision for this electron
					vx[i]=vx_n[i];
					vy[i]=vy_n[i];
					vz[i]=vz_n[i];
					Ek[i]=Ek_n[i];
				}
				Ekelt+=Ek[i];
			}
		} else {					//no electron collisions at all, update the velocities
			for (i=0;i<Nel;i++) {
				vx[i]=vx_n[i];
				vy[i]=vy_n[i];
				vz[i]=vz_n[i];
				Ek[i]=Ek_n[i];
			}
			Ekelt=Ekelt_n;
		}

		for (i=0;i<Nel;i++) {		//update the electron positions and accelerations
			x[i]=x_n[i];
			y[i]=y_n[i];
			z[i]=z_n[i];
			ax[i]=ax_n[i];
			ay[i]=ay_n[i];
			az[i]=az_n[i];
		}

		if (colls_cat) {
			Ekcatt=0.;
			for (j=0;j<Ncat;j++) {
				if (coll_cat[j]) {					//cation collisions
					vcatx[j]=kBTMcat*ran_gauss();
					vcaty[j]=kBTMcat*ran_gauss();
					vcatz[j]=kBTMcat*ran_gauss();
					Ekcat[j]=Mcat_half*(vcatx[j]*vcatx[j]+vcaty[j]*vcaty[j]+vcatz[j]*vcatz[j]);
				} else {							//no collision for this cation
					vcatx[j]=vcatx_n[j];
					vcaty[j]=vcaty_n[j];
					vcatz[j]=vcatz_n[j];
					Ekcat[j]=Ekcat_n[j];
				}
				Ekcatt+=Ekcat[j];
			}
		} else {							//no cation collisions at all, update the velocities
			for (j=0;j<Ncat;j++) {
				vcatx[j]=vcatx_n[j];
				vcaty[j]=vcaty_n[j];
				vcatz[j]=vcatz_n[j];
				Ekcat[j]=Ekcat_n[j];
			}
			Ekcatt=Ekcatt_n;
		}

		for (j=0;j<Ncat;j++) {		//update the cation positions and accelerations
			xcat[j]=xcat_n[j];
			ycat[j]=ycat_n[j];
			zcat[j]=zcat_n[j];
			acatx[j]=acatx_n[j];
			acaty[j]=acaty_n[j];
			acatz[j]=acatz_n[j];
		}
		for (kk=0;kk<kkmax;kk++)
			r_el_cat[kk]=r_el_cat_n[kk];
		Ep=Ep_n;
		Et=Ep+Ekelt+Ekcatt;
		Ekel_av=Ekelt/Nel;

		if (dt_increase) dt_inc();  //increase dt, if allowed
	}
	count_up();			//update the counters
	return;
}


// Set the initial time step and related variables
void dt_set()
{
	idt=DT_INIT;
	if (idt > ndt_max_cat)
		idt=ndt_max_cat;
	dt=dts[idt];
	dt_half=dts_half[idt];
	dt_sq_half=dts_sq_half[idt];
	r_tau_cat=r_tau_cats[idt];
	P_ncoll_cat=P_nccs[idt];
	dt_warn=0;
}


// Decrease the time step and related variables
void dt_dec()
{
	if (idt > 0) {
		idt--;
		dt=dts[idt];
		dt_half=dts_half[idt];
		dt_sq_half=dts_sq_half[idt];
		r_tau_cat=r_tau_cats[idt];
		P_ncoll_cat=P_nccs[idt];
	} else
		dt_warn++;
}


// Increase the time step and related variables
void dt_inc()
{
	if (idt < ndt_max_cat) {
		idt++;
		dt=dts[idt];
		dt_half=dts_half[idt];
		dt_sq_half=dts_sq_half[idt];
		r_tau_cat=r_tau_cats[idt];
		P_ncoll_cat=P_nccs[idt];
	}
}


// Check the reaction criterion
// Return values: GO_ON or FIN
int check_traj()
{
	int n_out,n_rec,n_gem,n_rm,n_canc;
	int i,j,kk,jm1;
	double r_sq;
	int el_out[NP],el_rm[NP],cat_rm[NP];
	int el_rec[N_EL_CAT],cat_rec[N_EL_CAT];
	char el_rm_t[NP];
	double xij,yij,zij,rij,rij_sq,axij,ayij,azij;
	double kmr3,Amcatr;
	double kMcatr3,AMccatr;

	n_out=0;
	n_rec=0;
	n_gem=0;
	for (i=0;i<Nel;i++) {
		r_sq=x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
		if (r_sq > Rmax_sq)						//check for escape
			el_out[n_out++]=i;
	}

	kk=0;
	for (i=0;i<Nel;i++) {
		if (Ek[i] < Ekcrit) {
			for (j=0;j<Ncat;j++) {
				if (r_el_cat[kk] < Rcrit) {	//check for reaction
					if (n_rec >= N_EL_CAT) 
						err_exit("check_traj(): Array size NP exceeded");
					el_rec[n_rec]=i;
					cat_rec[n_rec]=j;
					n_rec++;
				}
				kk++;
			}
		} else
			kk+=Ncat;
	}

	if (n_out==0 && n_rec==0) return GO_ON;		//no reactions or escapes - continue the simulation

	n_rm=0;
	for (i=0;i<Nel;i++)
		el_rm[i]=0;
	for (i=0;i<Ncat;i++)
		cat_rm[i]=0;

	if (n_out)
		for (kk=0;kk<n_out;kk++) {
			i=el_out[kk];
			el_rm[i]=1;		//mark the electrons to be removed due to escape
			el_rm_t[i]='E';
			n_rm++;
		}

	if (n_rec) {
		n_canc=0;
		for (kk=0;kk<n_rec;kk++) {
			i=cat_rec[kk];
			j=el_rec[kk];
			if (cat_rm[i] || el_rm[j])	//check for "double" reactions
				n_canc++;
			else {
				cat_rm[i]=1;		//mark the cations to be removed due to reaction
				el_rm[j]=1;			//mark the electrons to be removed due to reaction
				if (cat_id[i] == el_id[j]) {
					el_rm_t[j]='G';		//check if the reaction is geminate
					n_gem++;
				} else
					el_rm_t[j]='R';
				n_rm++;
			}
		}
		n_rec-=n_canc;
	}

	Nrec+=n_rec;
	Nout+=n_out;
	Ngem+=n_gem;

	for (i=0;i<Nel;i++)
		if (el_rm[i]) {

			printf("%c ",el_rm_t[i]);	//print info on the screen
			fprintf(fp_out,"%c ",el_rm_t[i]);	//print info in the output file

			t_rm[i_rm]=t;				//record data on reacting or escaping electrons
			type_rm[i_rm]=el_rm_t[i];
			Ek_rm[i_rm]=Ek[i];
			id_rm[i_rm]=el_id[i];
			x_rm[i_rm]=x[i];
			y_rm[i_rm]=y[i];
			z_rm[i_rm]=z[i];
			i_rm++;
		}

	if (n_rm==Nel) return FIN;			//end the simulation if no more electrons

	if (n_rec) {
		for (i=0;i<Ncat;i++)				//remove the reacting cations
			if (cat_rm[i]) {
				for (j=i+1;j<Ncat;j++) {
					jm1=j-1;
					xcat[jm1]=xcat[j];
					ycat[jm1]=ycat[j];
					zcat[jm1]=zcat[j];
					cat_id[jm1]=cat_id[j];
					vcatx[jm1]=vcatx[j];
					vcaty[jm1]=vcaty[j];
					vcatz[jm1]=vcatz[j];
					cat_rm[jm1]=cat_rm[j];
				}
				i--;
				Ncat--;
			}
		Ncatm1=Ncat-1;
	}

	for (i=0;i<Nel;i++)				//remove the reacting or escaping electrons
		if (el_rm[i]) {
			for (j=i+1;j<Nel;j++) {
				jm1=j-1;
				x[jm1]=x[j];
				y[jm1]=y[j];
				z[jm1]=z[j];
				el_id[jm1]=el_id[j];
				vx[jm1]=vx[j];
				vy[jm1]=vy[j];
				vz[jm1]=vz[j];
				el_rm[jm1]=el_rm[j];
			}
			i--;
			Nel--;
		}
	Nelm1=Nel-1;

// Calculate the accelerations and energies after the reactions or escapes
	Ep=0.;
	for (i=0;i<Nel;i++) ax[i]=ay[i]=az[i]=0.;
	for (j=0;j<Ncat;j++) acatx[j]=acaty[j]=acatz[j]=0.;

	kk=0;
	for (i=0;i<Nel;i++) {			//electron-cation contributions
		for (j=0;j<Ncat;j++) {
			xij=x[i]-xcat[j];
			yij=y[i]-ycat[j];
			zij=z[i]-zcat[j];
			rij_sq=xij*xij+yij*yij+zij*zij;
			r_el_cat[kk++]=rij=sqrt(rij_sq);
			if (rij > Rcat1) {
				kmr3=-km/(rij*rij_sq);
				ax[i]+=kmr3*xij;
				ay[i]+=kmr3*yij;
				az[i]+=kmr3*zij;
				kMcatr3=-kmr3*mMcat;
				acatx[j]+=kMcatr3*xij;
				acaty[j]+=kMcatr3*yij;
				acatz[j]+=kMcatr3*zij;
				Ep-=k/rij;
			} else if (rij > Rcat0) {
				Amcatr=Amcat*(Rcat0/rij-1.0);
				ax[i]+=Amcatr*xij;
				ay[i]+=Amcatr*yij;
				az[i]+=Amcatr*zij;
				AMccatr=-Amcatr*mMcat;
				acatx[j]+=AMccatr*xij;
				acaty[j]+=AMccatr*yij;
				acatz[j]+=AMccatr*zij;
				Ep+=Acat*rij_sq+Bcat*rij+Ccat;
			} else
				Ep+=EpRcat;
		}
	}
	kkmax=kk;

	for (i=0;i<Nelm1;i++) {			//electron-electron contributions
		for (j=i+1;j<Nel;j++) {
			xij=x[i]-x[j];
			yij=y[i]-y[j];
			zij=z[i]-z[j];
			rij_sq=xij*xij+yij*yij+zij*zij;
			rij=sqrt(rij_sq);
			kmr3=km/(rij*rij_sq);
			axij=kmr3*xij;
			ayij=kmr3*yij;
			azij=kmr3*zij;
			ax[i]+=axij;
			ay[i]+=ayij;
			az[i]+=azij;
			ax[j]-=axij;
			ay[j]-=ayij;
			az[j]-=azij;
			Ep+=k/rij;
		}
	}

	for (i=0;i<Ncatm1;i++) {			//cation-cation contributions
		for (j=i+1;j<Ncat;j++) {
			xij=xcat[i]-xcat[j];
			yij=ycat[i]-ycat[j];
			zij=zcat[i]-zcat[j];
			rij_sq=xij*xij+yij*yij+zij*zij;
			rij=sqrt(rij_sq);
			kMcatr3=kMcat/(rij*rij_sq);
			axij=kMcatr3*xij;
			ayij=kMcatr3*yij;
			azij=kMcatr3*zij;
			acatx[i]+=axij;
			acaty[i]+=ayij;
			acatz[i]+=azij;
			acatx[j]-=axij;
			acaty[j]-=ayij;
			acatz[j]-=azij;
			Ep+=k/rij;
		}
	}

	for (i=0;i<Nel;i++) {			//external electric field contributions
		ax[i]-=emFx;
		ay[i]-=emFy;
		az[i]-=emFz;
		Ep+=eFx*x[i]+eFy*y[i]+eFz*z[i];
	}
	for (j=0;j<Ncat;j++) {
		acatx[j]+=eMcatFx;
		acaty[j]+=eMcatFy;
		acatz[j]+=eMcatFz;
		Ep-=eFx*xcat[j]+eFy*ycat[j]+eFz*zcat[j];
	}

	for (i=0;i<Nel;i++) {			//magnetic field contributions (only for electrons)
		ax[i]-=emBz*vy[i]-emBy*vz[i];
		ay[i]-=emBx*vz[i]-emBz*vx[i];
		az[i]-=emBy*vx[i]-emBx*vy[i];
	}

	Ekelt=0.;
	for (i=0;i<Nel;i++) {							//kinetic energies for electrons
		Ek[i]=m_half*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
		Ekelt+=Ek[i];
	}
	Ekcatt=0.;
	for (j=0;j<Ncat;j++) {							//kinetic energies for cations
		Ekcat[j]=Mcat_half*(vcatx[j]*vcatx[j]+vcaty[j]*vcaty[j]+vcatz[j]*vcatz[j]);
		Ekcatt+=Ekcat[j];
	}
	Et=Ep+Ekelt+Ekcatt;
	Ekel_av=Ekelt/Nel;

	return GO_ON;
}


// Update the counters
void count_up()
{
	int i;

	if (Nrec+Nout!=i_rm) err_exit("count_up(): Error in recorded events");

	Nrec_tot+=Nrec;
	Nout_tot+=Nout;
	Ngem_tot+=Ngem;
	for (i=0;i<i_rm;i++) {
		if (i_re<N_RE) {
			run_re[i_re]=irun;
			t_re[i_re]=t_rm[i];
			type_re[i_re]=type_rm[i];
			Ek_re[i_re]=Ek_rm[i];
			id_re[i_re]=id_rm[i];
			x_re[i_re]=x_rm[i];
			y_re[i_re]=y_rm[i];
			z_re[i_re]=z_rm[i];
			i_re++;
		} else
			printf("count_up(): N_RE exceeded, no detailed output");
	}
}


/* Simulation input */
#define	MAX_LINE	120
void input()
{
	FILE *fp_in;
	char comment_line[MAX_LINE];
	int i;
	int ndtm_nset;
	double Fd2,Ffac,Bd2,Bfac;

	fp_in=fopen(IN_FNAME,"r");				//open input file
	fscanf(fp_in,"%d",&Npair);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&rk_A);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%d",&rk_dist);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&cat_mob_cm2Vs);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&cat_mass_MAr);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&F_Vcm);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Fd_x);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Fd_y);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Fd_z);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&B_T);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Bd_x);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Bd_y);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Bd_z);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%d",&Ek0_dist);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&r0_A);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Ek0_eV);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Rcrit_A);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Ekcrit_eV);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&Rmax_A);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%ld",&Nrun);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%ld",&Nav);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&g_E);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&h_E);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&g_tau);
	fgets(comment_line,MAX_LINE,fp_in);
	fscanf(fp_in,"%lf",&h_tau);
	fclose(fp_in);

	if (Npair > NP)
		err_exit("input(): Npair too large, increase NP");
	if (rk_A < R_MIN_CAT)
		err_exit("input(): rk is too short");
	rk=rk_A*1.e-10;
	cat_mob=cat_mob_cm2Vs*1.e-4;
	r0=r0_A*1.e-10;
	if (Ek0_dist == 0) {
		if (Ek0_eV > EK0_EV_MAX)
			err_exit("input(): Ek0 is too high");
		Ek0=Ek0_eV*EL_CHG;
	}
	Rcrit=Rcrit_A*1.e-10;
	Ekcrit=Ekcrit_eV*EL_CHG;
	Rmax=Rmax_A*1.e-10;

	F=F_Vcm*100.;
	if (F < 1.e-20) {
		Fx=0.;
		Fy=0.;
		Fz=0.;
	} else {
		Fd2=Fd_x*Fd_x+Fd_y*Fd_y+Fd_z*Fd_z;
		if (Fd2 < 1.e-6)
			err_exit("input(): wrong F direction");
		Ffac=F/sqrt(Fd2);
		Fx=Ffac*Fd_x;
		Fy=Ffac*Fd_y;
		Fz=Ffac*Fd_z;
	}

	B=B_T;
	if (B < 1.e-20) {
		Bx=0.;
		By=0.;
		Bz=0.;
	} else {
		Bd2=Bd_x*Bd_x+Bd_y*Bd_y+Bd_z*Bd_z;
		if (Bd2 < 1.e-6)
			err_exit("input(): wrong B direction");
		Bfac=B/sqrt(Bd2);
		Bx=Bfac*Bd_x;
		By=Bfac*Bd_y;
		Bz=Bfac*Bd_z;
	}

	m_half=0.5*EL_MASS;
	M=M_AR/N_AV*1.e-3;
	Mcat=cat_mass_MAr*M;
	Mcat_half=0.5*Mcat;
	mMcat=EL_MASS/Mcat;
	Msm=M/(EL_MASS+M);
	msm=EL_MASS/(EL_MASS+M);
	k=EL_CHG*EL_CHG/(4.0*PI*EPS_0*EPS);
	two_k=2.0*k;
	km=k/EL_MASS;
	kMcat=k/Mcat;
	Rcat=R_CAT*1.e-10;
	EpRcat=-k/Rcat;
	Rcat0=R_CAT_0*1.e-10;
	Rcat1=Rcat*(0.75+0.5*sqrt(2.25-2.0*Rcat0/Rcat));
	Acat=k/(2.0*Rcat1*Rcat1*(Rcat1-Rcat0));
	Bcat=-2.0*Acat*Rcat0;
	Ccat=Acat*Rcat0*Rcat0+EpRcat;
	Amcat=2.0*Acat/EL_MASS;
	E_Coul_r0=-k/r0;
	kBT=KB*TEMP;
	kBTM=sqrt(kBT/M);
	kBTMcat=sqrt(kBT/Mcat);
	dens_inv=1.0/DENSITY;
	tau_cat=cat_mob*Mcat/EL_CHG;
	eFx=EL_CHG*Fx;
	emFx=eFx/EL_MASS;
	eMcatFx=eFx/Mcat;
	eFy=EL_CHG*Fy;
	emFy=eFy/EL_MASS;
	eMcatFy=eFy/Mcat;
	eFz=EL_CHG*Fz;
	emFz=eFz/EL_MASS;
	eMcatFz=eFz/Mcat;
	emBx=EL_CHG*Bx/EL_MASS;
	emBy=EL_CHG*By/EL_MASS;
	emBz=EL_CHG*Bz/EL_MASS;
	Rmax_sq=Rmax*Rmax;
	Rmin_cat=R_MIN_CAT*1.e-10;

	if (g_tau*tau_cat < dts_fs[0]*1.e-15)
		err_exit("input(): Cation mobility is too low");
	ndt_max_cat=NDT-1;
	ndtm_nset=1;
	for (i=0;i<NDT;i++) {
		dts[i]=dts_fs[i]*1.e-15;
		dts_half[i]=0.5*dts[i];
		dts_sq_half[i]=0.5*dts[i]*dts[i];
		r_tau_cats[i]=dts[i]/tau_cat;
		if (ndtm_nset)
			if (r_tau_cats[i] > g_tau) {
				ndt_max_cat=i-1;
				ndtm_nset=0;
			}
		P_nccs[i]=exp(-r_tau_cats[i]);
	}
}


/* Simulation output */
void output(int mode)
{
	FILE *fp_res,*fp_det;
	int i;

	switch (mode) {
	case 1:
		printf("--------------------------------------------------------------------------------\n");
		printf("   Simulation of electron-ion recombination in liquid argon at 87 K\n");
		printf("   Npair=%2d  rk,A=%6.1f  rk_dist=%1d  Rmax,A=%8.1f  Nrun=%7d  Nav=%6d\n",Npair,rk_A,rk_dist,Rmax_A,Nrun,Nav);
		printf("   cat_mob,cm2/Vs=%10.4e  cat_mass/MAr=%6.2f\n",cat_mob_cm2Vs,cat_mass_MAr);
		printf("   External electric field: F,V/cm=%6.1f  Fd_x=%7.4f  Fd_y=%7.4f  Fd_z=%7.4f\n",F_Vcm,Fd_x,Fd_y,Fd_z);
		printf("   Magnetic field: B,T=%8.4f  Bd_x=%7.4f  Bd_y=%7.4f  Bd_z=%7.4f\n",B_T,Bd_x,Bd_y,Bd_z);
		printf("--------------------------------------------------------------------------------\n");
		char toAddFile[300];
		sprintf(toAddFile, "LAr_dm_B(%ld).out", currentSeed);
		fp_out=fopen(toAddFile,"w");
		fprintf(fp_out,"--------------------------------------------------------------------------------\n");
		fprintf(fp_out,"   Simulation of electron-ion recombination in liquid argon at 87 K\n");
		fprintf(fp_out,"   Npair=%2d  rk,A=%6.1f  rk_dist=%1d  Rmax,A=%8.1f  Nrun=%7d  Nav=%6d\n",Npair,rk_A,rk_dist,Rmax_A,Nrun,Nav);
		fprintf(fp_out,"   cat_mob,cm2/Vs=%10.4e  cat_mass/MAr=%6.2f\n",cat_mob_cm2Vs,cat_mass_MAr);
		fprintf(fp_out,"   External electric field: F,V/cm=%6.1f  Fd_x=%7.4f  Fd_y=%7.4f  Fd_z=%7.4f\n",F_Vcm,Fd_x,Fd_y,Fd_z);
		fprintf(fp_out,"   Magnetic field: B,T=%8.4f  Bd_x=%7.4f  Bd_y=%7.4f  Bd_z=%7.4f\n",B_T,Bd_x,Bd_y,Bd_z);
		fprintf(fp_out,"--------------------------------------------------------------------------------\n");
		n_av=0;
		sum_Prec=0.0;
		sum2_Prec=0.0;
		sum_Pgem=0.0;
		sum2_Pgem=0.0;
		Nrec_tot=0;
		Nout_tot=0;
		Ngem_tot=0;
		i_re=0;
		sprintf(toAddFile, "LAr_dm_B(%ld).det", currentSeed);
		fp_det=fopen(toAddFile,"w");
		fprintf(fp_det,"LAr_dm_B: Npair=%2d  rk,A=%6.1f  rk_dist=%1d  Rmax,A=%8.1f  Nrun=%7d  Nav=%6d\n",Npair,rk_A,rk_dist,Rmax_A,Nrun,Nav);
		fprintf(fp_det,"            cat_mob,cm2/Vs=%10.4e  cat_mass/MAr=%6.2f\n",cat_mob_cm2Vs,cat_mass_MAr);
		fprintf(fp_det,"            F,V/cm=%6.1f  Fd_x=%10.6f  Fd_y=%10.6f  Fd_z=%10.6f\n",F_Vcm,Fd_x,Fd_y,Fd_z);
		fprintf(fp_det,"            B,T=%8.4f  Bd_x=%10.6f  Bd_y=%10.6f  Bd_z=%10.6f\n",B_T,Bd_x,Bd_y,Bd_z);
		fprintf(fp_det,"   run \t          t,s \t \t   Ek,eV \tel_id \t    x,m \t    y,m \t    z,m \n");
		fclose(fp_det);
		break;
	case 2:
		Prec=(double)Nrec_tot/(Nav*Npair);
		if (n_av<NTAB) Prec_tab[n_av]=Prec;
		Pgem=(double)Ngem_tot/(Nav*Npair);
		if (n_av<NTAB) Pgem_tab[n_av]=Pgem;
		n_av++;
		sum_Prec+=Prec;
		sum2_Prec+=Prec*Prec;
		Prec_av=sum_Prec/n_av;
		if (n_av>1) Prec_stderr=sqrt((sum2_Prec/n_av-Prec_av*Prec_av)/(n_av-1));
		sum_Pgem+=Pgem;
		sum2_Pgem+=Pgem*Pgem;
		Pgem_av=sum_Pgem/n_av;
		if (n_av>1) Pgem_stderr=sqrt((sum2_Pgem/n_av-Pgem_av*Pgem_av)/(n_av-1));
		printf("-----------------------------------------------\n");
		printf(" RUN No. %7d\n",irun);
		printf(" Nrec=  %d       Ngem=  %d\n",Nrec_tot,Ngem_tot);
		printf(" Prec:  av / stderr\n");
		if (n_av>1)
			printf("     %9.6f / %9.6f\n",Prec_av,Prec_stderr);
		else
			printf("     %9.6f\n",Prec_av);
		printf(" Pgem:  av / stderr\n");
		if (n_av>1)
			printf("     %9.6f / %9.6f\n",Pgem_av,Pgem_stderr);
		else
			printf("     %9.6f\n",Pgem_av);
		printf("-----------------------------------------------\n");
		fprintf(fp_out,"-----------------------------------------------\n");
		fprintf(fp_out," RUN No. %7d\n",irun);
		fprintf(fp_out," Nrec=  %d       Ngem=  %d\n",Nrec_tot,Ngem_tot);
		fprintf(fp_out," Prec:  av / stderr\n");
		if (n_av>1)
			fprintf(fp_out,"     %9.6f / %9.6f\n",Prec_av,Prec_stderr);
		else
			fprintf(fp_out,"     %9.6f\n",Prec_av);
		fprintf(fp_out," Pgem:  av / stderr\n");
		if (n_av>1)
			fprintf(fp_out,"     %9.6f / %9.6f\n",Pgem_av,Pgem_stderr);
		else
			fprintf(fp_out,"     %9.6f\n",Pgem_av);
		fprintf(fp_out,"-----------------------------------------------\n");

		sprintf(toAddFile, "LAr_dm_B(%ld).res", currentSeed);
		fp_res=fopen(toAddFile,"w");
		fprintf(fp_res,"LAr_dm_B: Npair=%2d  rk,A=%6.1f  rk_dist=%1d  Rmax,A=%8.1f  Nrun=%7d  Nav=%6d\n",Npair,rk_A,rk_dist,Rmax_A,Nrun,Nav);
		fprintf(fp_res,"            cat_mob,cm2/Vs=%10.4e  cat_mass/MAr=%6.2f\n",cat_mob_cm2Vs,cat_mass_MAr);
		fprintf(fp_res,"            F,V/cm=%6.1f  Fd_x=%10.6f  Fd_y=%10.6f  Fd_z=%10.6f\n",F_Vcm,Fd_x,Fd_y,Fd_z);
		fprintf(fp_res,"            B,T=%8.4f  Bd_x=%10.6f  Bd_y=%10.6f  Bd_z=%10.6f\n",B_T,Bd_x,Bd_y,Bd_z);
		fprintf(fp_res,"Output at run no. %d\n",irun);
		fprintf(fp_res,"Recorded Prec values:\n");
		for (i=0;i<n_av;i++)
			fprintf(fp_res,"%9.6f\n",Prec_tab[i]);
		fprintf(fp_res,"Prec:  av / stderr\n");
		if (n_av>1)
			fprintf(fp_res,"%9.6f / %9.6f\n",Prec_av,Prec_stderr);
		else
			fprintf(fp_res,"%9.6f\n",Prec_av);
		fprintf(fp_res,"\nRecorded Pgem values:\n");
		for (i=0;i<n_av;i++)
			fprintf(fp_res,"%9.6f\n",Pgem_tab[i]);
		fprintf(fp_res,"Pgem:  av / stderr\n");
		if (n_av>1)
			fprintf(fp_res,"%9.6f / %9.6f\n",Pgem_av,Pgem_stderr);
		else
			fprintf(fp_res,"%9.6f\n",Pgem_av);
		fclose(fp_res);

		fp_det=fopen(DET_FNAME,"a");
		for (i=0;i<i_re;i++)
			fprintf(fp_det,"%7d\t%14.6e\t%c\t%9.6f\t%3d\t%14.6e\t%14.6e\t%14.6e\n",run_re[i],t_re[i],type_re[i],Ek_re[i]/EL_CHG,id_re[i],x_re[i],y_re[i],z_re[i]);
		fclose(fp_det);
				
		Nrec_tot=0;
		Nout_tot=0;
		Ngem_tot=0;
		i_re=0;
		break;
	case 3:
		printf("Simulation finished\n");
		fprintf(fp_out,"Simulation finished\n");
		fclose(fp_out);
		break;
	}
	
}


void err_exit(char err_txt[])
{
	fprintf(stderr,"%s\n",err_txt);
	fprintf(fp_out,"%s\n",err_txt);
	fclose(fp_out);
	exit(1);
}