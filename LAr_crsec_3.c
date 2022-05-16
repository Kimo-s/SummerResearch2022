#include "LAr_crsec.h"
#include <stdio.h>

#define	N_VSV	3265
#define	N_VSV1	3266
#define	DV		0.05
#define	UNIT_V	3.6313e4
#define	MDV_MS	1.0/(DV*UNIT_V)

double waE[N_VSV1],wbE[N_VSV1];
double wap[N_VSV1],wbp[N_VSV1];
double watot[N_VSV1],wbtot[N_VSV1];


/* Initializes the cross section data */
void crsec_init()
{
	FILE *sigma_fp;
	double sigmatot1,sigmatot2;
	double v2,sigmaE2,sigmap2;
	int i;

	sigma_fp=fopen("sigma_LArE_3.ab","r");
	for (i=0;i<N_VSV1;i++)
		fscanf(sigma_fp,"%lf%lf",&waE[i],&wbE[i]);
	fclose(sigma_fp);

	sigma_fp=fopen("sigma_LArp_3.ab","r");
	for (i=0;i<N_VSV1;i++)
		fscanf(sigma_fp,"%lf%lf",&wap[i],&wbp[i]);
	fclose(sigma_fp);

	printf("Cross section data read O.K.\n");

// calculate "total" cross section [=max(sigma_E,sigma_p)]
	sigmatot1=wbE[0];
	for(i=0;i<N_VSV1;i++) {
		v2=(i+1)*DV;
		sigmaE2=waE[i]*v2+wbE[i];
		sigmap2=wap[i]*v2+wbp[i];
		sigmatot2 = sigmaE2>sigmap2 ? sigmaE2 : sigmap2;
		watot[i]=(sigmatot2-sigmatot1)/DV;
		wbtot[i]=sigmatot2-watot[i]*v2;
		sigmatot1=sigmatot2;
	}

	for(i=0;i<N_VSV1;i++) {
		watot[i]/=UNIT_V;
		wap[i]/=UNIT_V;
	}
	
}

/* Returns the cross section sigma_tot [m^2] for a given velocity v_ms [m/s] */
double crsec_E(double v_ms)
{
	int ind;

	ind=(int)(v_ms*MDV_MS);
	if (ind>N_VSV) ind=N_VSV;
	return watot[ind]*v_ms+wbtot[ind];
}

/* Returns the cross section sigma_p [m^2] for a given velocity v_ms [m/s] */
double crsec_p(double v_ms)
{
	int ind;

	ind=(int)(v_ms*MDV_MS);
	if (ind>N_VSV) ind=N_VSV;
	return wap[ind]*v_ms+wbp[ind];
}
