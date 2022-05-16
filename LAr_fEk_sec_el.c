#include "LAr_fEk_sec_el.h"
#include <stdio.h>

#define	N_EK	1000
#define	N_EK1	999
#define	DEK		0.1
#define	MDEK	10.

double wase[N_EK],wbse[N_EK];

/* Initializes the Ek distribution data */
void fEk_sec_el_init()
{
	FILE *fEk_fp;
	int i;

	fEk_fp=fopen("fEk_sec_el_LAr.ab","r");
	for (i=0;i<N_EK;i++)
		fscanf(fEk_fp,"%lf%lf",&wase[i],&wbse[i]);
	fclose(fEk_fp);
	printf("Secondary electron energy distribution data read O.K.\n");
}


/* Returns the value of the secondary electron energy distribution for a given energy [eV] */
double fEk(double Ek_eV)
{
	int ind;

	ind=(int)(Ek_eV*MDEK);
	if (ind > N_EK1) ind=N_EK1;
	return wase[ind]*Ek_eV+wbse[ind];
}
