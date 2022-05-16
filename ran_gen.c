#include "ran_gen.h"
#include <math.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)


static long idum1,idum2;
static long iy;
static long iv[NTAB];


void ran_init(long iseed)
{
	int j;
	long k;

	idum2=idum1=iseed;
	for (j=NTAB+7;j>=0;j--) {
		k=idum1/IQ1;
		idum1=IA1*(idum1-k*IQ1)-k*IR1;
		if (idum1 < 0) idum1 += IM1;
		if (j < NTAB) iv[j] = idum1;
	}
	iy=iv[0];
	return;
}


double ran()
{
	int j;
	long k;

	k=idum1/IQ1;
	idum1=IA1*(idum1-k*IQ1)-k*IR1;
	if (idum1 < 0) idum1 += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum1;
	if (iy < 1) iy += IMM1;
	return AM*iy;
}


double ran_exp()
{
	return -log(ran());
}


double ran_gauss()
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (iset == 0) {
		do {
			v1=2.0*ran()-1.0;
			v2=2.0*ran()-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v2*fac;
		iset=1;
		return v1*fac;
	} else {
		iset=0;
		return gset;
	}
}
