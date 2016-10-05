/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "Ran1.h"
#include <math.h>
#include <stdio.h>
Ran1::Ran1(){
	init(1);
}
Ran1::Ran1(long seed){
	init(seed);
}
Ran1::~Ran1(){
}
void Ran1::init(long seed){
	// ran1
	IA = 16807;
	IM = 2147483647;
	AM = 1.0/2147483647.0;
	IQ = 127773;
	IR = 2836;
	NTAB = 32;
	NDIV = 1.0+(2147483647.0-1.0)/32.0;
	RNMX = 1.0-1.2e-7;
	iy=0;  
	idum = seed;
	// gasdev
	iset = 0;
}
double Ran1::nextDouble(){
	int j;
	long k;
	double temp;
	if (idum <= 0 || iy == 0){
		 if (-idum < 1) idum=1;  
	     else idum = -idum;
	     for (j=NTAB+7;j>=0;j--){   
			 k=idum/IQ;
	         idum=IA*(idum-k*IQ)-IR*k;
	         if (idum < 0) idum += IM;
	         if (j < NTAB) iv[j] = idum;
	     }
	     iy=iv[0];
	}
	k=idum/IQ;  
	idum=IA*(idum-k*IQ)-IR*k; 
	if (idum < 0) idum += IM;
	//printf("%ld, %d, %.8lf, %.8lf, %d, %ld, %ld\n", iy, NTAB, RNMX, NDIV, j, idum, k);
	j=((int)(iy/NDIV+0.5))%NTAB;//C++ not the same as C#
	//j=((int)(iy/NDIV))%NTAB;
	iy=iv[j];                  
	iv[j] = idum;
	//printf("%ld, %d, %.8lf, %.8lf, %d, %ld, %ld\n", iy, NTAB, RNMX, NDIV, j, idum, k);
	if ((temp=AM*iy) > RNMX) return RNMX;  
	else return temp;
}
double Ran1::nextDouble(double min, double max){
	return min + (max - min) * nextDouble();
}
double Ran1::nextGasdev(){
	double fac, rsq, v1, v2;
	if (iset == 0){
	    do{
	        v1 = 2.0 * nextDouble() - 1.0;
	        v2 = 2.0 * nextDouble() - 1.0;
	        rsq = v1 * v1 + v2 * v2;
	    } while (rsq >= 1.0 || rsq == 0.0);
	    fac = sqrt(-2.0 * log(rsq) / rsq);
	    gset = v1 * fac;
	    iset = 1;
	    return v2 * fac;
	} else {
	    iset = 0;
	    return gset;
	}
}
double Ran1::nextGasdev(double mean, double sigma){
	return nextGasdev() * sigma + mean;
}