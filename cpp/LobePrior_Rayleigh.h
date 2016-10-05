/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __LobePrior_Rayleigh__h__
#define __LobePrior_Rayleigh__h__
#include "LobePrior.h"
class LobePrior_Rayleigh:public LobePrior {
public:
	double sigma;
	double resigma2;
	double s;
	double mean;
public:
	LobePrior_Rayleigh(double sigma);
public:
	virtual double prior(double r2);
	virtual double variance();
	double getSigmaFromMean(double mean);
};
#endif