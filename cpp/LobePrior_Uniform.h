/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __LobePrior_Uniform__h__
#define __LobePrior_Uniform__h__
#include "LobePrior.h"
class LobePrior_Uniform:public LobePrior {
public:
	double min2;
	double max2;
	double min;
	double max;
	double mean;
	double uni;
	double s;//the integral is not 1, so have to divide the S.
public:
	LobePrior_Uniform(double min, double max);
public:
	virtual double prior(double x2);
	virtual double variance();
};
#endif