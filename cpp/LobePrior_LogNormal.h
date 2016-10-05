/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __LobePrior_LogNormal__h__
#define __LobePrior_LogNormal__h__
#include "LobePrior.h"
class LobePrior_LogNormal:public LobePrior {
public:
	double mu;
	double sigma;
	double sigma2;
	double sigma2_m2;
	double resigmapi2sqrt;
	double mean;
	double s;
public:
	LobePrior_LogNormal(double mu, double sigma);
	static LobePrior_LogNormal getLogNormal(double mean, double variance);
public:
	virtual double prior(double r2);
	virtual double variance();
};
#endif