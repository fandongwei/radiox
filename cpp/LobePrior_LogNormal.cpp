/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "LobePrior_LogNormal.h"
#include "Constants.h"
#include <math.h>

LobePrior_LogNormal::LobePrior_LogNormal(double mu, double sigma){
	this->mu = mu;
	this->sigma = sigma;
	this->sigma2_m2 = sigma * sigma * 2.0;
	this->resigmapi2sqrt = 1.0 / sqrt(sigma2_m2 * Constants::PI);
	this->mean = exp(mu + sigma * sigma / 2.0);
	this->s = 1 / 2.0 / Constants::PI / mean;
}
LobePrior_LogNormal LobePrior_LogNormal::getLogNormal(double mean, double variance){
	double sigma = sqrt(log(1 + variance / mean / mean));
	double mu = log(mean) - sigma * sigma / 2.0;
	return LobePrior_LogNormal(mu, sigma);
}
double LobePrior_LogNormal::prior(double r2){
	if(r2==0) return 0;
	
	double x = sqrt(r2);
	
	double e = log(x)-mu;
	e = e * e / sigma2_m2;
	e = exp(-e);
	return resigmapi2sqrt * e / x * s;
}
double LobePrior_LogNormal::variance(){
	return (exp(sigma * sigma) - 1.0) * exp(2 * mu + sigma * sigma);
}