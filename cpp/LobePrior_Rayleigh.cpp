/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "LobePrior_Rayleigh.h"
#include "Constants.h"
#include <math.h>

LobePrior_Rayleigh::LobePrior_Rayleigh(double sigma){
	this->sigma = sigma;
	this->mean = sigma * sqrt(Constants::PI / 2.0);
	this->resigma2 = pow(sigma, -2);
	this->s = 1 / 2.0 / Constants::PI / mean;
}
double LobePrior_Rayleigh::prior(double r2){
	double x = sqrt(r2);
	return x * resigma2 * exp(-x * x / 2 * resigma2)*s;
}
double LobePrior_Rayleigh::variance(){
	return (4.0 - Constants::PI) / 2.0 / resigma2;
}
double LobePrior_Rayleigh::getSigmaFromMean(double mean){
	return mean * sqrt(2.0 / Constants::PI);
}