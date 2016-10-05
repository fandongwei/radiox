/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "LobePrior_Uniform.h"
#include "Constants.h"

LobePrior_Uniform::LobePrior_Uniform(double min, double max){
	this->min = min;
	this->max = max;
	this->min2 = min*min;
	this->max2 = max*max;
	this->mean = (min+max)/2.0;
	this->s = 1.0/2.0/Constants::PI/this->mean;
	this->uni = 1.0/(max-min);
	this->uni *= s;
}
double LobePrior_Uniform::prior(double r2){
	if (r2 > min2 && r2 < max2) return this->uni;
	else return 0;
}
double LobePrior_Uniform::variance(){
	double sub = max - min;
	return sub * sub / 12.0;
}
