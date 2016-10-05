/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "RadioXResult.h"
#include<stdio.h>
#include <iostream>
#include <sstream>

std::string RadioXResult::toString() {
	char buf[500];
	sprintf(buf, "%ld,%d,%s,%s,%.17lf,%.17lf", opticalid, combitype, combi.c_str(), radioids.c_str(), likelihood, likelihood_err);
	return std::string(buf);
	//std::stringstream ss;
	//ss << opticalid << ',' << combitype << ',' << combi << ',' << radioids << ',' << likelihood << ',' << likelihood_err;
	//return ss.str();
}
bool RadioXResult::operator <(const RadioXResult& x) const{
	//return this->likelihood > x.likelihood;//large first
	return this->likelihood < x.likelihood;//small first
}