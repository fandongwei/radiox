/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __RadioXResult__h__
#define __RadioXResult__h__

#include "SpatialPoint.h"
#include <string>
#include <vector>

//one pair result of xmatch
class RadioXResult
{
public:
	//the two points
	long opticalid;
	int combitype;
	std::string combi;
	std::string radioids;
	double likelihood;
	double likelihood_err;
public:
	std::vector<long> radiocomponents;//only the radio components
public:
	std::string toString();
	bool operator <(const RadioXResult& p) const;
};
#endif