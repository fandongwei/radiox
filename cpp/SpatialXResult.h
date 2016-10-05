/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __SpatialXResult__h__
#define __SpatialXResult__h__

#include "SpatialPoint.h"
#include <string>

//one pair result of xmatch
class SpatialXResult
{
public:
	//the two points
	SpatialPoint p1, p2;
	//distance of the two points
	double dist;
public:
	std::string toString();
	bool operator <(const SpatialXResult& p) const;
};
#endif