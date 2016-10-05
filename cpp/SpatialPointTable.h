/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __SpatialPointTable__h__
#define __SpatialPointTable__h__

#include <vector>
#include <string>
#include "SpatialPoint.h"

//point list
class SpatialPointTable {
public:
	std::vector<SpatialPoint> data;
public:
	SpatialPointTable(std::string filename, bool sortbydec=false);
	SpatialPointTable(std::vector<SpatialPoint>& d, bool sortbydec=false);
	~SpatialPointTable();
};
#endif