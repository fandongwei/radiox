/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __SpatialXResultTable__h__
#define __SpatialXResultTable__h__

#include <vector>
#include <string>
#include "SpatialXResult.h"

//point list
class SpatialXResultTable {
public:
	std::vector<SpatialXResult> data;
public:
	SpatialXResultTable(std::string filename, bool sortbydec=false);
	SpatialXResultTable(std::vector<SpatialXResult>& d, bool sortbydec=false);
	~SpatialXResultTable();
private:
	std::vector<SpatialXResult>::iterator currentIndex;
public:
	std::vector<SpatialXResult> getNextCombination();
};
#endif