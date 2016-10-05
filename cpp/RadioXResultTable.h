/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __RadioXResultTable__h__
#define __RadioXResultTable__h__

#include <vector>
#include <string>
#include <map>
#include "RadioXResult.h"

//point list
class RadioXResultTable {
public:
	std::vector<RadioXResult> data;
public:
	std::map<long, double> corelikelihood;
	void init(bool sortbylikelihood);
public:
	RadioXResultTable(std::string filename, bool sortbylikelihood=false);
	RadioXResultTable(std::vector<RadioXResult>& d, bool sortbylikelihood=false);
	~RadioXResultTable();
public:
	bool nextTopRadioXResult(RadioXResult& rxr);
	bool nextTopRadioXResult(RadioXResult& rxr, std::map<long, int>& opticalremovedcomponents, std::map<long, int>& radioremovedcomponents);
	void removeTopRadioxResult();
	//void removeRelatedRadioXResult(RadioXResult& rxr);
	std::vector<RadioXResult> getResultsByType(int type);
	std::vector<RadioXResult> getTriplets();
	std::vector<RadioXResult> getDoublets();
	std::vector<RadioXResult> getCores();
};
#endif