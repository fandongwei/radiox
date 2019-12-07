/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __RadioXRunner__h__
#define __RadioXRunner__h__
#include "SpatialXResultTable.h"
#include "RadioXResultTable.h"
#include <string>
class LobePrior;
class RadioXRunner{
public:
	static void doRadioX(bool sym, LobePrior* pri, double* sigmas, SpatialXResultTable& tbl, std::string outfile, int sampleCount);
	static void doRadioX(bool sym, LobePrior* pri, double* sigmas, std::string infile, std::string outfile, int sampleCount);
	static void analysisResult(RadioXResultTable& tbl, std::string outfile);
	//static void analysisResult0(RadioXResultTable& tbl, std::string outfile);
	static void analysisResult(std::string infile, std::string outfile);
};
#endif