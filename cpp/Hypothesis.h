/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __Hypothesis__h__
#define __Hypothesis__h__
#include <vector>
#include <string>
#include <Eigen/Dense>
class LobePrior;
class SpatialPoint;
class Hypothesis {
public:
	Hypothesis();
	~Hypothesis();
public:
	//optical point
	SpatialPoint* optical;
	//radio core point
	SpatialPoint* rcore;
	//radio lobe point 1
	SpatialPoint* rlobes1;
	//radio lobe point 2
	SpatialPoint* rlobes2;
	//Eigen::VectorXd result;
	double likelihood, likelihood_err;
public:
	double get1stPrior();
	void likelihoodSym(LobePrior* pri, int* sampleCounts, double& likelihood, double& likelihood_err);
	void likelihoodAsym(LobePrior* pri, int* sampleCounts, double ksigma, double& likelihood, double& likelihood_err);
	//generate all possible combinations
	static std::vector<Hypothesis> generateHypList(SpatialPoint* optical, std::vector<SpatialPoint*>& radios);
	void print();
	std::string combiStr();
	int getType();//1:core,2:lobe1,4:lobe2. mask 7:core,lobe,lobe; 6:lobe,lobe; 3:core,lobe
	std::string toString();
	bool operator <(const Hypothesis& p);
};
#endif