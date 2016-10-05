/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __GaussianRandom2D__h__
#define __GaussianRandom2D__h__
#include <Eigen/Dense>
class Ran1;
class GaussianRandom2D{
public:
	GaussianRandom2D(Eigen::MatrixXd* mean, Eigen::MatrixXd* sigma, Ran1* ran);
public:
	void nextSample(double&x ,double& y);
	Eigen::VectorXd nextSample();
private:
	Eigen::MatrixXd* M;
	Eigen::MatrixXd* Sigma;
	Ran1* RAN;
	long Seed();
};
#endif