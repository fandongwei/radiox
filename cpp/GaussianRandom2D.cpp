/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "GaussianRandom2D.h"
#include <math.h>
#include <stdio.h>
#include "Ran1.h"
GaussianRandom2D::GaussianRandom2D(Eigen::MatrixXd* mean, Eigen::MatrixXd* sigma, Ran1* ran){
	//puts("Initing GaussianRandom2D");
	M = mean;
	Sigma = sigma;
	RAN = ran;
}
long GaussianRandom2D::Seed(){
	return RAN->idum;
}
void GaussianRandom2D::nextSample(double& x, double& y){
	double rannor1 = RAN->nextGasdev();
	//printf("%lf\n",rannor1);
	double rannor2 = RAN->nextGasdev();
	//printf("%lf\n",rannor2);
	double sigma1 = sqrt(Sigma->operator()(0, 0));
	x = M->operator()(0, 0) + sigma1 * rannor1;
	y = M->operator()(1, 0) + (Sigma->operator()(1, 0) * rannor1 + sqrt(Sigma->operator()(0, 0) * Sigma->operator()(1, 1) - Sigma->operator()(1, 0) * Sigma->operator()(1, 0)) * rannor2) / sigma1;
}
Eigen::VectorXd GaussianRandom2D::nextSample(){
	double rannor1 = RAN->nextGasdev();
	//printf("%lf\n",rannor1);
	double rannor2 = RAN->nextGasdev();
	//printf("%lf\n",rannor2);
	double sigma1 = sqrt(Sigma->operator()(0, 0));
	double y1 = M->operator()(0, 0) + sigma1 * rannor1;
	double y2 = M->operator()(1, 0) + (Sigma->operator()(1, 0) * rannor1 + sqrt(Sigma->operator()(0, 0) * Sigma->operator()(1, 1) - Sigma->operator()(1, 0) * Sigma->operator()(1, 0)) * rannor2) / sigma1;
	Eigen::VectorXd v(2);
	v<<y1,y2;
	return v;
}
