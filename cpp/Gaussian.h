/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __Gaussian__h__
#define __Gaussian__h__
#include <Eigen/Dense>
class Gaussian {
public:
	Eigen::MatrixXd M;
	Eigen::MatrixXd F;
	double N;
	double S;
	int D;
public:
	Gaussian();
	Gaussian(Eigen::VectorXd m, Eigen::MatrixXd V, bool invert=true);
	Gaussian(double x, double y, double sigma, bool invert=true);
private:
	void init(Eigen::VectorXd m, Eigen::MatrixXd V, bool invert=true);
public:
	static double sandwich(Eigen::MatrixXd a, Eigen::MatrixXd M, Eigen::MatrixXd b);
	double evaluate(Eigen::MatrixXd x);
	Gaussian multiply(Gaussian g, bool normalize);
public://Det
	double likelihood(Eigen::MatrixXd model);
	double sigma();
	double sigmaSquare();
	double x();
	double y();
};
#endif