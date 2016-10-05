/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __Det__h__
#define __Det__h__
#include "Gaussian.h"
#include <Eigen/Dense>
class Det : public Gaussian {
public:
	Det(Eigen::MatrixXd pos, Eigen::MatrixXd mat, bool invert);
	Det(double x, double y, double sigma);
};
#endif