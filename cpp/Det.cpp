/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "Det.h"
#include <math.h>
#include <stdio.h>
Det::Det(Eigen::MatrixXd pos, Eigen::MatrixXd mat, bool invert)
	:Gaussian(pos, mat, invert){
	
}
Det::Det(double x, double y, double sigma)
	:Gaussian(x,y,sigma){
	//puts("initted Det");
}