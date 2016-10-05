/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __Constants__h__
#define __Constants__h__
//useful constants
class Constants {
public:
	//PI
	static const double PI;
	//multiply with this number to convert degree to radian
	static const double Degree2Radian;
	static const double Arcmin2Radian;
	static const double Arcsec2Radian;
	//multiply with this number to convert radian to degree
	static const double Radian2Degree;
	static const double Radian2Arcmin;
	//multiply with this number to convert radian to arc second
	static const double Radian2ArcSecond;
	static const double SquareRadian2SquareDegree;
	static const double WholeSphereInSquareDegree;
	static const double WholeSphereInSquareArcsec;
};
#endif