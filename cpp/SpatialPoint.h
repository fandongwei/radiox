/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __SpatialPoint__h__
#define __SpatialPoint__h__
class Gaussian;
#include <string>
//one point
class SpatialPoint
{
public:
	//id of the point
	long id;
	//ra and ra_err of the point
	double ra, ra_err;
	//dec and dec_err of the point
	double dec, dec_err;
	//Cartesian x,y,z coordinate of the point
	double x, y, z;
	//Zones Algorithm's alpha value
	double alpha;
public:
	SpatialPoint();
	SpatialPoint(long id, double ra, double dec, double ra_err, double dec_err);
public:
	//for sort() function
	bool operator <(const SpatialPoint& p) const;
	//return as [p,q,psigma,qsigma]
	SpatialPoint posSigma(SpatialPoint& centerobj);
	Gaussian getGaussian(SpatialPoint& centerobj, bool invert=true);
	std::string toString();
};
#endif