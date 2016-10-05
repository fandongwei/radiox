/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "Constants.h"
#include "SpatialPoint.h"
#include <math.h>
#include <algorithm>
#include <string>
#include <Eigen/Dense>
#include "Gaussian.h"
#include <stdio.h>
SpatialPoint::SpatialPoint(){
}
SpatialPoint::SpatialPoint(long id, double ra, double dec, double ra_err, double dec_err){
	this->id = id;
	this->ra = ra;
	this->dec = dec;
	this->ra_err = ra_err;
	this->dec_err = dec_err;
}
bool SpatialPoint::operator <(const SpatialPoint& p) const{
	//return this->id < p.id;
	if (dec == p.dec) return ra < p.ra;
	else return dec < p.dec;
}
SpatialPoint SpatialPoint::posSigma(SpatialPoint& centerobj) {
	double cra = centerobj.ra * Constants::Degree2Radian;
	double cdec = centerobj.dec * Constants::Degree2Radian;
	double sinRA = sin(cra);
	double cosRA = cos(cra);
	double sinDec = sin(cdec);
	double cosDec = cos(cdec);
	///////////
	double tra = ra * Constants::Degree2Radian;
	double tdec = dec * Constants::Degree2Radian;
	double sinRAp = sin(tra);
	double cosRAp = cos(tra);
	double sinDecp = sin(tdec);
	double cosDecp = cos(tdec);
	///////////
	double p = Constants::Radian2ArcSecond * (
		-sinDec * cosRA * (cosDecp * cosRAp - cosDec * cosRA)
		- sinDec * sinRA * (cosDecp * sinRAp - cosDec * sinRA)
		+ cosDec * (sinDecp - sinDec)
		);
	double q = Constants::Radian2ArcSecond * (
		sinRA * (cosDecp * cosRAp - cosDec * cosRA)
		- cosRA * (cosDecp * sinRAp - cosDec * sinRA)
		);
	//////////
	double era2 = ra_err * ra_err;
	double edec2 = dec_err*dec_err;
	/////////
	double pgRA = sinDec * cosRA * cosDecp * sinRAp - sinDec * sinRA * cosDecp * cosRAp;
	double pgDec = sinDec * cosRA * cosRAp * sinDecp + sinDec * sinRA * sinRAp * sinDecp + cosDec * cosDecp;
	double psigma = pgRA * pgRA * era2 + pgDec * pgDec * edec2;
	///////
	double qgRA = -sinRA * cosDecp * sinRAp - cosRA * cosDecp * cosRAp;
	double qgDec = -sinRA * cosRAp * sinDecp + cosRA * sinRAp * sinDecp;
	double qsigma = qgRA * qgRA * era2 + qgDec * qgDec * edec2;
	///////
	SpatialPoint res;
	res.ra = p;
	res.ra_err = psigma;
	res.dec = q;
	res.dec_err = qsigma;
	///////
	//double res[4];
	//res[0] = p, res[1] = q, res[2] = psigma, res[3] = qsigma;
	return res;
}
Gaussian SpatialPoint::getGaussian(SpatialPoint& centerobj, bool invert){
	SpatialPoint p = posSigma(centerobj);
	Eigen::VectorXd pos(2);
	pos<<p.ra,p.dec;
	Eigen::MatrixXd s(2,2);
	s.setZero();
	s(0,0) = p.ra_err;
	s(1,1) = p.dec_err;
	return Gaussian(pos, s, invert);
}
std::string SpatialPoint::toString(){
	char buf[255];
	sprintf(buf, "%ld(%lf %lf)", id, ra, dec);
	return std::string(buf);
}