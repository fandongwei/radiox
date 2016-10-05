/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "SpatialXResult.h"
#include<stdio.h>
#include <iostream>
#include <sstream>

std::string SpatialXResult::toString() {
	char buf[500];
	sprintf(buf, "%ld,%.17lf,%.17lf,%ld,%.17lf,%.17lf,%.17lf", p1.id, p1.ra, p1.dec, p2.id, p2.ra, p2.dec, dist);
	return std::string(buf);
	//std::stringstream ss;
	//ss << p1.id << ',' << p1.ra << ',' << p1.dec << ',' << p2.id << ',' << p2.ra << ',' << p2.dec << ',' << dist;
	//return ss.str();
}
bool SpatialXResult::operator <(const SpatialXResult& x) const{
	//return this->id < p.id;
	if (p1.dec == x.p1.dec){
		if(p1.ra == x.p1.ra){
			return dist < x.dist;
		}
		return p1.ra < x.p1.ra;
	}
	else return p1.dec < x.p1.dec;
}