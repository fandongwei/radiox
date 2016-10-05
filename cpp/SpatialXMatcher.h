/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/

#ifndef __SpatialXMatcher__h__
#define __SpatialXMatcher__h__
#include "SpatialPoint.h"
#include "SpatialPointTable.h"
#include "SpatialXResult.h"
#include <vector>
#include <string>

//functions that assist to do the spatial xmatch
class SpatialXMatcher{
public:
	//calculate two points' distance by trigonometric functions, ra1,dec1,ra2,dec2 are all in degree, return value also in degree
	static double distance(double ra1, double dec1, double ra2, double dec2);
public:
	//convert (ra,dec) to Cartesian (x,y,z), ra & dec in degree
	static void xyz(double ra, double dec, double& x, double& y, double& z);
	//calculat the Zones Algorithm's alpha value, dec & threshold in degree
	static double zone_alpha(double dec, double threshold);
	//calculate two points' distance in Cartesian system, ra1,dec1,ra2,dec2 are all in degree, return value also in degree
	static double distance_xyz(double ra1, double dec1, double ra2, double dec2);
public:
	//in a sorted 'arr', to find the smallest value which is larger than 'bottom'
	static int findStartIndex(SpatialPoint* arr, int size, double bottom);
	//using the binary search method to locate the position of the 'findStartIndex'
	static int binarySearch(SpatialPoint* arr, int low, int high, double bottom);
	//xmatch two tables, degthreshold in degree
	static std::vector<SpatialXResult> matchTables(SpatialPointTable& tbl1, SpatialPointTable& tbl2, double degthreshold, bool issortresult=true);
	//write result to file
	static void outputSpatialXResult(std::vector<SpatialXResult>& res, std::string filename);
};
#endif