/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "SpatialXMatcher.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include "Constants.h"

double SpatialXMatcher::distance(double ra1, double dec1, double ra2, double dec2) {
	ra1 = Constants::Degree2Radian * ra1;
	dec1 = Constants::Degree2Radian * dec1;

	ra2 = Constants::Degree2Radian * ra2;
	dec2 = Constants::Degree2Radian * dec2;

	double rad = 2.0 * asin(sqrt(
			sin((dec1 - dec2) / 2) * sin((dec1 - dec2) / 2)
			+ cos(dec2) * cos(dec1) * sin((ra1 - ra2) / 2) * sin((ra1 - ra2) / 2)
		)
	);
	return Constants::Radian2Degree * rad;
}
void SpatialXMatcher::xyz(double ra, double dec, double& x, double& y, double& z) {
	ra = Constants::Degree2Radian * ra;
	dec = Constants::Degree2Radian * dec;

	z = sin(dec);
	double cdec = cos(dec);
	x = cdec * cos(ra);
	y = cdec * sin(ra);
}
double SpatialXMatcher::zone_alpha(double dec, double threshold) {
	if (fabs(dec) + threshold > 89.9) return 180;
	dec = dec * Constants::Degree2Radian;
	threshold = threshold * Constants::Degree2Radian;

	return Constants::Radian2Degree * (fabs(atan(
		sin(threshold) / sqrt(fabs(cos((dec - threshold)) * cos(dec + threshold)))
	)));
}
double SpatialXMatcher::distance_xyz(double ra1, double dec1, double ra2, double dec2) {
	double x1, y1, z1;
	xyz(ra1, dec1, x1, y1, z1);
	double x2, y2, z2;
	xyz(ra2, dec2, x2, y2, z2);
	double m2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
	double rad = 2.0 * asin(0.5 * sqrt(m2));
	//rad = acos(1.0-0.5*m2);//resolved by Tangent_half-angle_formula, will suffer something like [-1.0000000000000004], so not recommended
	return Constants::Radian2Degree * rad;
}
static const double halfpi = 1.570796326794896619231321691639751442099;
double dec2theta(double dec)
{
	if (dec>halfpi) dec = halfpi;
	else if (dec<-halfpi) dec = -halfpi;
	return halfpi - dec;
}
static const double twopi = 6.283185307179586476925286766559005768394;
double ra2phi(double ra)
{
	if (ra < 0) return  twopi + ra;
	return ra;
}
int SpatialXMatcher::findStartIndex(SpatialPoint* arr, int size, double bottom)
{
	int end = size - 1;
	if (arr[0].dec >= bottom) return 0;
	else if (arr[end].dec < bottom) return end + 1;
	return binarySearch(arr, 0, end, bottom);
}
int SpatialXMatcher::binarySearch(SpatialPoint* arr, int low, int high, double bottom) {
	if (arr[0].dec >= bottom) return 0;

	int mid = (low + high) / 2;
	if (low > high)
		return -1;
	else
	{
		if (arr[mid].dec >= bottom && arr[mid - 1].dec < bottom)
			return mid;
		else if (arr[mid].dec > bottom && arr[mid - 1].dec > bottom)
			return binarySearch(arr, low, mid - 1, bottom);
		else
			return binarySearch(arr, mid + 1, high, bottom);
	}
}
std::vector<SpatialXResult> SpatialXMatcher::matchTables(SpatialPointTable& tbl1, SpatialPointTable& tbl2, double degthreshold, bool issortresult) {
	int tbl1len = tbl1.data.size();
	SpatialPoint* ptbl1data = &tbl1.data[0];
	for (int i = 0; i < tbl1len; i++) {
		xyz(ptbl1data[i].ra, ptbl1data[i].dec, ptbl1data[i].x, ptbl1data[i].y, ptbl1data[i].z);
		ptbl1data[i].alpha = zone_alpha(ptbl1data[i].dec, degthreshold);
	}
	/*for(std::vector<SpatialPoint>::iterator it=tbl1.data.begin(); it<tbl1.data.end();it++)//too slow
	{
		xyz(it->ra,it->dec, it->x, it->y, it->z);
		it->alpha = alpha(it->dec, degthreshold);
	}*/

	int tbl2len = tbl2.data.size();
	SpatialPoint* ptbl2data = &tbl2.data[0];
	for (int i = 0; i < tbl2len; i++) {
		xyz(ptbl2data[i].ra, ptbl2data[i].dec, ptbl2data[i].x, ptbl2data[i].y, ptbl2data[i].z);
	}
	/*for (std::vector<SpatialPoint>::iterator it = tbl2.data.begin(); it<tbl2.data.end(); it++)//too slow
	{
		xyz(it->ra, it->dec, it->x, it->y, it->z);
	}*/
	double threshold = degthreshold * Constants::Degree2Radian / 2.0;
	{
		threshold = sin(threshold);
		threshold = 4 * threshold * threshold;
	}
	std::vector<SpatialXResult> al;
	SpatialXResult xres;
	//for (std::vector<SpatialPoint>::iterator it = tbl1.data.begin(); it<tbl1.data.end(); it++)//too slow
	for(int i=0;i<tbl1len;i++)
	{
		SpatialPoint* p1 = &ptbl1data[i];
		double bottom = p1->dec - degthreshold;
		double top = p1->dec + degthreshold;
		double left = p1->ra - p1->alpha;
		double right = p1->ra + p1->alpha;
		bool warparound = false;
		if (left < 0)
		{
			warparound = true;
			left += 360;
		}
		else if (right > 360)
		{
			warparound = true;
			right -= 360;
		}
		std::vector<SpatialXResult> al2;
		//foreach(var p2 in tbl2.data)//3.5s
		for (int j = findStartIndex(ptbl2data, tbl2len, bottom); j < tbl2len; j++)//0.35s
		{
			SpatialPoint* p2 = &ptbl2data[j];
			//filter dec, tbl2 is sorted by dec asc, so if p2->dec > top, do not need to compare the rest
			if (p2->dec > top) break;
			if (p2->dec < bottom || p2->dec > top) continue;

			if (warparound)//filter ra by alpha
			{
				if (p2->ra > right && p2->ra < left) continue;
			}
			else
			{
				if (p2->ra < left || p2->ra > right) continue;
			}
			//check two points in the circle
			double m2 = (p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y) + (p1->z - p2->z) * (p1->z - p2->z);
			if (m2 < threshold)
			{
				//caculate distance
				//should not use the Math.Acos(1-thredhost) formula, due to the precision of double, sometimes we will get something like -1.0000000000000004.
				//-1.0000000000000004 is actually -1 acoording to the precision of double, but it's out of the valid input range of Math.Acos
				xres.dist = 2.0 * asin(0.5 * sqrt(m2)) * Constants::Radian2Degree;
				xres.p1 = *p1;
				xres.p2 = *p2;
				al2.push_back(xres);
			}
		}
		if(issortresult){
			std::sort(al2.begin(), al2.end());
		}
		al.insert(al.end(), al2.begin(), al2.end());
	}
	return al;
}
void SpatialXMatcher::outputSpatialXResult(std::vector<SpatialXResult>& res, std::string filename) {
	int reslen = res.size();
	SpatialXResult* pres = &res[0];

	//slowest
	//ofstream stream(filename);
	//for (int i = 0; i < reslen; i++) {
	//	//stream << pres[i].toString() << endl;
	//	stream << pres[i].p1.id << ',' << pres[i].p1.ra << ',' << pres[i].p1.dec << ',';
	//	stream << pres[i].p2.id << ',' << pres[i].p2.ra << ',' << pres[i].p2.dec << ',';
	//	stream << pres[i].dist << '\n';
	//}
	//stream.close();

	FILE* pf = fopen(filename.c_str(), "w");

	////fast
	//for (int i = 0; i < reslen; i++) {
	//	fprintf(pf, "%ld,%.17lf,%.17lf,%ld,%.17lf,%.17lf,%.17lf\n",
	//		pres[i].p1.id, pres[i].p1.ra, pres[i].p1.dec, pres[i].p2.id, pres[i].p2.ra, pres[i].p2.dec, pres[i].dist);
	//}
	//fclose(pf);

	////fwrite a little faster than fprintf
	//char * buf = new char[1500];
	//*buf = 0;
	//int currentlen = 0;
	//for (int i = 0, j = 0; i < reslen; i++, j++) {
	//	currentlen = strlen(buf);
	//	if (10 == j) {
	//		fwrite(buf, 1, currentlen, pf);
	//		currentlen = 0;
	//		j = 0;
	//		*buf = 0;
	//	}
	//	sprintf(buf + currentlen, "%ld,%.17lf,%.17lf,%ld,%.17lf,%.17lf,%.17lf\n",
	//		pres[i].p1.id, pres[i].p1.ra, pres[i].p1.dec, pres[i].p2.id, pres[i].p2.ra, pres[i].p2.dec, pres[i].dist);
	//}
	//currentlen = strlen(buf);
	//fwrite(buf, 1, currentlen, pf);
	//fclose(pf);
	//delete[] buf;

	//strcat is more faster then tricky buf+currentlen
	char * buf = new char[2*1024];
	char * subbuf = new char[300];
	*buf = 0;
	for (int i = 0, j = 0; i < reslen; i++, j++) {
		if (10 == j) {
			fwrite(buf, 1, strlen(buf), pf);
			j = 0;
			*buf = 0;
		}
		sprintf(subbuf, "%ld,%.17lf,%.17lf,%ld,%.17lf,%.17lf,%.17lf\n", pres[i].p1.id, pres[i].p1.ra, pres[i].p1.dec, pres[i].p2.id, pres[i].p2.ra, pres[i].p2.dec, pres[i].dist);
		//sprintf(subbuf, "%ld,%lf,%lf,%ld,%lf,%lf,%lf\n", pres[i].p1.id, pres[i].p1.ra, pres[i].p1.dec, pres[i].p2.id, pres[i].p2.ra, pres[i].p2.dec, pres[i].dist);//faster
		strcat(buf, subbuf);
	}
	fwrite(buf, 1, strlen(buf), pf);
	fclose(pf);
	delete[] buf;
	delete[] subbuf;
}
/*
importing C:/Users/Dongwei/Projects/Radio Catalog Crossmatching/data/catalogs/cstar.csv, 21845 points, time:00:00:00.0325065
importing C:/Users/Dongwei/Projects/Radio Catalog Crossmatching/data/catalogs/cstar.csv, 21845 points, time:00:00:00.0475082
start matching
matched 31867 lines, time:00:00:00.1690233
write file takes time:00:00:00.1060122
total time:00:00:00.3610502
*/