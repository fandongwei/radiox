/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __testSpatialX__h__
#define __testSpatialX__h__
#include "SpatialXMatcher.h"
#include "Constants.h"
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <math.h>
using namespace std;

void testConstants() {
	cout << "--------testConstants" << endl;
	Constants s;
	printf("PI=%.17lf \ndeg2rad=%.17lf \nrad2deg=%.17lf\n", s.PI, s.Degree2Radian, s.Radian2Degree);
	cout << endl;
}
void test2Points()
{
	cout << "--------test two points" << endl;
	double ra1 = 10.1111, dec1 = 20.22222;
	double ra2 = 30.33333, dec2 = 40.44444;

	double dist = SpatialXMatcher::distance(ra1, dec1, ra2, dec2);
	printf("dist =%.17lf\n", dist);

	dist = SpatialXMatcher::distance_xyz(ra1, dec1, ra2, dec2);
	printf("dist2=%.17lf\n", dist);

	SpatialXResult r;
	r.dist = 0.111;
	r.p1.id = 0;
	r.p1.alpha = 0;
	r.p1.ra = 0;
	r.p1.dec = 0;
	r.p1.x = 0;
	r.p1.y = 0;
	r.p1.z = 0;
	r.p2.id = 0;
	r.p2.alpha = 0;
	r.p2.ra = 0;
	r.p2.dec = 0;
	r.p2.x = 0;
	r.p2.y = 0;
	r.p2.z = 0;
	cout << r.toString() << endl;
	cout << endl;
}
void testSpatialPointTable1() {
	cout << "--------test Data Table 1" << endl;
	vector<SpatialPoint> data;
	SpatialPoint p;
	{
		p.dec = 0; p.ra = 0;
		data.push_back(p);
		p.dec = -10; p.ra = -10;
		data.push_back(p);
		p.dec = -10; p.ra = 0;
		data.push_back(p);
		p.dec = 0; p.ra = -10;
		data.push_back(p);
	}
	SpatialPointTable d(data, true);
	printf("size %d\n", d.data.size());
	for (vector<SpatialPoint>::iterator it = d.data.begin(); it < d.data.end(); it++) {
		printf("%lf\n", it->ra);
	}
	cout << endl;
}
void testSpatialPointTable2() {
	cout << "--------test Data Table 2" << endl;
	clock_t t = clock();
	SpatialPointTable d("../data/cstar/cstar.csv", true);
	cout <<"time:"<<float(clock()-t)/CLOCKS_PER_SEC<<" s"<< endl;
	cout <<"rows:"<<d.data.size();
}
void testMatch2Tables() {
	cout << "--------test xmatch two tables" << endl;
	clock_t totalt = clock();
	clock_t t = clock();
	SpatialPointTable d1("../data/cstar/cstar.csv", true);
	printf("import  %d lines from table1 in %fs\n", d1.data.size(), float(clock() - t) / CLOCKS_PER_SEC);
	t = clock();
	SpatialPointTable d2("../data/cstar/cstar.csv", true);
	printf("import  %d lines and sort from table2 in %fs\n", d2.data.size(), float(clock() - t) / CLOCKS_PER_SEC);
	t = clock();
	vector<SpatialXResult> res = SpatialXMatcher::matchTables(d1, d2, 1.0 / 60.0,  true);
	printf("matched %ld pairs in %fs\n", res.size(), float(clock() - t) / CLOCKS_PER_SEC);
	t = clock();
	SpatialXMatcher::outputSpatialXResult(res, "../data/cstar_res_cpp.csv");
	printf("write file in %fs\n", float(clock() - t) / CLOCKS_PER_SEC);
	printf("time performance in total: %f\n", float(clock() - totalt) / CLOCKS_PER_SEC);
	cout << endl;
}

void testSpatialX() {
	testConstants();
	test2Points();
	testSpatialPointTable1();
	testSpatialPointTable2();
	testMatch2Tables();
}

#endif