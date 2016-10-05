/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "SpatialXResultTable.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Constants.h"
using namespace std;

SpatialXResultTable::SpatialXResultTable(string filename, bool sortbydec) {
	FILE* pf = fopen(filename.c_str(), "r");
	SpatialXResult x;

	//method 1
	while (!feof(pf)) {
		fscanf(pf, "%ld,%lf,%lf,%ld,%lf,%lf,%lf\n", &x.p1.id, &x.p1.ra, &x.p1.dec, &x.p2.id, &x.p2.ra, &x.p2.dec, &x.dist);
		data.push_back(x);
		//cout<<x.toString()<<endl;
	}
	fclose(pf);
	/*
	//method 2
	//read 50MB one time, if find ',','\n' and fill SpatialXResult
	size_t buflen = 64 * 1024;// 50 * 1024 * 1024;
	char* buf = new char[buflen + 1];
	char* subbuf = new char[300];
	char* subbufcurrent = subbuf;
	fseek(pf, 0, 0);
	int numberindex = 0;
	while (!feof(pf)) {
		size_t readlen = fread(buf, sizeof(char), buflen, pf);
		char* bufend = buf + readlen;
		for (char * bufcurrent = buf; bufcurrent < bufend; bufcurrent++) {
			if (*bufcurrent == '\r') continue;

			*subbufcurrent = *bufcurrent;
			subbufcurrent++;
			if (*bufcurrent != ',' && *bufcurrent != '\n') continue;

			//find a number string
			*(subbufcurrent-1) = '\0';//terminate number string
			numberindex++;
			switch (numberindex)
			{
			case 1:
				x.p1.id = atoi(subbuf);
				break;
			case 2:
				x.p1.ra = atof(subbuf);
				break;
			case 3:
				x.p1.dec = atof(subbuf);
			case 4:
				x.p2.id = atoi(subbuf);
				break;
			case 5:
				x.p2.ra = atof(subbuf);
				break;
			case 6:
				x.p2.dec = atof(subbuf);
				break;
			case 7:
				x.dist = atof(subbuf);
				data.push_back(x);
				numberindex = 0;
				break;
			}
			//reset the number string
			subbufcurrent = subbuf;
		}
	}
	delete[] buf;
	delete[] subbuf;
	fclose(pf);
	*/
	//sort data
	if (sortbydec) {
		std::sort(data.begin(), data.end());
	}
	currentIndex = data.begin();
}
SpatialXResultTable::SpatialXResultTable(std::vector<SpatialXResult>& d, bool sortbydec) {
	//copy method 1
	//int len = d.size();
	//SpatialXResult* ps = &d[0];
	//for (int i = 0; i < len; i++) {
	//	data.push_back(ps[i]);
	//}

	//copy method 2
	data = d;
	
	if (sortbydec) {
		sort(data.begin(), data.end());
	}
	currentIndex = data.begin();
}
SpatialXResultTable::~SpatialXResultTable() {
	data.clear();
}

vector<SpatialXResult> SpatialXResultTable::getNextCombination(){
	vector<SpatialXResult> radios;
	if(currentIndex == data.end()) return radios;
	
	long oid = currentIndex->p1.id;
	do
	{
		radios.push_back(*currentIndex);
		currentIndex++;
	}
	while(currentIndex < data.end() && oid==currentIndex->p1.id);
	return radios;
}