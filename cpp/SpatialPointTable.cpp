/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "SpatialPointTable.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "Constants.h"

SpatialPointTable::SpatialPointTable(std::string filename, bool sortbydec) {
	FILE* pf = fopen(filename.c_str(), "r");
	SpatialPoint p;

	
	//method 1
	while (!feof(pf)) {
		fscanf(pf, "%ld,%lf,%lf\n", &p.id, &p.ra, &p.dec);
		data.push_back(p);
	};
	fclose(pf);
	
	/*
	//method 2
	//read 50MB one time, if find ',','\n' and fill SpatialPoint
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
			else *(subbufcurrent-1) = 0;//terminate number string

			numberindex++;
			switch (numberindex)
			{
			case 1:
				p.id = atoi(subbuf);
				break;
			case 2:
				p.ra = atof(subbuf);
				break;
			case 3:
				p.dec = atof(subbuf);
				data.push_back(p);
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
}
SpatialPointTable::SpatialPointTable(std::vector<SpatialPoint>& d, bool sortbydec) {
	//copy method 1
	//int len = d.size();
	//SpatialPoint* ps = &d[0];
	//for (int i = 0; i < len; i++) {
	//	data.push_back(ps[i]);
	//}

	//copy method 2
	data = d;
	
	if (sortbydec) {
		std::sort(data.begin(), data.end());
	}
}
SpatialPointTable::~SpatialPointTable() {
	data.clear();
}
