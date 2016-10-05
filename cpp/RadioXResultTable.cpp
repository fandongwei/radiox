/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "RadioXResultTable.h"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Constants.h"
using namespace std;

RadioXResultTable::RadioXResultTable(string filename, bool sortbylikelihood) {
	FILE* pf = fopen(filename.c_str(), "r");
	RadioXResult x;
	
	/*
	//method 1, not work
	char buf0[255], buf1[255];
	while (!feof(pf)) {
		fscanf(pf, "%ld,%d,%s,%s,%lf,%lf\n", &x.opticalid, &x.combitype, buf0, buf1, &x.likelihood, &x.likelihood_err);
		x.combi = string(buf0);
		x.radioids = string(buf1);
		if(x.combitype!=2 && x.combitype!=0){
			data.push_back(x);
		}
	}
	fclose(pf);
	*/
	//method 2
	//read 50MB one time, if find ',','\n' and fill RadioXResult
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
			//printf("subbuf=%s\t", subbuf);
			switch (numberindex)
			{
			case 1:
				x.opticalid = atoi(subbuf);
				break;
			case 2:
				x.combitype = atoi(subbuf);
				break;
			case 3:
				x.combi = subbuf;
			case 4:
				x.radioids = subbuf;
				break;
			case 5:
				x.likelihood = atof(subbuf);
				break;
			case 6:
				x.likelihood_err = atof(subbuf);
				//if(x.combitype!=2 && x.combitype!=0){//remove the 'lobe' && 'none' type data
					data.push_back(x);
				//}
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
	//sort data
	init(sortbylikelihood);
}
RadioXResultTable::RadioXResultTable(std::vector<RadioXResult>& d, bool sortbylikelihood) {
	//copy method 1
	//int len = d.size();
	//RadioXResult* ps = &d[0];
	//for (int i = 0; i < len; i++) {
	//	data.push_back(ps[i]);
	//}

	//copy method 2
	data = d;
	init(sortbylikelihood);
}
void RadioXResultTable::init(bool sortbylikelihood){
	if (sortbylikelihood) {
		std::sort(data.begin(), data.end());
	}
	//cout<<"init corelikelihood"<<endl;
	for(vector<RadioXResult>::iterator it = data.begin(); it!=data.end(); it++){
		{
			//it->components.push_back(it->opticalid);//wrong
			stringstream ss;
			ss.str(it->radioids);
			string item;
			while (getline(ss, item, '|')) {
	        	it->radiocomponents.push_back(atoi(item.c_str()));
	    	}
	    	//if(it->components.size()>2) printf(">3:%s\n", it->radioids.c_str());
    	}
		if(it->combitype == 1){
			if(corelikelihood.find(it->opticalid)==corelikelihood.end() || corelikelihood[it->radiocomponents[0]] < it->likelihood){
			//if(corelikelihood[it->opticalid] < it->likelihood){//if item doesn't exists corelikelihood[item]==0, it will be a risk at here.... probably
				corelikelihood[it->radiocomponents[0]] = it->likelihood;//save the largest core-likelihood
			}
		}
	}
	//cout<<"init corelikelihood done!"<<endl;
}
RadioXResultTable::~RadioXResultTable() {
	data.clear();
	corelikelihood.clear();
}

bool RadioXResultTable::nextTopRadioXResult(RadioXResult& rxr){
	if(data.size()==0) return false;
	//rxr = this->data[0];//top
	rxr = *(data.end()-1);//bottom
	return true;
}
void RadioXResultTable::removeTopRadioxResult(){
	//this->data.erase(this->data.begin());
	if(data.size()==0) return;
	this->data.erase(this->data.end()-1);
}
/*
void RadioXResultTable::removeRelatedRadioXResult(RadioXResult& rxr){
	//cout<<"removing "<<rxr.toString()<<endl;
	if(rxr.radiocomponents.size()==0) return;
	for(vector<RadioXResult>::iterator it = data.begin(); it!=data.end(); ){
		//cout<<"it->component size="<<it->components.size()<<", rxr.component size="<<rxr.components.size()<<endl;
		bool iserase = false;
  		for(vector<long>::iterator it2 = it->radiocomponents.begin(); it2 != it->radiocomponents.end(); it2++){
  			//cout<<"looking for "<<*it2<<endl;
  			if (find (rxr.radiocomponents.begin(), rxr.radiocomponents.end(), *it2) != rxr.radiocomponents.end() || rxr.opticalid == *it2){
  				iserase = true;
  				break;
  			}
  		}
  		if (find (rxr.radiocomponents.begin(), rxr.radiocomponents.end(), it->opticalid) != rxr.radiocomponents.end() || rxr.opticalid == it->opticalid){
  			iserase = true;
  		}
  		if(iserase){
  			//cout<<"erase "<<it->toString()<<endl;
  			data.erase(it);
  		}else{
  			 it++;
  		}
	}
}
*/
bool RadioXResultTable::nextTopRadioXResult(RadioXResult& rxr, std::map<long, int>& opticalremovedcomponents, std::map<long, int>& radioremovedcomponents){
	//if contain removed components, erase current top item, and get new top 1
	//vector<RadioXResult>::iterator it=data.begin();
	if(data.size()==0) return false;
	vector<RadioXResult>::iterator it=data.end()-1;
	//while(it!=data.end()){
	while(it!=data.begin()){
		//clock_t t = clock();
		bool iserase = false;
		std::map<long, int>::iterator it3 = opticalremovedcomponents.find(it->opticalid);
		if(it3 != opticalremovedcomponents.end()) iserase = 666==it3->second;
		//clock_t nt = clock(); cout<<"optical removedcomponents detection:"<<float(nt-t)/CLOCKS_PER_SEC<<endl;t = nt;
		if(!iserase){
	  		for(vector<long>::iterator it2 = it->radiocomponents.begin(); it2 != it->radiocomponents.end(); it2++){
	  			//cout<<"looking for "<<*it2<<endl;
	  			it3 = radioremovedcomponents.find(*it2);
	  			if(it3 != radioremovedcomponents.end()){
	  				iserase = 666==it3->second;
	  				if (iserase) break;
	  			}
	  		}
  		}
  		
  		//nt = clock(); cout<<"radio removedcomponents detection:"<<float(nt-t)/CLOCKS_PER_SEC<<endl;
  		if(iserase){
  			//cout<<"erase "<<it->toString()<<endl;
  			data.erase(it);
  		}else{
  			rxr=*it;
			return true;
  		}
		it = data.end()-1;
  	}
  	return false;
}
vector<RadioXResult> RadioXResultTable::getResultsByType(int type){
	vector<RadioXResult> rrv;
	for(vector<RadioXResult>::iterator it = data.begin(); it!=data.end(); it++){
		if(type==it->combitype)
			rrv.push_back(*it);
	}	
	return rrv;
}
vector<RadioXResult> RadioXResultTable::getTriplets(){
	return getResultsByType(1|2|4);
}
vector<RadioXResult> RadioXResultTable::getDoublets(){
	vector<RadioXResult> rrv;
	for(vector<RadioXResult>::iterator it = data.begin(); it!=data.end(); it++){
		if(3==it->combitype || 6==it->combitype)
			rrv.push_back(*it);
	}	
	return rrv;
}
vector<RadioXResult> RadioXResultTable::getCores(){
	return getResultsByType(1);
}