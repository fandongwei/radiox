/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "CommandParser.h"
#include <math.h>
#include <stdio.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>


using namespace std;

CommandParser::CommandParser(int argc, char* argv[]){
	for(int i=0;i<argc;i++){
		stringstream ss;
		ss.str(argv[i]);
		string key,val;
		getline(ss, key, '=');
		getline(ss, val, '=');
		if('"'==val[0] || '\''==val[0]) val=val.substr(1, val.size()-2);
		params[key]=val;
	}
}
string CommandParser::get(string name){
	return params[name];
}
bool CommandParser::getBoolean(string name){
	if(params[name].compare("0")==0) return false;
	if(params[name].compare("1")==0) return true;
	string val = params[name];
	transform (val.begin(), val.end(), val.begin(), ::tolower);
	return val.compare("true")==0;
}
int CommandParser::getInt(string name){
	return atoi(params[name].c_str());
}
double CommandParser::getDouble(string name){
	return atof(params[name].c_str());
}
void CommandParser::getCharArr(std::string name, char* arr){
	string s = params[name];
	memcpy(arr,s.c_str(),s.size());
}
void CommandParser::getIntArr(string name, int* arr){
	stringstream ss;
	ss.str(params[name]);
	string item;
	int i=0;
	while (getline(ss, item, ',')) {
    	arr[i]=atoi(item.c_str());
    	i++;
	}
}
void CommandParser::getDoubleArr(string name, double* arr){
	stringstream ss;
	ss.str(params[name]);
	string item;
	int i=0;
	while (getline(ss, item, ',')) {
    	arr[i]=atof(item.c_str());
    	i++;
	}	
}

void CommandParser::printAll(){
	for(map<string, string>::iterator it=params.begin(); it!=params.end(); it++){
		cout<<it->first<<"="<<it->second<<endl;
	}
}