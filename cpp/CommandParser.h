/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __CommandParser__h__
#define __CommandParser__h__
#include <map>
#include <string>
class CommandParser{
public:
	CommandParser(int argc, char* argv[]);
private:
	std::map<std::string, std::string> params;
public:
	void printAll();
	std::string get(std::string name);
	void getCharArr(std::string name, char* arr);
	bool getBoolean(std::string name);
	int getInt(std::string name);
	double getDouble(std::string name);
	void getIntArr(std::string name, int* arr);
	void getDoubleArr(std::string name, double* arr);
};
#endif