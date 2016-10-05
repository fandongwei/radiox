/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __Ran1__h__
#define __Ran1__h__
class Ran1{
public:
	Ran1();
	Ran1(long seed);
	~Ran1();
public:
	double nextDouble();
	double nextDouble(double min, double max);
	double nextGasdev();
	double nextGasdev(double mean, double sigma);
private:
	void init(long seed=1);
	long IA;
	long IM;
	double AM;
	long IQ;
	long IR;
	int NTAB;
	double NDIV;
	double RNMX;
	long iy;
	long iv[32];
	int iset;
	double gset;
public:
	long idum;
};
#endif