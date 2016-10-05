/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __LobePrior__h__
#define __LobePrior__h__
class LobePrior{
	public:
		virtual double prior(double x2){return -1;}
		virtual double variance(){return -1;}
};
#endif