/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "RadioXRunner.h"
#include <math.h>
#include <stdio.h>
#include <vector>
#include <string>
#include "Hypothesis.h"
#include "LobePrior.h"
#include "RadioXResult.h"
#include <iostream>
#include <time.h>

using namespace std;

void RadioXRunner::doRadioX(bool sym, LobePrior* pri, double* sigmas, SpatialXResultTable& tbl, std::string outfile, int sampleCounts[]){
	vector<SpatialXResult> com;
	
	FILE* pf = fopen(outfile.c_str(), "w");
	char buf[500];
	//clock_t t1 = clock(), t2;
	while((com=tbl.getNextCombination()).size()>0){
		//t2 = clock(); printf("next combination, size=%d, takes %f seconds\n", com.size(), float(t2-t1)/CLOCKS_PER_SEC); t1=t2;
		SpatialPoint* optical = &(com[0].p1);
		optical->ra_err = sigmas[0];
		optical->dec_err = sigmas[1];
		vector<SpatialPoint*> radios;
		for(vector<SpatialXResult>::iterator it=com.begin(); it!=com.end(); it++){
			it->p2.ra_err = sigmas[0];
			it->p2.dec_err = sigmas[1];
			radios.push_back(&(it->p2));
		}
		//t2 = clock(); printf("\tinited sigmas, takes %f seconds\n", float(t2-t1)/CLOCKS_PER_SEC); t1=t2;
		vector<Hypothesis> hyps = Hypothesis::generateHypList(optical, radios);
		//t2 = clock(); printf("\tgenerated hyps, size=%d, takes %f seconds\n", hyps.size(), float(t2-t1)/CLOCKS_PER_SEC); t1=t2;
		if(sym){
			for(std::vector<Hypothesis>::iterator it=hyps.begin();it!=hyps.end();it++){
				double likelihood, likelihood_err;
				it->likelihoodSym(pri, sampleCounts, likelihood, likelihood_err);
				sprintf(buf, ",%lf,%lf\n", likelihood, likelihood_err);
				string str = it->combiStr() + buf;
				fwrite(str.c_str(), 1, str.size(), pf);
			}
		}else{
			for(std::vector<Hypothesis>::iterator it=hyps.begin();it!=hyps.end();it++){
				//it->print();
				double likelihood, likelihood_err;
				it->likelihoodAsym(pri, sampleCounts, sigmas[2], likelihood, likelihood_err);
				sprintf(buf, ",%lf,%lf\n", likelihood, likelihood_err);
				string str = it->combiStr() + buf;
				fwrite(str.c_str(), 1, str.size(), pf);
			}
		}
		//fflush(pf);
		//t2 = clock(); printf("calculated %d hyps for combination size: %d, takes %f seconds\n", hyps.size(), com.size(), float(t2-t1)/CLOCKS_PER_SEC); t1=t2;
	}
	
	fclose(pf);
}

void RadioXRunner::doRadioX(bool sym, LobePrior* pri, double* sigmas, string infile, string outfile, int sampleCounts[]){
	SpatialXResultTable tbl(infile, false);
	doRadioX(sym, pri, sigmas, tbl, outfile, sampleCounts);
}
/*
void RadioXRunner::analysisResult0(RadioXResultTable& tbl, string outfile){
	FILE* pf = fopen(outfile.c_str(), "w");
	RadioXResult rxr;
	
	while(tbl.nextTopRadioXResult(rxr)){
		double smf = 0;
		for(vector<string>::iterator it = rxr.components.begin();it!=rxr.components.end();it++){
			smf += tbl.corelikelihood[*it];
		}
		if(rxr.likelihood<smf){
			tbl.removeTopRadioxResult();
		}else{
			string str = rxr.toString()+"\n";
			fwrite(str.c_str(), 1, str.size(), pf);
			tbl.removeRelatedRadioXResult(rxr);//slow part
		}
		fflush(pf);
	};
	fclose(pf);
}
*/
void RadioXRunner::analysisResult(RadioXResultTable& tbl, string outfile){
	//cout<<"radiox result lines:"<<tbl.data.size()<<endl;
	FILE* pf = fopen(outfile.c_str(), "w");
	RadioXResult rxr;
	map<long, int> opticalremovedcomponents;
	map<long, int> radioremovedcomponents;
	//cout<<"start analysing"<<endl;
	//clock_t totalt = clock(), t = clock(), nt;
	while(tbl.nextTopRadioXResult(rxr, opticalremovedcomponents, radioremovedcomponents)){
		//nt = clock(); cout<<"next combination : "<<float(nt-t)/CLOCKS_PER_SEC<<endl; t = nt;
		double smf = 0;
		for(vector<long>::iterator it = rxr.radiocomponents.begin();it!=rxr.radiocomponents.end();it++){
			smf += tbl.corelikelihood[*it];
		}
		//nt = clock(); cout<<"corelikelihood = "<<smf<<" : "<<float(nt-t)/CLOCKS_PER_SEC<<endl; t = nt;
		//cout<<"rxr.likelihood="<<rxr.likelihood<<", smf="<<smf<<endl;
		if(rxr.likelihood>=smf){
			string str = rxr.toString()+"\n";
			//cout<<"pick "<<str;
			fwrite(str.c_str(), 1, str.size(), pf);
			opticalremovedcomponents[rxr.opticalid]=666;
			for(vector<long>::iterator itx=rxr.radiocomponents.begin();itx!=rxr.radiocomponents.end();itx++){
				radioremovedcomponents[*itx]=666;
			}
			//fflush(pf);
		}
		tbl.removeTopRadioxResult();
		//nt = clock(); cout<<"erased : "<<float(nt-t)/CLOCKS_PER_SEC<<endl; t = nt;
	}
	//cout<<"finish analysing!"<<endl;
	fclose(pf);
}

void RadioXRunner::analysisResult(string infile, string outfile){
	RadioXResultTable tbl(infile, true);
	analysisResult(tbl, outfile);
}
