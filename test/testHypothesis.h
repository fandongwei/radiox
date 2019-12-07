/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#ifndef __testHypothesis__h__
#define __testHypothesis__h__
#include "Hypothesis.h"
#include "Gaussian.h"
#include "Ran1.h"
#include "GaussianRandom2D.h"
#include "Det.h"
#include "LobePrior.h"
#include "LobePrior_Uniform.h"
#include "LobePrior_Rayleigh.h"
#include "LobePrior_LogNormal.h"
#include "SpatialXResultTable.h"
#include "RadioXResult.h"
#include "RadioXResultTable.h"
#include "RadioXRunner.h"
#include <map>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

void testHypsCombination() {
	cout << "--------test hypotheses combinations" << endl;
	SpatialPoint op;
	op.id = 1;
	vector<SpatialPoint*> rs;
	SpatialPoint r1,r2,r3;
	r1.id=101; rs.push_back(&r1);
	r2.id=102; rs.push_back(&r2);
	r3.id=103; rs.push_back(&r2);
	vector<Hypothesis> hs = Hypothesis::generateHypList(&op, rs);
	int hypssize = hs.size();
	printf("size=%d\n", hypssize);
	for(int i=0;i<hypssize;i++){
		hs[i].print();
		printf(" %ld\n", hs[i].rcore->id);
	}
	cout << endl;
}
void testMatrix(){
	cout<< "--------test Matrix" << endl; 
	MatrixXd m(2,2);
	cout << "m = \n" << m << endl;
	m.setZero();
	cout << "m.setZero = \n" << m << endl;
	m(0, 0) = 3;
	m(0, 1) = 2.5;
	m(1, 0) = -1;
	m(1, 1) = m(1,0) + m(0,1);
	cout << "m = \n" << m << endl;
	cout << "m.inverse = \n" << m.inverse() << endl;
	cout << "m * m.inverse = \n" << m * m.inverse() << endl;
	cout << "m.norm = " << m.norm() << endl;
	cout << "m.norm 2= " << sqrt(m(0,0)*m(0,0)+m(0,1)*m(0,1)+m(1,0)*m(1,0)+m(1,1)*m(1,1)) << endl;
	cout << "m.determinant = " << m.determinant() << endl;
	cout << "m.transpose = \n" << m.transpose() << endl;
	cout << "m.rows = " << m.rows() << endl;
	cout << "m.cols = " << m.cols() << endl;
	VectorXd v(2);
	v << 2, 3;
	cout << "v =" << endl << v << endl;
	cout << "v.transpose =" << endl << v.transpose() << endl;
	cout << "m * v =" << endl << m * v << endl;
	cout << "v.norm = " << v.norm()<<endl;
	cout << "v.squaredNorm = " << v.squaredNorm()<<endl;
}
void testPGSigma(){
	cout<< "--------test PGSigma" << endl; 
	SpatialPoint optical(0, 51.87185, -27.73129, 0.0301, 0.0249);//"SWIRE3_J032729.24-274352.6"
	SpatialPoint c0(85, 51.8700249989828, -27.7393055555556, 0.26, 0.62);//C085
	SpatialPoint c1(87, 51.8718583345413, -27.7312777773539, 0.32, 0.85);//C087
	SpatialPoint c2(88, 51.8749375025431, -27.7223333332274, 0.14, 0.23);//C088
	
	SpatialPoint res = optical.posSigma(optical);
	printf("%lf,%lf,%lf,%lf\n", res.ra, res.dec, res.ra_err, res.dec_err);
	res = c0.posSigma(optical);
	printf("%lf,%lf,%lf,%lf\n", res.ra, res.dec, res.ra_err, res.dec_err);
	res = c1.posSigma(optical);
	printf("%lf,%lf,%lf,%lf\n", res.ra, res.dec, res.ra_err, res.dec_err);
	res = c2.posSigma(optical);
	printf("%lf,%lf,%lf,%lf\n", res.ra, res.dec, res.ra_err, res.dec_err);
	/*
	0.000000,0.000000,0.000620,0.000710
	-28.856043,5.814943,0.384400,0.052955
	0.044002,-0.026558,0.722500,0.080228
	32.243877,-9.839144,0.052900,0.015359
	*/
}
void testGuassian(){
	cout<< "--------test Gaussian" << endl; 
	Gaussian op(0.01, 0.01, 0.2, true);
	Gaussian g0(0.1, 0.1, 0.2, true);
	VectorXd mx(2);
	mx<< -0.033448604850929028, -0.033448604850929028;
	double v0 = g0.evaluate(mx) * op.evaluate(mx);
	Gaussian g1 = g0.multiply(op, false);
	double v1 = g1.evaluate(mx);
	printf("v0=%.14lf \nv1=%.14lf \ng1.S=%.14lf\n", v0, v1, g1.S);
	//v1=9.675430 v2=9.675430
}
void testRandom(){
	cout<< "--------test Random" << endl; 
	Ran1 ran;
	for(int i=0;i<100;i++){
		//ran.nextDouble();
		printf("%.16lf,", ran.nextDouble());
	}
	
	puts("\n\n\n");
	Det o(0, 0, 0.4);
	
	MatrixXd ofi = o.F.inverse();
	GaussianRandom2D gr(&o.M, &ofi, &ran);
	
	//Eigen::MatrixXd * pm = &o.M;
	//printf("pm->(0,0)=%lf\n", (*pm)(0,0));
	
	for(int i=0;i<100;i++){
		Eigen::VectorXd v = gr.nextSample();
		printf("%.16lf, %.16lf\t%.16lf\n", v(0,0), v(1,0), o.evaluate(v));
	}
}
void testLikelihood(){
	cout<< "--------test Likelihood" << endl; 
	clock_t ttotal = clock(), t1 = clock(), t2;
	SpatialPoint optical(0, 51.87185, -27.73129, 0.0301, 0.0249);//"SWIRE3_J032729.24-274352.6"
	SpatialPoint c0(85, 51.8700249989828	,	-27.7393055555556	,	0.26	,	0.62);//"C085"
	SpatialPoint c1(87, 51.8718583345413	,	-27.7312777773539	,	0.32	,	0.85);//"C087"
	SpatialPoint c2(88, 51.8749375025431	,	-27.7223333332274	,	0.14	,	0.23);//"C088"
	vector<SpatialPoint*> radios;
	radios.push_back(&c0);
	radios.push_back(&c1);
	radios.push_back(&c2);
	t2 = clock(); printf("inited data, takes %f second\n",float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	
	/////////////////////////////////////////////////////
	Gaussian g = optical.getGaussian(optical); printf("%.14lf, %.14lf, %.14lf\n", g.x(), g.y(), g.sigma());
	g = c0.getGaussian(optical); printf("%.14lf, %.14lf, %.14lf\n", g.x(), g.y(), g.sigma());
	g = c1.getGaussian(optical); printf("%.14lf, %.14lf, %.14lf\n", g.x(), g.y(), g.sigma());
	g = c2.getGaussian(optical); printf("%.14lf, %.14lf, %.14lf\n", g.x(), g.y(), g.sigma());
	
	t2 = clock(); printf("generated Gaussian, takes %f second\n",float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	/////////////////////////////////////////////////////
	
	LobePrior* pri = new LobePrior_Rayleigh(44);
	printf("prior variance = %.14lf\n", pri->variance());
	
	int sampleCount=2000;
	printf("sampleCount = %d\n", sampleCount);
	
	
	vector<Hypothesis> hyps = Hypothesis::generateHypList(&optical, radios);
	int hypssize = hyps.size();
	t2 = clock(); printf("generated %d hyps, takes %f second\n", hypssize, float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	
	for(std::vector<Hypothesis>::iterator it=hyps.begin();it!=hyps.end();it++){
		it->print();
		double likelihood, likelihood_err;
		//it->likelihoodSym(pri, sampleCounts, likelihood, likelihood_err);
		it->likelihoodAsym(pri, sampleCount, 0.2, likelihood, likelihood_err);
		printf("\t%.4lf\t%lf", likelihood, likelihood_err);
		t2 = clock(); printf("\t takes %f second\n", float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	}
	
	/*
	Hypothesis hyp;
	hyp.optical = &optical;
	hyp.rcore = &c1;
	//hyp.rlobes1 = &c1;
	hyp.rlobes1 = &c0;
	hyp.rlobes2 = &c2;
	Eigen::VectorXd v = hyp.likelihoodSym(pri, sampleCounts);
	//Eigen::VectorXd v = hyp.likelihoodAsym(pri, sampleCounts, 0.2);
	hyp.print();
	printf("\t%.4lf\t%lf\n", v(0,0), v(1,0));
	*/
	delete pri;
	printf("totally takes: %f seconds\n", float(clock()-ttotal)/CLOCKS_PER_SEC);
	/*
	--------test Likelihood
0.00000000000000, 0.00000000000000, 0.02490000000000
-28.85604299963637, 5.81494290559149, 0.61999999387413
0.04400152506837, -0.02655804127065, 0.84999999999998
32.24387689160699, -9.83914356563251, 0.22999999713845
prior variance = 830.93831132508024
sampleCounts = 80, 50, 180
hyps count=19
None None None  0.0000  0.000000
None Lobe None  7.2730  0.000089
None Lobe Lobe  -16.2148        0.431562
None Lobe Lobe  2.5215  0.067948
None Lobe None  5.7832  0.004826
None Lobe Lobe  -165.5097       0.000000
None Lobe None  7.3018  0.000026
Core None None  -inf    0.000000
Core Lobe None  -inf    nan
Core Lobe Lobe  -inf    nan
Core Lobe None  -inf    nan
Core None None  11.5437 0.000000
Core Lobe None  18.8167 0.000089
Core Lobe Lobe  14.0671 0.067734
Core Lobe None  18.8455 0.000026
Core None None  -inf    0.000000
Core Lobe None  -inf    nan
Core Lobe Lobe  -inf    nan
Core Lobe None  -inf    nan
	*/
	/*C# output 
	SWIRE3_J032729.24-274352.6 C085 C087 C088
0,0,0.0249
-28.8560429996333,5.81494290560204,0.619999993874127
0.0440015250773424,-0.0265580412642843,0.849999999999979
32.2438768916012,-9.83914356560657,0.229999997138455
prior variance = 830.93831132508
sampleCounts = 80, 50, 180
hyps count=19
 None None None :       0.000         0.000      takes 00:00:00.0030034
 None None Core :          -¡Þ         0.000     takes 00:00:00.0020002
 None None Lobe :       7.302         0.000      takes 00:00:00.0090047
 None Core None :      11.544         0.000      takes 00:00:00
 None Core Lobe :      18.846         0.000      takes 00:00:00.0060043
 None Lobe None :       5.783         0.005      takes 00:00:00.0070049
 None Lobe Core :          -¡Þ           NaN     takes 00:00:00.0060047
 None Lobe Lobe :    -165.510         0.000      takes 00:00:01.2209308
 Core None None :          -¡Þ         0.000     takes 00:00:00
 Core None Lobe :          -¡Þ           NaN     takes 00:00:00.0050022
 Core Lobe None :          -¡Þ           NaN     takes 00:00:00.0040024
 Core Lobe Lobe :          -¡Þ           NaN     takes 00:00:01.1023320
 Lobe None None :       7.273         0.000      takes 00:00:00.0049967
 Lobe None Core :          -¡Þ           NaN     takes 00:00:00.0050010
 Lobe None Lobe :       2.522         0.068      takes 00:00:01.0267926
 Lobe Core None :      18.817         0.000      takes 00:00:00.0060042
 Lobe Core Lobe :      14.067         0.068      takes 00:00:00.9272274
 Lobe Lobe None :     -16.215         0.432      takes 00:00:01.0127659
 Lobe Lobe Core :          -¡Þ           NaN     takes 00:00:00.9742296
Time Performance:00:00:06.3263128
	*/
}
void testSpatialXResultTable(){
	cout<< "--------test SpatialXResultTable" << endl; 
	SpatialXResultTable tbl("../data/cstar/cstar_res_cpp.csv", false);
	cout<<"rows = "<<tbl.data.size()<<endl;
	cout<<tbl.data[0].toString()<<endl;
}
void testRadioXResultTable(){
	cout<< "--------test RadioXResultTable" << endl; 
	RadioXResultTable tbl("../data/cstar/radiox.csv", false);
	cout<<"rows = "<<tbl.data.size()<<endl;
	cout<<tbl.data[0].toString()<<endl;
	cout<<tbl.data[1].toString()<<endl;
	cout<<tbl.data[2].toString()<<endl;
}
void testNextRadios(){
	cout<< "--------test NextRadios" << endl; 
	SpatialXResultTable tbl("../data/cstar/cstar_res_cpp.csv", false);
	vector<SpatialXResult> com;
	
	map<int, int> res;
	int count=0;
	
	while((com=tbl.getNextCombination()).size()>0){
		if(res.find(com.size())==res.end()){
			res[com.size()] = 1;
		}else{
			res[com.size()] = res[com.size()] + 1;
		}
		count+=com.size();
	}
	for(map<int, int>::iterator it=res.begin(); it!=res.end(); it++){
		cout<<"size="<<it->first<<", count="<<it->second<<endl;
	}
	cout<<"total rows="<<count<<endl;
}

void testRunRadioX(){
	cout<< "--------test RunRadioX" << endl; 
	LobePrior* pri = new LobePrior_Rayleigh(44);
	printf("prior variance = %.14lf\n", pri->variance());
	
	int sampleCount = 2000;
	printf("sampleCounts = %d\n", sampleCount);
	
	SpatialXResultTable tbl("../data/cstar/cstar_res_cpp.csv", false);
	cout<<"reading ../data/cstar/cstar_res_cpp.csv"<<endl;
	vector<SpatialXResult> com;
	
	FILE* pf = fopen("../data/cstar/radiox.csv", "w");
	cout<<"doing radiox"<<endl;
	char buf[255];
	while((com=tbl.getNextCombination()).size()>0){
		SpatialPoint* optical = &(com[0].p1);
		optical->ra_err = 0.2;
		optical->dec_err = 0.788035;
		vector<SpatialPoint*> radios;
		for(vector<SpatialXResult>::iterator it=com.begin(); it!=com.end(); it++){
			it->p2.ra_err = 0.2;
			it->p2.dec_err = 0.788035;
			radios.push_back(&(it->p2));
		}
		
		vector<Hypothesis> hyps = Hypothesis::generateHypList(optical, radios);
		for(std::vector<Hypothesis>::iterator it=hyps.begin();it!=hyps.end();it++){
			double likelihood, likelihood_err;
			//it->likelihoodAsym(pri, sampleCounts, 0.2, likelihood, likelihood_err);
			it->likelihoodSym(pri, sampleCount, likelihood, likelihood_err);
			sprintf(buf, ",%lf,%lf\n", likelihood, likelihood_err);
			string str = it->combiStr();
			str += buf;
			fwrite(str.c_str(), 1, str.size(), pf);
		}
		fflush(pf);
	}
	
	fclose(pf);
}

void testRunRadioX2(){
	cout<< "--------test RunRadioX2" << endl; 
	LobePrior* pri = new LobePrior_Rayleigh(44);
	printf("prior variance = %.14lf\n", pri->variance());
	
	int sampleCount = 2000;
	printf("sampleCounts = %d\n", sampleCount);
	
	double sigmas[3]={0.2, 0.788035, 0.5};
	
	SpatialXResultTable tbl("../data/cstar/cstar_res_cpp.csv", false);
	cout<<"reading ../data/cstar/cstar_res_cpp.csv"<<endl;
	vector<SpatialXResult> com;
	cout<<"doing radiox"<<endl;
	RadioXRunner::doRadioX(true, pri, sigmas, tbl, "../data/cstar/radiox.csv", sampleCount);
	cout<<"done"<<endl;
}

void testResultAnalysis(){
	cout<< "--------test ResultAnalysis" << endl; 
	cout<<"analysing "<<endl;
	RadioXRunner::analysisResult("../data/cstar/radiox.csv", "../data/cstar/radiox_pick.csv");
	cout<<"done!"<<endl;
}

void testHypothesis() {
	//testSpatialXResultTable();
	//testRadioXResultTable();
	//testHypsCombination();
	//testMatrix();
	//testPGSigma();
	//testGuassian();
	//testRandom();
	testLikelihood();
	//testNextRadios();
	////testRunRadioX();
	////testRunRadioX2();
	//testResultAnalysis();
}

#endif