/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "Hypothesis.h"
#include "Constants.h"
#include "SpatialPoint.h"
#include "Gaussian.h"
#include "Det.h"
#include "Ran1.h"
#include "GaussianRandom2D.h"
#include "LobePrior.h"
#include <math.h>
#include <stdio.h>
#include <iostream>

Hypothesis::Hypothesis()
{
	optical = NULL;
	rcore = NULL;
	rlobes1 = NULL;
	rlobes2 = NULL;
}


Hypothesis::~Hypothesis()
{
}

std::vector<Hypothesis> Hypothesis::generateHypList(SpatialPoint* optical, std::vector<SpatialPoint*>& radios) {
	std::vector<Hypothesis> hyps;
	int n = radios.size();

	Hypothesis hyp;
	hyp.optical = optical;
	for (int k = -1; k < n; k++) {
		hyp.rcore = -1 == k ? NULL : radios[k];
		for (int i = -1; i < n; i++) {
			if (i == k && -1 != i)continue;
			if (-1 == i) {
				hyp.rlobes1 = hyp.rlobes2 = NULL;
				hyps.push_back(hyp);
				continue;
			}
			hyp.rlobes1 = radios[i];
			hyp.rlobes2 = NULL;
			hyps.push_back(hyp);
			for (int j = i + 1; j < n; j++) {
				if (j == k) continue;
				hyp.rlobes2 = radios[j];
				hyps.push_back(hyp);
			}
		}
	}

	return hyps;
}

void Hypothesis::print(){
	printf("%s %s %s ", rcore==NULL?"None":"Core", rlobes1==NULL?"None":"Lobe", rlobes2==NULL?"None":"Lobe");
}

double Hypothesis::get1stPrior(){
	int ex=0;
	if(NULL != rcore) ex++;
	if(NULL != rlobes1) ex++;
	if(NULL != rlobes2) ex++;
	return pow(Constants::WholeSphereInSquareArcsec, ex);
}
void Hypothesis::likelihoodSym(LobePrior* pri, int sampleCount, double& likelihood, double& likelihood_err) {
	//using namespace std;
	Gaussian goptical = optical->getGaussian(*optical);

	int lobescount = 0;
	if (NULL != rlobes1) lobescount++;
	if (NULL != rlobes2) lobescount++;

	if (0 == lobescount) {
		double int1 = 1;
		if (NULL != rcore) {
			Gaussian grcore = rcore->getGaussian(*optical);
			Gaussian gm = goptical.multiply(grcore, false);
			int1 = get1stPrior() * gm.S;
		}
		this->likelihood = likelihood = log10(int1);
		this->likelihood_err = likelihood_err = 0;
		return;
	}
	Gaussian first = goptical;
	if (NULL != rcore) {
		Gaussian grcore = rcore->getGaussian(*optical);
		first = goptical.multiply(grcore, false);
	}
	Eigen::MatrixXd sigma0 = first.F.inverse(), sigma1;

	Gaussian second;
	if (1 == lobescount) {
		second = rlobes1->getGaussian(*optical);
		sigma1 = second.F.inverse();
	}
	else if (2 == lobescount) {
	}

	Ran1 random;
	GaussianRandom2D rand0(&first.M, &sigma0, &random);

	double sum0 = 0;
	double sum0s = 0;
	for (int i = 0; i < sampleCount; i++) {
		Eigen::VectorXd m0 = rand0.nextSample();
		double p = 0;
		if (1 == lobescount) {
			GaussianRandom2D rand1(&second.M, &sigma1, &random);
			Eigen::VectorXd m1 = rand1.nextSample();
			double r2 = (m0 - m1).squaredNorm();
			p = pri->prior(r2);
		}
		else if (2 == lobescount) {
			Eigen::VectorXd m02 = m0 * 2;
			Gaussian grl1 = rlobes1->getGaussian(*optical);
			Gaussian grl2 = rlobes2->getGaussian(*optical);
			second = Det(m02(0, 0) - grl2.x(), m02(1, 0) - grl2.y(), grl2.sigma());
			second = second.multiply(grl1, false);
			//printf("second.S=%.14e,",second.S);
			sigma1 = second.F.inverse();
			GaussianRandom2D rand1(&second.M, &sigma1, &random);
			Eigen::VectorXd m1 = rand1.nextSample();
			double r2 = (m0 - m1).squaredNorm();
			p = pri->prior(r2);
		}
		sum0 += p;
		sum0s += p * p;
	}
	double firstprior = get1stPrior();
	sum0 = sum0 * first.S * firstprior / sampleCount;
	sum0s = firstprior * firstprior * first.S * first.S * sum0s / sampleCount;
	sum0s = sqrt(abs(sum0s - sum0 * sum0) / sampleCount);
	this->likelihood_err = likelihood_err = sum0s / sum0 / log(10);
	this->likelihood = likelihood = log10(sum0);
}
void Hypothesis::likelihoodAsym(LobePrior* pri, int sampleCount, double ksigma, double& likelihood, double& likelihood_err) {
	using namespace std;
	Gaussian goptical = optical->getGaussian(*optical);

	int lobescount = 0;
	if (NULL != rlobes1) lobescount++;
	if (NULL != rlobes2) lobescount++;

	if (0 == lobescount) {
		double int1 = 1;
		if (NULL != rcore) {
			Gaussian grcore = rcore->getGaussian(*optical);
			Gaussian gm = goptical.multiply(grcore, false);
			int1 = get1stPrior() * gm.S;
		}
		this->likelihood = likelihood = log10(int1);
		this->likelihood_err = likelihood_err = 0;
		return;
	}
	Gaussian first = goptical;
	if (NULL != rcore) {
		Gaussian grcore = rcore->getGaussian(*optical);
		first = goptical.multiply(grcore, false);
	}
	Eigen::MatrixXd sigma0 = first.F.inverse(), sigma1;

	Gaussian second;
	if (1 == lobescount) {
		second = rlobes1->getGaussian(*optical);
		sigma1 = second.F.inverse();
	}
	else if (2 == lobescount) {
	}

	Ran1 random;
	GaussianRandom2D rand0(&first.M, &sigma0, &random);

	double sum0 = 0;
	double sum0s = 0;
	for (int i = 0; i < sampleCount; i++) {
		double m0x, m0y;
		rand0.nextSample(m0x, m0y);

		double p = 0;
		if (1 == lobescount) {
			GaussianRandom2D rand1(&second.M, &sigma1, &random);
			{
				double m1x, m1y;
				rand1.nextSample(m1x, m1y);
				double r2 = (m0x - m1x)*(m0x - m1x) + (m0y - m1y)*(m0y - m1y);
				p = pri->prior(r2);
			}
		}
		else if (2 == lobescount) {
			Gaussian grl1 = rlobes1->getGaussian(*optical);
			Gaussian grl2 = rlobes2->getGaussian(*optical);
			double grl2x = grl2.x();
			double grl2y = grl2.y();
			double grl2sigma = grl2.sigma();

			double k = random.nextGasdev(0, ksigma);

			double sum2 = 0;
			double int2 = 0;
			if (k < -1)continue;//tophat for the k

			double m0_2kx = m0x * (2 + k);
			double m0_2ky = m0y * (2 + k);

			double k1 = k + 1;
			second = Det((m0_2kx - grl2x) / k1, (m0_2ky - grl2y) / k1, grl2sigma / k1).multiply(grl1, false);
			sigma1 = second.F.inverse();
			GaussianRandom2D rand1(&second.M, &sigma1, &random);
			double m1x, m1y;
			rand1.nextSample(m1x, m1y);
			double r2 = (m0x - m1x)*(m0x - m1x) + (m0y - m1y)*(m0y - m1y);
			p = pri->prior(r2);
			p *= second.S;//S is the integral of two Gaussian PDF multiplication, it is not 1
		}
		sum0 += p;
		sum0s += p * p;
	}
	double firstprior = get1stPrior();
	sum0 = sum0 * first.S * firstprior / sampleCount;
	sum0s = firstprior * firstprior * first.S * first.S * sum0s / sampleCount;
	sum0s = sqrt(abs(sum0s - sum0 * sum0) / sampleCount);
	this->likelihood_err = likelihood_err = sum0s / sum0 / log(10);
	this->likelihood = likelihood = log10(sum0);
}
std::string Hypothesis::combiStr(){
	std::string rid;
	char buf[255];
	sprintf(buf, "%ld,%d,", optical->id, getType());
	std::string str(buf);
	if(rcore!=NULL){
		str += "core";
		sprintf(buf, "%ld", rcore->id); rid += buf;
	}
	if(rlobes1!=NULL){
		if(rid.size()>0){
			str += "|";
			rid += "|";
		}
		str += "lobe";
		sprintf(buf, "%ld", rlobes1->id); rid += buf;
	}
	if(rlobes2!=NULL){
		if(rid.size()>0){
			str += "|";
			rid += "|";
		}
		str += "lobe";
		sprintf(buf, "%ld", rlobes2->id); rid += buf;
	}
	return str + "," + rid;
}
int Hypothesis::getType(){
	int type=0;
	if(rcore != NULL) type |= 1;
	if(rlobes1 != NULL) type |= 2;
	if(rlobes2 != NULL) type |= 4;
	return type;
}

std::string Hypothesis::toString(){
	std::string str(optical->toString());
	if(rcore != NULL){
		str += ",";
		str += rcore->toString();
	}
	if(rlobes1 != NULL){
		str += ",";
		str += rlobes1->toString();
	}
	if(rlobes2 != NULL){
		str += ",";
		str += rlobes2->toString();
	}
	return str;
}

bool Hypothesis::operator <(const Hypothesis& p){
	//return this->result(0,0) < p.result(0,0);
	return this->likelihood < p.likelihood;
}