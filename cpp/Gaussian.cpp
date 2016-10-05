/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "Gaussian.h"
#include "Constants.h"
#include <stdio.h>
#include <iostream>
Gaussian::Gaussian(){
}
void Gaussian::init(Eigen::VectorXd m, Eigen::MatrixXd V, bool invert){
	//puts("initted Gaussian");
	this->D = m.rows();
	if(this->D != V.rows())
		throw("Gaussian..ctor(): dim mismatch");
	if(V.cols() == 1){
		Eigen::MatrixXd nV(this->D, this->D);
		nV.setZero();
		for(int i=0;i<this->D;i++){
			nV(i, i) = V(i, 0);
		}
		V = nV;
	}
	if(V.rows() != V.cols()) throw("Gaussian..ctor(): matrix not valid");
	this->M = m;
	if(invert){
		this->F = V.inverse();
	}else{
		this->F = V;
	}
	this->S = 1;
	this->N = sqrt(pow(2*Constants::PI, this->D)/this->F.determinant());
}
Gaussian::Gaussian(Eigen::VectorXd m, Eigen::MatrixXd V, bool invert){
	init(m, V, invert);
}
Gaussian::Gaussian(double x, double y, double sigma, bool invert){
	Eigen::VectorXd m(2);
	m<<x,y;
	Eigen::MatrixXd V(2,2);
	V.setZero();
	V(0,0)=V(1,1)=sigma*sigma;
	init(m, V, invert);
}
double Gaussian::sandwich(Eigen::MatrixXd a, Eigen::MatrixXd M, Eigen::MatrixXd b){
	Eigen::MatrixXd X = M*b;
	Eigen::MatrixXd q = a.transpose() * X;
	if(q.rows()!=1 || q.cols() !=1) throw("Gaussian.Sandwich(): not scalar?");
	return q(0, 0);
}
double Gaussian::evaluate(Eigen::MatrixXd x){
	Eigen::MatrixXd t = M - x;
	double q = sandwich(t, F, t);
	//printf("q=%lf\n",q);
	return S / N * exp(-0.5 * q);
}
Gaussian Gaussian::multiply(Gaussian g, bool normalize){
	Gaussian r;
	r.D = this->D;
	r.F = this->F + g.F;
	
	if(r.F.rows()==r.F.cols()){
		//LuDecomposition
		//r.M = r.F.colPivHouseholderQr().solve(this->F * this->M + g.F * g.M);//colPivHouseholderQr
		//r.M = r.F.fullPivLu().solve(this->F * this->M + g.F * g.M);//fullPivLu
		r.M = r.F.lu().solve(this->F * this->M + g.F * g.M);//PartialPivLU
		//r.M = r.F.colPivHouseholderQr().solve(this->F * this->M + g.F * g.M);//colPivHouseholderQr
		//r.M = r.F.llt().solve(this->F * this->M + g.F * g.M);//LLT
	}else{
		//QrDecomposition
		//puts("\nQrDecomposition");
		r.M = r.F.colPivHouseholderQr().solve(this->F * this->M + g.F * g.M);//colPivHouseholderQr
	}
	
	r.N = sqrt(pow(2 * Constants::PI, r.D)/r.F.determinant());
	if(normalize) r.S = 1;
	else{
		double A = sandwich(r.M, r.F, r.M);
		//printf("A=%.14lf, ", A);
		A -= sandwich(this->M, this->F, this->M);
		A -= sandwich(g.M, g.F, g.M);
		r.S = r.N / this->N / g.N * exp(A / 2.0);
	}
	return r;
}
double Gaussian::likelihood(Eigen::MatrixXd model){
	return evaluate(model);
}
double Gaussian::sigma(){
	return sqrt(sigmaSquare());
}
double Gaussian::sigmaSquare(){
	Eigen::MatrixXd fin = F.inverse();
	return fin(0, 0);
}
double Gaussian::x(){
	return M(0,0);
}
double Gaussian::y(){
	return M(1,0);
}