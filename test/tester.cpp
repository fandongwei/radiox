/*
@author Dongwei Fan <fandongwei@nao.cas.cn> http://www.lamost.org/~dwfan/
sponsored by NSFC (National Natural Science Foundation of China) No. 11503051
Institute: National Astronomical Observatories, Chinese Academy of Sciences
and also from Chinese Virtual Observatory http://www.china-vo.org/
please cite http://mnras.oxfordjournals.org/content/451/2/1299
*/
#include "testSpatialX.h"
#include "testHypothesis.h"
#include "CommandParser.h"
#include <time.h>

void testCommandParser(int argc, char* argv[]);
void runCstar(int argc, char* argv[]);
void runSwireAtlas(int argc, char* argv[]);
void testMap();

int main(int argc, char* argv[]) {
	//testCommandParser(argc, argv);
	//testMap();
	//testSpatialX();
	//testHypothesis();
	//runCstar(argc, argv);
	runSwireAtlas(argc, argv);
	return 0; 
}

void testMap(){
	cout << "--------test Map" << endl;
	map<string, int> xxx;
	cout<<"xxx[\"a\"]="<<xxx["a"]<<endl;
}

void testCommandParser(int argc, char* argv[]){
	cout << "--------test CommandParser" << endl;
	CommandParser cmd(argc, argv);
	cout<<"parameters:"<<endl;
	cmd.printAll();
	cout<<endl<<endl<<"test type"<<endl;
	cout<<"a="<<cmd.get("a")<<endl;
	cout<<"a="<<cmd.get("a")<<endl;
	cout<<"b="<<cmd.getInt("b")<<endl;
	cout<<"c="<<cmd.getDouble("c")<<endl;
	int ai[3];
	cmd.getIntArr("d", ai);
	printf("d=%d,%d,%d\n", ai[0], ai[1], ai[2]);
	double ad[3];
	cmd.getDoubleArr("e", ad);
	printf("e=%lf,%lf,%lf\n", ad[0], ad[1], ad[2]);
	printf("f=%d, g=%d, h=%d, i=%d\n", cmd.getBoolean("f"), cmd.getBoolean("g"), cmd.getBoolean("h"), cmd.getBoolean("i"));
}

void runcstar(int argc, char* argv[]){
	puts("\n\nstart running ...");
	CommandParser cmd(argc, argv);
	cmd.printAll();
	//spatialx
	if(!cmd.getBoolean("skipspatialx")){
		//load optical & radio catalogs
		cout<<"reading optical:"<<cmd.get("opticaltbl")<<endl;
		SpatialPointTable opticaltbl(cmd.get("opticaltbl"), true);
		cout<<"reading radio:"<<cmd.get("radiotbl")<<endl;
		SpatialPointTable radiotbl(cmd.get("radiotbl"), true);
		//do spatialx
		cout<<"doing spatialx..."<<endl;
		vector<SpatialXResult> res = SpatialXMatcher::matchTables(opticaltbl, radiotbl, cmd.getDouble("spatialxradius"),  true);
		cout<<"write spatialx result to "<<cmd.get("spatialxresult")<<endl;
		//write spatialx result
		SpatialXMatcher::outputSpatialXResult(res, cmd.get("spatialxresult"));
	}else{
		cout<<"skip spatialx"<<endl;
	}
	//radiox
	if(!cmd.getBoolean("skipradiox")){
		cout<<"doing radiox"<<endl;
		//prepare params
		LobePrior* pri = NULL;
		double lobepriorparam[3]; cmd.getDoubleArr("lobepriorparam", lobepriorparam);
		if(cmd.get("lobepriortype").compare("rayleigh")==0) pri = new LobePrior_Rayleigh(lobepriorparam[0]);
		int samplecounts[3]; cmd.getIntArr("samplecounts", samplecounts);
		double sigmas[3]; cmd.getDoubleArr("sigma", sigmas);
		
		//do radiox
		RadioXRunner::doRadioX(true, pri, sigmas, cmd.get("spatialxresult"), cmd.get("radioxreslut"), samplecounts);
		delete pri;
	}else{
		cout<<"skip radiox"<<endl;
	}
	
	cout<<"analysing radiox data"<<endl;
	//do analysis
	RadioXRunner::analysisResult(cmd.get("radioxreslut"), cmd.get("analysisresult"));
	cout<<"analysing done!"<<endl;
	
	//parse analysis result
	RadioXResultTable tbl(cmd.get("analysisresult"), false);
	//get triplets,doublets,cores
	vector<RadioXResult> triplets = tbl.getTriplets();	
	cout<<"size of triplets="<<triplets.size()<<endl;
	
	vector<RadioXResult> doublets = tbl.getDoublets();
	cout<<"size of doublets="<<doublets.size()<<endl;
	
	vector<RadioXResult> cores = tbl.getCores();
	cout<<"size of cores="<<cores.size()<<endl;
}
vector<SpatialPoint> readSwire(string filename){
	vector<SpatialPoint> optical;
	FILE* pf = fopen(filename.c_str(), "r");
	char buf[1024];
	fscanf(pf, "%s\n", buf);//skip first line
	long id;
	double ra, dec;
	while(!feof(pf)){
		fscanf(pf, "%s\n", buf);
		stringstream ss;
		ss.str(buf);
		string item;
		getline(ss, item, ','); id = atoi(item.c_str());
		getline(ss, item, ',');
		getline(ss, item, ','); ra = atof(item.c_str());
		getline(ss, item, ','); dec = atof(item.c_str());
		//cout<<id<<","<<ra<<","<<dec<<endl;
		SpatialPoint sp(id, ra, dec, 0.2, 0.2);
		optical.push_back(sp);
	}
	fclose(pf);
	return optical;
}
vector<SpatialPoint> readAtlas(string filename){
	vector<SpatialPoint> radio;
	FILE* pf = fopen(filename.c_str(), "r");
	char buf[1024];
	fscanf(pf, "%s\n", buf);//skip first line
	long id=0;
	double ra, dec;
	while(!feof(pf)){
		id++;
		fscanf(pf, "%s\n", buf);
		stringstream ss;
		ss.str(buf);
		string item;
		getline(ss, item, ',');//skip ID
		getline(ss, item, ',');//skip name
		getline(ss, item, ','); ra = atof(item.c_str()); ra *= Constants::Radian2Degree;
		getline(ss, item, ',');//skip ra_degress
		getline(ss, item, ','); dec = atof(item.c_str()); dec *= Constants::Radian2Degree;
		//cout<<id<<","<<ra<<","<<dec<<endl;
		SpatialPoint sp(id, ra, dec, 0, 0);
		radio.push_back(sp);
	}
	fclose(pf);
	return radio;
}
void runSwireAtlas(int argc, char* argv[]){
	cout<<"\n\nparsing command  ......"<<endl;
	clock_t ttotal = clock(), t1 = clock();
	CommandParser cmd(argc, argv);
	cmd.printAll();
	clock_t t2 = clock(); printf("command parsed, takes %f second\n",float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	
	//spatialx
	if(!cmd.getBoolean("skipspatialx")){
		cout<<"\ndoing spatialx  ......"<<endl;
		//load swire
		string filename = cmd.get("opticaltbl");
		cout<<"reading optical: "<<filename;
		vector<SpatialPoint> opticaldata = readSwire(filename);
		SpatialPointTable opticaltbl(opticaldata, true);
		cout<<", imported "<<opticaltbl.data.size()<<" rows"<<endl;
		//load atlas
		filename=cmd.get("radiotbl");
		cout<<"reading radio: "<<filename;
		vector<SpatialPoint> radiodata = readAtlas(filename);
		SpatialPointTable radiotbl(radiodata, true);
		cout<<", imported "<<radiotbl.data.size()<<" rows"<<endl;
		//do spatialx
		vector<SpatialXResult> res = SpatialXMatcher::matchTables(opticaltbl, radiotbl, cmd.getDouble("spatialxradius"),  true);
		filename = cmd.get("spatialxresult");
		cout<<"write spatialx result to "<<filename<<", with "<<res.size()<<" rows"<<endl;
		//write spatialx result
		SpatialXMatcher::outputSpatialXResult(res, filename);
		t2 = clock(); printf("spatialx done! takes %f second\n",float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	}else{
		cout<<"\nskip spatialx\n\n";
	}
	//radiox
	if(!cmd.getBoolean("skipradiox")){
		cout<<"\ndoing radiox ......"<<endl;
		//prepare params
		LobePrior* pri = NULL;
		double lobepriorparam[3]; cmd.getDoubleArr("lobepriorparam", lobepriorparam);
		if(cmd.get("lobepriortype").compare("rayleigh")==0) pri = new LobePrior_Rayleigh(lobepriorparam[0]);
		int samplecounts[3]; cmd.getIntArr("samplecounts", samplecounts);
		double sigmas[3]; cmd.getDoubleArr("sigma", sigmas);
		bool sym = cmd.getBoolean("symmodel");
		string infile = cmd.get("spatialxresult");
		string outfile = cmd.get("radioxreslut");

		//do radiox
		RadioXRunner::doRadioX(sym, pri, sigmas, infile, outfile, samplecounts);
		delete pri;
		t2 = clock(); printf("radiox done! takes %f second\n",float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	}else{
		cout<<"\nskip radiox\n";
	}

	cout<<"\nanalysing radiox data"<<endl;
	//do analysis
	RadioXRunner::analysisResult(cmd.get("radioxreslut"), cmd.get("analysisresult"));
	t2 = clock(); printf("analysised result, takes %f second\n",float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;

	//parse analysis result
	RadioXResultTable tbl(cmd.get("analysisresult"), false);
	//pick triplets
	vector<RadioXResult> triplets = tbl.getTriplets();	
	t2 = clock(); printf("picked triplets, get %ld rows, takes %f second\n", triplets.size(), float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	//pick doublets
	vector<RadioXResult> doublets = tbl.getDoublets();
	t2 = clock(); printf("picked doublets, get %ld rows, takes %f second\n", doublets.size(), float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	//pick cores
	vector<RadioXResult> cores = tbl.getCores();
	t2 = clock(); printf("picked singles, get %ld rows, takes %f second\n", cores.size(), float(t2 - t1) / CLOCKS_PER_SEC);t1=t2;
	
	//everything done
	printf("\neverything done! takes %f second\n",float(clock() - ttotal) / CLOCKS_PER_SEC);
}