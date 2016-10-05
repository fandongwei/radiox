CXX = g++

###### RADIOX
RADIOXPATH = ../cpp
RADIOXLIB = -L$(RADIOXPATH) $(RADIOXPATH)/libradiox.a

OTHERLIBS = -L. -lm -lstdc++ 

OBJS = tester.o


all: tester

tester : $(OBJS)
	$(CXX) -o tester $(OBJS) $(RADIOXLIB) $(OTHERLIBS)

tester.o : tester.cpp testHypothesis.h testSpatialX.h $(RADIOXPATH)/libradiox.a
	$(CXX) -c tester.cpp -I$(RADIOXPATH)  -I/usr/local/include/eigen3


.PHONY : clean test run
clean:
	-rm -f tester *.o

test:
	-./tester a="ac;asda sdf adf" \
	 b=1 c=1.2 \
	 d=1,2,3 e=1.1,2.2,3.3 \
	 f=false g=true h=1 i=0

run:
	#-./tester opticaltbl='../data/cstar/cstar.csv' radiotbl='../data/cstar/cstar.csv' \
	#skipspatialx=true spatialxradius=0.03333333333333333333333333333333 spatialxresult='../data/cstar/spatialx.csv' \
	#skipradiox=false symmodel=false samplecounts=80,50,180 sigma=0.2,0.788035,0.5 lobepriortype=rayleigh lobepriorparam=44 radioxreslut='../data/cstar/radiox.csv' \
	#analysisresult='../data/cstar/analysisresult.csv'
	
	-./tester opticaltbl='../data/atlas/swire_ori.csv' radiotbl='../data/atlas/atlas_ori.csv' \
	skipspatialx=false spatialxradius=0.03333333333333333333333333333333 spatialxresult='../data/atlas/spatialx.csv' \
	skipradiox=false symmodel=false sigma=0.2,0.788035,0.5 lobepriortype=rayleigh lobepriorparam=44 radioxreslut='../data/atlas/radiox.csv' \
	analysisresult='../data/atlas/analysisresult.csv' samplecounts=40,30,80