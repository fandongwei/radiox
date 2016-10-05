CXX = g++

OBJS = Constants.o SpatialPoint.o SpatialPointTable.o SpatialXResultTable.o SpatialXResult.o RadioXResult.o RadioXResultTable.o SpatialXMatcher.o Gaussian.o Det.o Hypothesis.o LobePrior_Uniform.o LobePrior_Rayleigh.o LobePrior_LogNormal.o Ran1.o GaussianRandom2D.o RadioXRunner.o CommandParser.o

all: libradiox.a

libradiox.a : $(OBJS)
	ar rv libradiox.a $(OBJS)
	
Constants.o : Constants.cpp Constants.h
	$(CXX) -c Constants.cpp

SpatialPoint.o : SpatialPoint.cpp SpatialPoint.h Constants.h Gaussian.h
	$(CXX) -c SpatialPoint.cpp Gaussian.h -I/usr/local/include/eigen3

SpatialPointTable.o : SpatialPointTable.cpp SpatialPointTable.h SpatialPoint.h Constants.h
	$(CXX) -c SpatialPointTable.cpp

SpatialXResult.o : SpatialXResult.cpp SpatialXResult.h SpatialPoint.h
	$(CXX) -c SpatialXResult.cpp

SpatialXResultTable.o : SpatialXResultTable.cpp SpatialXResultTable.h SpatialXResult.h
	$(CXX) -c SpatialXResultTable.cpp

SpatialXMatcher.o : SpatialXMatcher.cpp SpatialXMatcher.h SpatialXResult.h SpatialPointTable.h SpatialPoint.h Constants.h
	$(CXX) -c SpatialXMatcher.cpp

Gaussian.o : Gaussian.cpp Gaussian.h Constants.h
	$(CXX) -c Gaussian.cpp -I/usr/local/include/eigen3
	
Det.o : Det.cpp Det.h Gaussian.h
	$(CXX) -c Det.cpp -I/usr/local/include/eigen3
	
LobePrior_Uniform.o : LobePrior_Uniform.cpp LobePrior_Uniform.h Constants.h LobePrior.h
	$(CXX) -c LobePrior_Uniform.cpp Constants.h

LobePrior_Rayleigh.o : LobePrior_Rayleigh.cpp LobePrior_Rayleigh.h Constants.h LobePrior.h
	$(CXX) -c LobePrior_Rayleigh.cpp Constants.h

LobePrior_LogNormal.o : LobePrior_LogNormal.cpp LobePrior_LogNormal.h Constants.h LobePrior.h
	$(CXX) -c LobePrior_LogNormal.cpp Constants.h

Ran1.o : Ran1.cpp Ran1.h
	$(CXX) -c Ran1.cpp
	
GaussianRandom2D.o : GaussianRandom2D.cpp GaussianRandom2D.h Ran1.h 
	$(CXX) -c GaussianRandom2D.cpp -I/usr/local/include/eigen3
	
Hypothesis.o : Hypothesis.cpp Hypothesis.h Gaussian.h SpatialPoint.h LobePrior.h Ran1.h GaussianRandom2D.h
	$(CXX) -c Hypothesis.cpp -I/usr/local/include/eigen3
	
RadioXResult.o : RadioXResult.cpp RadioXResult.h
	$(CXX) -c RadioXResult.cpp
	
RadioXResultTable.o : RadioXResultTable.cpp RadioXResultTable.h RadioXResult.h
	$(CXX) -c RadioXResultTable.cpp

RadioXRunner.o : RadioXRunner.cpp RadioXRunner.h RadioXResultTable.h SpatialXResultTable.h Hypothesis.h LobePrior.h RadioXResult.h
	$(CXX) -c RadioXRunner.cpp -I/usr/local/include/eigen3
	
CommandParser.o : CommandParser.cpp CommandParser.h
	$(CXX) -c CommandParser.cpp
	
.PHONY : clean
clean:
	-rm -f *.a *.o
