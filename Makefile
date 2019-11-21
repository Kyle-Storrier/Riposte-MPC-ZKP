CXX=g++
CXXFLAGS=--std=c++17 -march=native -O3 -g -fopenmp -Wignored-attributes -Wno-ignored-attributes  #-DBOOST_ERROR_CODE_HEADER_ONLY #-DDEBUG #-DVERBOSE
LIBS=-lbsd

default: test dpflowmc dpfaes

test: LowMC.cpp LowMC.h block.h key.h test.cpp
	$(CXX) $(CXXFLAGS) LowMC.h block.h key.h LowMC.cpp -o test test.cpp $(LIBS)

dpflowmc: dpf.cpp dpf.h prg.h LowMC.h  aes.h block.h 
	$(CXX) $(CXXFLAGS) -DLOWMC dpf.h LowMC.h prg.h aes.h block.h key.h LowMC.cpp -o dpflowmc dpf.cpp $(LIBS)

dpfaes: dpf.cpp dpf.h prg.h LowMC.h  aes.h block.h 
	$(CXX) $(CXXFLAGS) -DAES dpf.h LowMC.h prg.h aes.h block.h key.h LowMC.cpp -o dpfaes dpf.cpp $(LIBS)

zkplowmc: mpcZkp.cpp dpf.h prg.h LowMC.h  aes.h block.h 
	$(CXX) $(CXXFLAGS) -DLOWMC dpf.h LowMC.h prg.h aes.h block.h key.h LowMC.cpp -o mpczkplowmc mpcZkp.cpp $(LIBS)

clean:
	rm -f test dpflowmc dpfaes mpczkplowmc

