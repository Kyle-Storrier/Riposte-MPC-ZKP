CXX=g++
CXXFLAGS=--std=c++17 -march=native -O3 -g -Wno-ignored-attributes
LIBS=-lbsd

default: dpflowmc dpfaes mpc zkplowmc

mpc: MPC.cpp dpf.h prg.h LowMC.h LowMC.cpp block.h simulator.h verifier.h 
	$(CXX) $(CXXFLAGS) -o mpc LowMC.cpp MPC.cpp -DLOWMC $(LIBS)

dpflowmc: dpf.cpp dpf.h prg.h LowMC.h LowMC.cpp block.h
	$(CXX) $(CXXFLAGS) -o dpflowmc LowMC.cpp dpf.cpp -DLOWMC $(LIBS)

dpfaes: dpf.cpp dpf.h prg.h aes.h block.h 
	$(CXX) $(CXXFLAGS) -o dpfaes dpf.cpp -DAES $(LIBS)

zkplowmc: mpcZkp.cpp dpf.h prg.h LowMC.h  aes.h block.h 
	$(CXX) $(CXXFLAGS) -DLOWMC dpf.h LowMC.h prg.h aes.h block.h LowMC.cpp -o mpczkplowmc mpcZkp.cpp $(LIBS)

clean:
	rm -f dpflowmc dpfaes

