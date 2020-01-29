CXX=g++
CXXFLAGS=--std=c++17 -march=native -O3 -g -Wno-ignored-attributes -mavx2
LIBS=-lbsd 

default: matrix

matrix: matrix.cpp
	$(CXX) $(CXXFLAGS) matrix.cpp -o matrix.out

clean:
	rm -f  matrix.out

