programName = hu_code
CXX =  mpicxx -DOPENMP -fopenmp#icpc g++
OPT = -std=c++0x -g -O3 -mtune=native -msse4.2#-dynamic -ffast-math -funroll-loops
INC = -I../include -I./include
ARMA = -L../lib/ -larmadillo
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

%.o : %.cpp
	$(CXX) -c -o $@ $^ $(OPT) $(INC) $(ARMA)

$(programName) : $(OBJECTS)
	$(CXX) -o $@ $^ $(OPT) $(INC) $(ARMA)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out

