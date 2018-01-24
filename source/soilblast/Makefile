programName = hu_code
CXX = g++
#CXX = CC -DOPENMP -fopenmp
OPT = -dynamic -std=c++0x -g -O3 -mtune=native -msse4.2 
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

