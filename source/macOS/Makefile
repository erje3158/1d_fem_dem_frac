programName = hu_code
CXX = g++ 
OPT = -dynamic -std=c++0x -g -O3 -mtune=native -msse4.2 -ggdb3
ARMA = -larmadillo
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

%.o : %.cpp
	$(CXX) -c -o $@ $^ $(OPT) $(ARMA)

$(programName) : $(OBJECTS)
	$(CXX) -o $@ $^ $(OPT) $(ARMA)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out

