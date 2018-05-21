programName = ellip3d
CXX = CC -DOPENMP -fopenmp
OPT = -dynamic -std=c++0x -g -O3 -mtune=native -msse4.2 #-ffast-math -funroll-loops
INC = -I../include -I./include
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

%.o : %.cpp
	$(CXX) -c -o $@ $^ $(OPT) $(INC)

$(programName) : $(OBJECTS)
	$(CXX) -o $@ $^ $(OPT) $(INC)

clean :
	rm -rf *.o *.mod $(programName) *~ *.out

