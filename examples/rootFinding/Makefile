CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2
LDFLAGS = 
INCLUDE = ../../include/
LIB =
SRC = ../../src/
#OBJ = $(SRC:%.cpp=%.o)
OBJ = rootFinding.o main.o
OUT = main.exe

#%.o: $(INCLUDE) %.hpp $(SRC) %.cpp
#	 $(CXX) $(CXXFLAGS) $(SRC) -c %.cpp 
rootFinding.o : $(INCLUDE)rootFinding.hpp $(SRC)rootFinding.cpp
				$(CXX) $(CXXFLAGS) -c $(SRC)rootFinding.cpp

main.o : $(INCLUDE)rootFinding.hpp main.cpp
		 $(CXX) $(CXXFLAGS) -c main.cpp

main : $(OBJ)
	   $(CXX) $(CXXFLAGS) $(OBJ) -o $(OUT)

clean : 
	del rootFinding.o
	del main.o 
	del main.exe