CPP = g++
CC = gcc 
CFlag = -Wall
main.o: main.cpp
	$(CPP) $(CFlag) $^ -o $@
  
