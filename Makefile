CPP = g++
CC = gcc 
CFlag = -Wall -std=c++11
main.o: main.cpp
	$(CPP) $(CFlag) $^ -o $@
clean:
	rm main.o 
