default:myprog

flags = -std=c++11 -Wall -O2

main.o: main.cpp
	g++ $(flags) -o main.o -c main.cpp

CommandLineOptions.o: CommandLineOptions.cpp CommandLineOptions.h
	g++ $(flags) -o CommandLineOptions.o -c CommandLineOptions.cpp
	
myprog: main.o CommandLineOptions.o
	g++ -o myprog main.o CommandLineOptions.o -L. -lblas -lboost_program_options