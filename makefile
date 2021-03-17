default:myprog

flags = -std=c++11 -Wall 

main.o: main.cpp
	g++ -g $(flags) -o main.o -c main.cpp

CommandLineOptions.o: CommandLineOptions.cpp CommandLineOptions.h
	g++ -g $(flags) -o CommandLineOptions.o -c CommandLineOptions.cpp
	
SPH.o: SPH.cpp SPH.h
	g++ -g $(flags) -o SPH.o -c SPH.cpp

myprog: main.o CommandLineOptions.o SPH.o
	g++ -g -o myprog main.o CommandLineOptions.o SPH.o -L. -lblas -lboost_program_options