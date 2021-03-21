default:myprog

flags = -std=c++11 -Wall 

main.o: main.cpp
	mpicxx -g $(flags) -o main.o -c main.cpp

CommandLineOptions.o: CommandLineOptions.cpp CommandLineOptions.h
	mpicxx -g $(flags) -o CommandLineOptions.o -c CommandLineOptions.cpp
	
SPH.o: SPH.cpp SPH.h
	mpicxx -g $(flags) -o SPH.o -c SPH.cpp

myprog: main.o CommandLineOptions.o SPH.o
	mpicxx -g -o myprog main.o CommandLineOptions.o SPH.o -L. -lblas -lboost_program_options