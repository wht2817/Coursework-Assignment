#include <iostream>
#include <cmath>
#include "cblas.h"
#include <ctime>
#include <random>
#include "CommandLineOptions.h"
#include <fstream>
#include "mpi.h"
#include "SPH.h"
using namespace std;



int main(int argc, char *argv[])
{
    
	//START OF  BOOST STUFF
	//Initialize Variables
	//string InitialCondition;
	
	//double dt, T, h;
	
	//po::variables_map vm;
	
	//check for valid inputs
	//bool checkvalid;
	
	//checkvalid = Options(argc, argv, vm, InitialCondition,dt,T,h);
	
	//if (!checkvalid){
		//End programme if there are invalid inputs
		//return 0;
	//}
	
	//SPH test(InitialCondition, dt, T, h, MPI_COMM_WORLD, rank, size);
	
	//END OF BOOST STUFF
	
	MPI_Init(&argc, &argv);
	int rank, size;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	
	SPH test("ic-two-particles", 0.0001, 0.0002, 0.01, MPI_COMM_WORLD, rank, size);
	
	if (!test.checksize()){
		MPI_Finalize();
		return 0;
		
	}
	//test.printV();
	test.solver();
	
	MPI_Finalize();
	return 0;
}
