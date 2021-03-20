#include <iostream>
#include <cmath>
#include "cblas.h"
#include <ctime>
#include <random>
#include "CommandLineOptions.h"
#include <fstream>
#include "mpi.h"
#include "SPH.h"
#include <iomanip>
#include <chrono>
using namespace std;
using namespace std::chrono;


int main(int argc, char *argv[])
{
	
	
	MPI_Init(&argc, &argv);
	int rank, size;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    
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
	
	
	auto begin = high_resolution_clock::now();
	
	try{
		
		SPH test("ic-droplet", 0.0001, 0.4, 0.01, MPI_COMM_WORLD, rank, size);
		
		test.solver();
	
		
	
		if (rank == 0){
			auto stop = high_resolution_clock::now();
	
			auto duration = duration_cast<seconds>(stop - begin);
	
			cout << duration.count() << "seconds" << endl;
		}
		
		MPI_Finalize();
	}
		
	catch (const std::bad_alloc& e) {
		
		cout << "An error occured: " << e.what() << endl;
		
		MPI_Finalize();
	}
	
//	
//	if (test.checksize() == false){
//		
//		MPI_Finalize();
//		return 0;
//		
//	}
//	
//
	
	return 0;
}
