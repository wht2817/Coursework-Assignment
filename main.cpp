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
	string InitialCondition;
	
	double dt, T, h;
	
	try{
		
		auto begin = high_resolution_clock::now();
	
		po::variables_map vm;
	
		Options(argc, argv, vm, InitialCondition,dt,T,h);
	
		SPH SPHsolver(InitialCondition, dt, T, h, MPI_COMM_WORLD, rank, size);
	
		SPHsolver.solver();
		
		if (rank == 0){

			auto stop = high_resolution_clock::now();
	
			auto duration = duration_cast<seconds>(stop - begin);

			cout << duration.count() << "seconds" << endl;
		}
		
		MPI_Finalize();
	}
	
	catch (const std::bad_alloc& e) {
		
		if (rank == 0){
		
			cout << "An error occured: " << e.what() << endl;
		
		}
		MPI_Finalize();
		
	}
	
	catch (const std::logic_error& e){
		
		if (rank == 0){
		
			cout << "An error occured: " << e.what() << endl;
		
		}
		MPI_Finalize();
	}
	
	
//	auto begin = high_resolution_clock::now();
//	
//	try{
//		
//		SPH test("ic-dam-break", 0.0001, 1, 0.01, MPI_COMM_WORLD, rank, size);
//		
//		test.solver();
//	
//		
//	
//		if (rank == 0){
//			auto stop = high_resolution_clock::now();
//	
//			auto duration = duration_cast<seconds>(stop - begin);
//	
//			cout << duration.count() << "seconds" << endl;
//		}
//		
//		MPI_Finalize();
//	}
//		
//	catch (const std::bad_alloc& e) {
//		
//		cout << "An error occured: " << e.what() << endl;
//		
//		MPI_Finalize();
//	}
	

	
	return 0;
}
