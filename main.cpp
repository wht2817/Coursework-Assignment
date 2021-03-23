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
	
	//Initialize MPI
	
	int err = MPI_Init(&argc, &argv);
	
	//End code if MPI failed to initialize
	
	if (err != MPI_SUCCESS){
		
		cout << "MPI failed to initialize." << endl;
		
		return -1;
	}
	
	int rank, size;
	
	//Get size and ranks
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    
	
	//Initialize Variables for input into class
	
	string InitialCondition;
	
	double dt, T, h;
	
	bool help = false;
	
	try{
		
		//Start recording of time
		
		auto begin = high_resolution_clock::now();
		
		//Parse Command Line Inputs
		
		po::variables_map vm;
	
		Options(argc, argv, vm, InitialCondition, dt, T, h, rank, help);
		
		if (help == true){
			
			MPI_Finalize();
			
			return 0;
		}
	
		//Initalize SPH object
		
		SPH SPHsolver(InitialCondition, dt, T, h, MPI_COMM_WORLD, rank, size);
	
		/* Solve SPH problem with solver() member function.
		 * Solver() utilizes other member functions and private variables within
		 * the SPH class to solve the SPH problem
		 */
		 
		SPHsolver.solver();
		
		//Output time taken for code to run
		if (rank == 0){

			auto stop = high_resolution_clock::now();
	
			auto duration = duration_cast<seconds>(stop - begin);

			cout << duration.count() << "seconds" << endl;
		}
		
		//Finalize MPI
		MPI_Finalize();
	}
	
	
	//Catch bad allocation error where no. of proccesses > no. of particles
	catch (const std::bad_alloc& e) {
		
		if (rank == 0){
		
			cout << "An error occured: " << e.what() << endl;
		
		}
		
		
		MPI_Finalize();
		
	}
	
	/*Catch logic errors for invalid inputs such as too many initial conditions
	 *or negative T, dt and h.
	 */
	 
	catch (const std::logic_error& e){
		
		if (rank == 0){
		
			cout << "An error occured: " << e.what() << endl;
		
		}
		
		
		MPI_Finalize();
	}
	

	
	return 0;
}
