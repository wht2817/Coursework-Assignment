#include <iostream>
#include <cmath>
#include "cblas.h"
#include <ctime>
#include <random>
#include "CommandLineOptions.h"
#include <fstream>
#include "mpi.h"
#include <iomanip>

#include "SPH.h"

using namespace std;
/////
//Constructor
/////

SPH::SPH(string configuration, double pdt, double pT, double ph, MPI_Comm pcomm, int prank, int psize){
	
	
	
	//Assign values to private variables
	
	dt = pdt;							// Time step
	
	T  = pT;							// Total simulation time
	
	h = ph;								// Radius of influence
	
	hprime = 1.0 / h; 					// 1/h precalculated to be used to calculate q
	
	particles = configuration;			// Initial condition
				
	this -> comm = pcomm;				//Communicator variable (MPI_COMM_WORLD)
	
	this -> rank = prank;				//Rank of process
	
	this -> size = psize;				//No. of processes
	
	//Assign number of particles based on option chosen by user
	
	if (particles == "ic-one-particle"){
		N = 1;
	}
	else if (particles == "ic-two-particles"){
		N = 2;
	}
	else if (particles == "ic-three-particles"){
		N = 3;
	}
	else if (particles == "ic-four-particles"){
		N = 4;
	}
	else if (particles == "ic-dam-break"){
		N = 400;
	}
	else if (particles == "ic-block-drop"){
		N = 651;
	}
	else if (particles == "ic-droplet"){
		N = 331;
	}
	
	//Throw exception if number of processes exceeds number of particles. 
	
	if (size > N) {
		
		if (rank == 0){
			
			cout << "Number of proccesses cannot exceed " << N << " for "<< particles << " case." << endl;
			
		}
		
		throw std::bad_alloc();
		
	}
	
			
				
}

SPH::~SPH()
{
	//deallocate memory
	delete[] x[0];
	
	delete[] x;
	
	delete[] x_root[0];
	
	delete[] x_root;
	
	delete[] v[0];
	
	delete[] v;
	
	delete[] v_root[0];
	
	delete[] v_root;
	
	delete[] q;
	
	delete[] rho;
	
	delete[] rho_root;
	
	delete[] p;
	
	delete[] p_root;

	delete[] Fp[0];
	
	delete[] Fg[0];
	
	delete[] Fv[0];
	
	delete[] Fp;
	
	delete[] Fv;
	
	delete[] Fg;
	
	delete[] a[0];
	
	delete[] a;
	
	delete[] r;

	delete[] sendsizex;
	
	delete[] stridex;
	
	delete[] sendsize;
	
	delete[] stride;
	
}

			
/*
	Member Functions
*/
			
/////
//Iteration Functions
/////

/*Initialize parameters to aid in distributing work among processors
 * r        : Remainder of N divided by no. of processes
 * k        : Dummy variable to determine no. of particles distributed to process
 * lengthloc: No. of particles distributed to process
 * sendsize : Array of sendcounts where sendsize[i] is the send count of rank[i]
 * stride   : Array of displacement from start of send/receive buffer where stride[i] is the displacement of rank[i]
 * */
void SPH::initMPIVariables(){
	
	//Setting up parameters to distribute data evenly between proccesses
	
	int r = N % size;
	int k = (N - r) / size;
				
	if (rank < r) {
		
		k++;
		start = k * rank;
		finish   = k * (rank + 1);
	
	}
	
	else{
		
		start = (k + 1) * r + k * (rank -r);
		finish   = (k + 1) * r + k * (rank -r + 1);
	
	}

	lengthloc = finish - start;
	
	//Setting up arrays to use in MPI_Scatterv and MPI_Allgatherv
	
	sendsizex = new int[size](); 			// Array of sendcounts for x and v

	stridex   = new int[size](); 			// Array of displ for x and v


	sendsize  = new int[size](); 			// Array of sendcounts for p and rho
	
	stride    = new int[size](); 			// Array of displ for p and rho
	
	int k_scatter = 0;
	
	int scatter_start;
	
	int scatter_finish;
	
	/* Calculate send counts and displacement from start of array for each rank
	 * for MPI_Scatterv and MPI_Gatherv
	 */
	
	for (int i = 0; i < size; i++){
		
		k_scatter = (N - r) / size;
		
		if (i < r){
			
			k_scatter++;
			
			scatter_start = k_scatter * i;
			
			scatter_finish = k_scatter * (i + 1);
		
		}
		
		else{
			
			scatter_start = (k_scatter + 1) * r + k_scatter * (i - r);
		
			scatter_finish = (k_scatter + 1) * r + k_scatter * (i - r + 1);
			
		}
		
		sendsizex[i] = (scatter_finish - scatter_start) * 2;          //Multiplied by 2 to account for x and y coordinate 
		
		stridex[i]   = scatter_start * 2;							  //Multiplied by 2 to account for x and y coordinate
		
		sendsize[i]  = scatter_finish - scatter_start;
		
		stride[i]    = scatter_start;
		
		
	}	
	
	
}
/* Initialize root arrays and particle position. Data will be gathered to and scattered from these arrays.
 * Even though they are named root, they exist in all processes and correspond to x_j, v_j etc...
 * x_root    : Array of pointers pointing to x_rootpool.
 * x_rootpool: Pointer to contiguous block of memory containing particle positions.
 * v_root    : Array of pointers pointing to v_rootpool.
 * v_rootpool: Pointer to contiguous block of memory containing particle velocities.
 * rho_root  : Array which will contain all particle densities.
 * p_root    : Array which will contain all particle pressures.
 */ 

void SPH::initRootParticles(){
	
	//Allocate memory for root(shared) variables in all proccesses
	
	x_root     = new double*[N];
	
	x_rootpool = new double[N*2]();
	
	v_root     = new double*[N];
	
	v_rootpool = new double[N*2]();
	
	rho_root   = new double[N]();
	
	p_root     = new double[N]();
	
	//Forming contiguous 2D array for x_root and v_root
	
	for (int i = 0; i < N; ++i, x_rootpool += 2, v_rootpool += 2){
		x_root[i] = x_rootpool;
		v_root[i] = v_rootpool;
	}
	
	
	//Initialize Particle positions in x_root
	
	if (particles == "ic-one-particle"){
		
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.5;
	}
	
	else if (particles == "ic-two-particles"){ 
	
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.5;
		x_root[1][0] = 0.5;
		x_root[1][1] = h;
			
	}
	
	else if (particles == "ic-three-particles"){
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.5;
		x_root[1][0] = 0.495;
		x_root[1][1] = h;
		x_root[2][0] = 0.505;
		x_root[2][1] = h;
	}

	else if (particles == "ic-four-particles"){

		x_root[0][0] = 0.505;
		x_root[0][1] = 0.5;
		x_root[1][0] = 0.515;
		x_root[1][1] = 0.5;
		x_root[2][0] = 0.51;
		x_root[2][1] = 0.45;
		x_root[3][0] = 0.5;
		x_root[3][1] = 0.45;

	}

	else if (particles == "ic-dam-break"){
		
		/* Seeding with time was avoided so that results could be validated 
		 * in case changes needed to be made to code.
		 */
		
		for (int i = 0; i < 20; ++i){
			for (int j = 0; j < 20; ++j){
				
				x_root[i*20+j][0] = 0.01 +  i*0.01 + (double)rand()/(RAND_MAX*1000.0); //Psuedo-random values for noise 
				x_root[i*20+j][1] = 0.01 +  j*0.01 + (double)rand()/(RAND_MAX*1000.0); //Psuedo-random values for noise

			}
		}
					
	}
	
	else if (particles == "ic-block-drop"){

		for (int i = 0; i < 21; ++i){
			for (int j = 0; j < 31; ++j){

				x_root[i*31+j][0] = 0.1 + i * 0.01 + (double)rand()/(RAND_MAX*1000.0); //Psuedo-random values for noise
				x_root[i*31+j][1] = 0.3 + j * 0.01 + (double)rand()/(RAND_MAX*1000.0); //Psuedo-random values for noise
				
			}
		}
	}

	else if (particles == "ic-droplet"){
				
		//Initialize plotting parameters
		int Npoints[11] = {1, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60}; 					  //No. of points on each circumference to give evenly spaced points
		double R[11]    = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1}; // Radii at which the circle will be formed.
		
		//initialize dummy variable to aid in populating array with evenly spaced points
		int loopdummy = 0;
		double angledummy = 0.0;
					
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.7;
					
		//No noise added as particles are not directly above or beside each other
		
		for (int i = 1; i < 11; ++i){
			for (int j = Npoints[i-1]+loopdummy; j < Npoints[i] + Npoints[i-1] + loopdummy; j++){
				
				x_root[j][0] = 0.5 + R[i]*cos(2.0 * (j-angledummy)*M_PI/Npoints[i]);
				x_root[j][1] = 0.7 + R[i]*sin(2.0 * (j-angledummy)*M_PI/Npoints[i]);
			}

			loopdummy += Npoints[i-1];
			angledummy += Npoints[i];
				
		}
					
					
	}
	
	
}


/* Initialize local arrays and particle positions. Data of particles allocated to a process
 * will be stored in these arrays. These corresponding to x_i, v_i etc...
 * x     : Array of pointers pointing to xpool.
 * xpool : Pointer to contiguous block of memory containing particle positions.
 * v     : Array of pointers pointing to vpool.
 * vpool : Pointer to contiguous block of memory containing particle velocities.
 * rho   : Array which will contain all particle densities.
 * p     : Array which will contain all particle pressures.
 * Fp    : Array of pointers pointing to Fppool.
 * Fppool: Pointer to contiguous block of memory containing local pressure forces.
 * Fv    : Array of pointers pointing to Fvpool.
 * Fvpool: Pointer to contiguous block of memory containing local viscous forces.
 * Fg    : Array of pointers pointing to Fgpool.
 * Fgpool: Pointer to contiguous block of memory containing local gravitational forces.
 * a     : Array of pointers pointing to apool.
 * apool : Pointer to contiguous block of memory containing local acceleration.
 * q	 : Array containing all q values.
 */ 

void SPH::initLocParticles(){
	
	//Allocating Memory to local variable arrays
				
	x     = new double*[lengthloc];					//Position of particles
	
	xpool = new double[lengthloc*2]();				//Particle position pool for contiguous 2D position array
			
	v     = new double*[lengthloc];               	//Velocity of particle
	
	vpool = new double[lengthloc*2]();				//Velocity pool for contiguous 2D velocity array
			
	q     = new double[lengthloc*N]();	   			//q array
				
	rho   = new double[lengthloc]();                //Density
			
	p     = new double[lengthloc]();                //Pressure
			
	//Forces
	Fp    = new double*[lengthloc];          	    //Pressure Force
	
	Fppool = new double[lengthloc*2]();
			
	Fv    = new double*[lengthloc];         	    //Viscous Force
	
	Fvpool = new double[lengthloc*2]();

	Fg    = new double*[lengthloc];             	//Gravitational Force
	
	Fgpool = new double[lengthloc*2]();
			
	a     = new double*[lengthloc];              	//Acceleration
	
	apool = new double[lengthloc*2]();				//Acceleration pool for contiguous 2D acceleration array
	
				
	//Form contiguous 2D arrays so that they can be used in BLAS, MPI and exploit cache locality in loops.
	
	for (int i = 0; i < lengthloc; ++i, xpool += 2, vpool += 2, apool += 2, Fppool += 2, Fvpool += 2, Fgpool += 2){
		x[i]  = xpool;
		v[i]  = vpool;
		a[i]  = apool;
		Fp[i] = Fppool;
		Fv[i] = Fvpool;
		Fg[i] = Fgpool;
	}
				
}

			
/*Scatter sendsizex[rank] elements from x_root to x.
 */ 

void SPH::initX(){
	 
	MPI_Scatterv(&x_root[0][0], sendsizex, stridex, MPI_DOUBLE, &x[0][0], sendsizex[rank], MPI_DOUBLE, 0, comm);
	
}
			
			
/*Precalculate q values, storing them in the q array
 */ 

void SPH::calcQRVIJ(){
	
	
	
	for (int i = 0; i < lengthloc; ++i){
		for(int j = 0; j < N; ++j){
					
			//Calculate r first as it is a repeated value
			
			r[0] = x[i][0] - x_root[j][0];
			
			r[1] = x[i][1] - x_root[j][1];
						
			//Calculate qij once per timestep as it is a repeated value
			
			q[i*N+j] = sqrt(r[0] * r[0] + r[1] * r[1]);
						
						
		}
	}
	
	//Scale q matrix by hprime (1/h) using BLAS
	
	cblas_dscal(N * lengthloc, hprime, q, 1);
}



/* Calculate Density
 * q_coeff  : Precalculated repeated expression in rho calculation.
 * coeff_rho: Repeated expression calculated during construction and after mass scaling
 */ 

void SPH::calcRho(){
				
	double q_coeff;
				
				
	for (int i = 0; i < lengthloc; ++i){
		
		for (int j = 0; j < N; ++j){
						
			if (q[i*N+j] < 1){
				
				q_coeff = 1.0 - (q[i* N + j]*q[i* N + j]);
			
				rho[i] += coeff_rho * q_coeff * q_coeff * q_coeff;
							
				}
			}
		}
					
}


/* Calculate pressure
 */ 
	
void SPH::calcP(){
				
	for (int i = 0; i < lengthloc; ++i){

		p[i] = k*(rho[i] - rho_0);

	}
}
			
/* Scale Mass
 */ 

void SPH::scaleMass(){

	m = N * rho_0/cblas_dasum(N, rho_root, 1);

}
			
/* Calculate forces and acceleration
 * qcoeff : Precalculated repeated expression in both Fv and Fp
 * Fp_calc: Precalculated repeated expression in Fp
 * Fv_calc: Precalculated repeated expression in Fv
 */ 
			
void SPH::calcFA(){
				
	
	//Initialize repeated coefficients
	
	double qcoeff;
	double FP_calc;
	double FV_calc;
	
	
	for (int i = 0; i < lengthloc; ++i){
		for (int j = 0; j < N; ++j){
				
			if ((q[i*N+j] < 1) && ((i + start) != j)){
				
				//Precalculate repeated coefficients qcoeff, FP_calc and FV_calc
				qcoeff = (1.0-q[i* N + j]);
								
				FP_calc =  (coeff_p * (p[i]+p_root[j]) * qcoeff * qcoeff/(rho_root[j]*q[i* N + j]));
								
				FV_calc = (coeff_v * qcoeff / rho_root[j]);
				
				//Calculate Fp
				
				Fp[i][0] += FP_calc * (x[i][0] - x_root[j][0]);

				Fp[i][1] += FP_calc * (x[i][1] - x_root[j][1]);
						
				//Calculate Fv
								
				Fv[i][0] += FV_calc * (v[i][0] - v_root[j][0]);
				
				Fv[i][1] += FV_calc * (v[i][1] - v_root[j][1]);
								
			}
						
		}
		
		//Calculate Fg
						
		Fg[i][1] = -rho[i]*g;
						
		//Calculate Acceleration

		a[i][0] = (Fp[i][0] + Fv[i][0] + Fg[i][0])/rho[i];
		
		a[i][1] = (Fp[i][1] + Fv[i][1] + Fg[i][1])/rho[i];
		
	}
}
					
			
/* Enforce boundary conditions
 */ 

void SPH::calcBC(){

	for (int i = 0; i < lengthloc; ++i){
					
		//Right bc
		
		if (x[i][0] >= (1.0 - h)){
			
			x[i][0] = 1.0 - h;
			v[i][0] = -e * v[i][0];
		
		}

		//Left bc
		
		else if (x[i][0] <= h){
						
			x[i][0] = h;
			v[i][0] = -e * v[i][0];
		
		}
					
		//Top bc
		
		if (x[i][1] >= (1.0 - h)){
						
			x[i][1] = 1.0 - h;
			v[i][1] = -e * v[i][1];
		
		}

		//Bottom bc
		else if (x[i][1] <= h){
						
			x[i][1] = h;
			v[i][1] = -e * v[i][1];
		
		}
			
	}
		
}
			
/*Calculate Energy
*/

void SPH::calcE(){
				
	//Reset Energy values everytime they are calculated
	Ek = 0;
	Ep = 0;
	Et = 0;
	
	//Calculate Energy
	
	for (int i = 0; i < lengthloc; ++i){
					
		Ek += (v[i][0] * v[i][0] + v[i][1] * v[i][1]);
		
		Ep += x[i][1];
					
	}
				
	Ek *= 0.5*m;
	Ep *= m*g;
	Et  = Ek + Ep;
}
			
	
/////
//Main Solver
/////			
			
void SPH::solver(){
	
	//Initialize variables to aid MPI scattering and gathering
	
	initMPIVariables();
	
	//Initialize Particle Positions and allocate memory into root (shared) arrays
	
	initRootParticles();
	
	//Allocate Memory for local arrays
	
	initLocParticles();
	
	//Distribute particle positions to x arrays in individual processes
	
	initX();
				
	//Open energy file for output at every time interval
	
	ofstream EnergyOut("energy.txt", ios::out | ios::trunc);
				
		if(!EnergyOut.good()){
					
			cout << "Energy file could not be opened." << endl;
				
		}

				
	//Initialize time loop
	
	double t = 0;
	
	if (rank == 0){
		
		cout << "Evaluating..."<< endl; //Friendly message to let user know it is running
	
	}
	
	//Start time loop iteration
	
	while (t < T){
					
					
		//Calculate Q
			
		calcQRVIJ();
													
		//If First step, scale mass and recalculate density and force coefficients
		
		if (t == 0){
			
			//Calculate Density
			
			calcRho();

			//Gather individual density matrices in each process to rho_root array for use in force calculations

			MPI_Allgatherv(rho, sendsize[rank], MPI_DOUBLE, rho_root, sendsize, stride, MPI_DOUBLE, comm);
						
			//Scale Mass (inlined as it is a simple function)		
			
			m = N * rho_0/cblas_dasum(N, rho_root, 1);
						
			//Reset Rho and recalculate coefficient for rho calculation
						
			coeff_rho = m*4.0/(M_PI*h*h);
						
			cblas_dscal(lengthloc, 0.0, rho, 1);
						
			//Recalculate Rho and P and Allgather to root arrays
			
			calcRho();
			
			MPI_Allgatherv(rho, sendsize[rank], MPI_DOUBLE, rho_root, sendsize, stride, MPI_DOUBLE, comm);					
						
			calcP(); 
		
			MPI_Allgatherv(p, sendsize[rank], MPI_DOUBLE, p_root, sendsize, stride, MPI_DOUBLE, comm);
						
			//Rescale calculate Fp and Fv coefficients after mass has been calculated
						
			coeff_p   = (-m/2.0)*(-30.0/(M_PI*h*h*h));

			coeff_v   = -40.0*mu*m/(M_PI*h*h*h*h);
						     
						
		}

		else{
			
			//Calculate Density
			
			calcRho();
			
			MPI_Allgatherv(rho, sendsize[rank], MPI_DOUBLE, rho_root, sendsize, stride, MPI_DOUBLE, comm);
														
			//Calculate Pressure
			
			calcP();
			
			MPI_Allgatherv(p, sendsize[rank], MPI_DOUBLE, p_root, sendsize, stride, MPI_DOUBLE, comm);
						
		}
		
		//Calculate Forces and Acceleration
			
		calcFA();

		//Update values of x and v (inlined as they are simple functions)
	
		if (t == 0){
			
			//&a[0][0] etc are pointers to the first element in the contiguous 2D array
			
			cblas_daxpy((lengthloc)*2, 0.5 * dt, &a[0][0], 1, &v[0][0], 1);
			
			cblas_daxpy((lengthloc)*2, dt, &v[0][0], 1, &x[0][0], 1);
		
		}
	
		else {

			cblas_daxpy((lengthloc)*2, dt, &a[0][0], 1, &v[0][0], 1);
			
			cblas_daxpy((lengthloc)*2, dt, &v[0][0], 1, &x[0][0], 1);
		}
							
		//Validate Boundary conditions
		
		calcBC();
					
		//Calculate energy
		
		calcE();
	
		//Gather individual position and velocity contributions into x_root and v_root on all processes for use in next iteration
		
		MPI_Allgatherv(&x[0][0], sendsizex[rank], MPI_DOUBLE, &x_root[0][0], sendsizex, stridex, MPI_DOUBLE, comm);
		
		MPI_Allgatherv(&v[0][0], sendsizex[rank], MPI_DOUBLE, &v_root[0][0], sendsizex, stridex, MPI_DOUBLE, comm);

		//Sum energy contributions from each process onto root process via reduce and MPI_SUM
				
		MPI_Reduce(&Ek, &Ek_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
		
		MPI_Reduce(&Ep, &Ep_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
		
		MPI_Reduce(&Et, &Et_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
										
		//Only output energy values to file only if rank 0
		
		if (rank == 0){

					
			//Output Energy values into file

			EnergyOut << setprecision(5) << setw(15) << t << " "<< setprecision(5) << setw(15) << Ek_root << " " << setprecision(5) << setw(15) << Ep_root << " " << setprecision(5) << setw(15) << Et_root << endl;
	
		}
		
		

		//Reset Values for next iteration			
		
		//Reset density
		
		cblas_dscal(lengthloc, 0.0, rho, 1);
		
		//Reset Forces
		
		cblas_dscal(lengthloc * 2, 0.0, &Fv[0][0], 1);
		
		cblas_dscal(lengthloc * 2, 0.0, &Fp[0][0], 1);

		//Increment Time step
					
		t += dt;
	
				
					
	}
				
				
	//Output particle positions into file only in rank 0	
			
	if (rank == 0){

		
		
		ofstream PositionOut("output.txt", ios::out | ios::trunc);
				
		if(!PositionOut.good()){
					
			cout << "output.txt file could not be opened." << endl;
				
		}

		else {

			for (int i = 0; i < N; ++i){

				PositionOut << setprecision(5) << setw(20) << x_root[i][0] << setprecision(5) <<setw(20) << x_root[i][1]<< endl;

			}
			
			PositionOut.close();
		}
				
		cout << "Completed" << endl;
	}
	
		EnergyOut.close();
			
}