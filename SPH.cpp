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
/////
//Constructor
/////
SPH::SPH(string configuration, double pdt, double pT, double ph, MPI_Comm pcomm, int prank, int psize){
				
	//Assign values to private variable
	dt = pdt;
	T  = pT;
	h = ph;
	particles = configuration;
				
	this -> comm = pcomm;
	this -> rank = prank;
	this -> size = psize;
	
	//Assign number of particles based on option chosen by user
	if (particles == "ic-one-particle"){
		N = 1;
	}
	else if (particles == "ic-two-particles"){
		N = 2;
	}
	else if (particles == "ic-four-particles"){
		N = 4;
	}
	else if (particles == "ic-dam-break"){
		N = 25;
	}
	else if (particles == "ic-block-drop"){
		N = 35;
	}
	else if (particles == "ic-droplet"){
		N = 51;
	}
	
	

	//Allocate memory for root(shared) variables in all proccesses
	
	x_root     = new double*[N];
	
	x_rootpool = new double[N*2]();
	
	v_root     = new double*[N];
	
	v_rootpool = new double[N*2]();
	
	q_root     = new double[N*N]();
	
	rho_root   = new double[N]();
	
	p_root     = new double[N]();
	
	//Forming contiguous 2D array for x_root and v_root
	
	for (int i = 0; i < N; ++i, x_rootpool += 2, v_rootpool += 2){
		x_root[i] = x_rootpool;
		v_root[i] = v_rootpool;
	}
				
				
	//Allocating Memory to local variable arrays
				
	x = new double*[N];					//Position of particles
	
	xpool = new double[N*2]();	
			
	v = new double*[N];               	//Velocity of particle
	
	vpool = new double[N*2]();
			
	q = new double[N*N](); 				//q array
				
	rho = new double[N]();                //Density
			
	p   = new double[N]();                //Pressure
			
	//Forces
	Fp = new double*[N];          	     //Pressure Force
			
	Fv = new double*[N];         	     //Viscous Force

	Fg = new double*[N];             	 //Gravitational Force
			
	a  = new double*[N];              	 //Acceleration
				
	//Form contiguous 2D array for x and v
	for (int i = 0; i < N; ++i, xpool += 2, vpool += 2){
		x[i] = xpool;
		v[i] = vpool;
	}
				
	for (int i = 0; i < N; i++){
			
	Fp[i] = new double[2]();
				
	Fv[i] = new double[2]();
					
	Fg[i] = new double[2]();
					
	a[i]  = new double[2]();
	}
				
	//Calculate start and end points of loops for each rank, and consequently, positions where values are allocated in x from x_root
		
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
				
	//Initialize Particle positions in x_root
	
	if (particles == "ic-one-particle"){
		
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.5;
	}
	
	else if (particles == "ic-two-particles"){ //REMEMBER TO CHANGE BACK TO ORIGINAL TEST CASE
	
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.5;
		x_root[1][0] = 0.5;
		x_root[1][1] = h;
			
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

		for (int i = 0; i < 5; ++i){
			for (int j = 0; j < 5; ++j){

				x_root[i*5+j][0] = i*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise 
				x_root[i*5+j][1] = j*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise

			}
		}
					
	}
	
	else if (particles == "ic-block-drop"){

		for (int i = 0; i < 5; ++i){
			for (int j = 0; j < 7; ++j){

				x_root[i*7+j][0] = 0.1 + i*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise
				x_root[i*7+j][1] = 0.3 + j*0.05 + rand()/(RAND_MAX*100.0);
				
			}
		}
	}

	else if (particles == "ic-droplet"){
				
		//Initialize plotting parameters
		int Npoints[5] = {1, 5, 10, 15, 20}; //No. of points on each circumference to give evenly spaced points
		double R[5]    = {0, 0.025, 0.05, 0.075, 0.1}; // Radii at which the circle will be formed.
		
		//initialize dummy variable to aid in populatin array with evenly spaced points
		int loopdummy = 0;
		double angledummy = 0.0;
					
		x_root[0][0] = 0.5;
		x_root[0][1] = 0.7;
					
					
		for (int i = 1; i < 5; ++i){
			for (int j = Npoints[i-1]+loopdummy; j < Npoints[i] + Npoints[i-1] + loopdummy; j++){
				
				x_root[j][0] = 0.5 + R[i]*cos(2.0 * (j-angledummy)*M_PI/Npoints[i]);
				x_root[j][1] = 0.7 + R[i]*sin(2.0 * (j-angledummy)*M_PI/Npoints[i]);
			}

			loopdummy += Npoints[i-1];
			angledummy += Npoints[i];
				
		}
					
					
				}
				
}

SPH::~SPH()
{
	//deallocate memory
//	delete[] x;
//	delete[] xpool;
//	delete[] x_root;
//	delete[] x_rootpool;
//	delete[] v;
//	delete[] vpool;
//	delete[] v_rootpool;
//	delete[] v_root;
	delete[] q;
	delete[] q_root;
	delete[] rho;
	delete[] rho_root;
	delete[] p;
	delete[] p_root;
//	
	for (int i = 0; i < N; ++i){
		
		delete[] Fp[i];
		delete[] Fv[i];
		delete[] Fg[i];
		delete[] a[i];
		

	}
	
	delete[] Fp;
	delete[] Fv;
	delete[] Fg;
	delete[] a;
	
}

			
/////
//Member Functions
/////
			
			
////////////////////
//Print Functions
////////////////////
			
			
//Print Positions
void SPH::printX(){

	for (int i = 0; i < N; ++i){
	
		cout << "Particle "<< i+1 << endl;
	
		cout << x_root[i][0] << " " << x_root[i][1] << endl;

	}
}

//Print Velocity			
void SPH::printV(){

	for (int i = 0; i < N; ++i){

		cout << "Particle "<< i+1 << endl;
		cout << v_root[i][0] << " " << v_root[i][1] << endl;

	}
}
			
//Print q array

void SPH::printQ(){
				
	cout << "Q Array: " << endl;

	for (int i = 0; i < N; ++i){
		for (int j = 0; j < N; ++j){

			cout << q_root[i*N + j] << " ";
		}

		cout << endl;
	}
}
			
//Print Rho
void SPH::printRho(){
	
	cout << "Density Array:" << endl;
	
	for (int i = 0; i < N; ++i){

		cout << rho[i] << " ";

	}
				
	cout << endl;
}

//Print Pressure
void SPH::printP(){
	cout << "Pressure Array: " << endl;
	
	for (int i = 0; i < N; ++i){

		cout << p_root[i] << " ";
	}

	cout << endl;
}

//Print Mass
void SPH::printMass(){
	
	cout << "Mass: " << m << endl;
}

//Print Forces
void SPH::printForce(){
				
	cout << "Pressure Force: " << endl;

	for (int i = 0; i < N; ++i){

		cout << Fp[i][0] << " " << Fp[i][1] << endl;

	}
				
	cout << "Viscous Force: " << endl;

	for (int i = 0; i < N; ++i){

					cout << Fv[i][0] << " " << Fv[i][1] << endl;

	}
				
	cout << "Gravitational Force: " << endl;

	for (int i = 0; i < N; ++i){

		cout << Fg[i][0] << " " << Fg[i][1] << endl;

	}

	cout << endl;

}

//Print acceleration
			
void SPH::printA(){

	cout << "Acceleration:" << endl;

	for (int i = 0; i < N; ++i){

		cout << a[i][0] << " " << a[i][1] << endl;

	}

cout << endl;

}
			
			
/////
//Iteration Functions
/////
			
//Initialize values in local x matrix from x_root matrix

void SPH::initX(){
	
	for (int i = start; i < finish; ++i){

		cblas_dcopy(2, x_root[i], 1, x[i], 1);

	}
}
			
			
//Precalculate array of q values

void SPH::calcQRVIJ(){
	
	double hprime = 1.0 / h; 
	
	for (int i = start; i < finish; ++i){
		for(int j = 0; j < N; ++j){
					
			//copy value from x[i] to r
			cblas_dcopy(2, x[i], 1, r, 1);
						
			//Subtract xj from xi using blas, which records the value into r
			cblas_daxpy(2, -1, x_root[j], 1, r, 1);
						
			//Calculate qij
			q[i*N+j] = cblas_dnrm2(2, r, 1) * hprime;
						
						
		}
	}
}

//Calculate Density

void SPH::calcRho(){
				
	//Define q_coeff
				
	double q_coeff;
				
				
	for (int i = start; i < finish; ++i){
		for (int j = 0; j < N; ++j){
						
			if (q_root[i*N+j] < 1){
				
				q_coeff = 1.0 - (q_root[i* N + j]*q_root[i* N + j]);
			
				rho[i] += coeff_rho * q_coeff * q_coeff * q_coeff;
							
				}
			}
		}
					
					
}
//Calculate Pressure
	
void SPH::calcP(){
				
	for (int i = start; i < finish; ++i){

		p[i] = k*(rho_root[i] - rho_0);

	}
}
			
//Scale Mass

void SPH::scaleMass(){

	m = N * rho_0/cblas_dasum(N, rho_root, 1);

}
			
//Calculate Pressure Force
			
void SPH::calcFA(){
				
	
	double qcoeff;
	
	for (int i = start; i < finish; ++i){
		for (int j = 0; j < N; ++j){
				
			if ((q_root[i*N+j] < 1) && (i != j)){
				
				//Precalculate repeated coefficients
				qcoeff = (1.0-q_root[i* N + j]);
				
				//Calculate Pressure Force
				
				//copy value from x[i] to r
				cblas_dcopy(2, x[i], 1, r, 1);
						
				//Subtract xj from xi using blas, which records the value into r[i][j]
				cblas_daxpy(2, -1, x_root[j], 1, r, 1);
				
				//Perform scaling on r vectors. Can be overwrriten as will not be using them anymore
				cblas_daxpy(2, coeff_p * (p_root[i]+p_root[j]) * qcoeff * qcoeff/(rho_root[j]*q_root[i* N + j]), r , 1, Fp[i], 1);
						
				//Calculate Viscous Force
								
				//copy value from v[i] to v[i][j]
				cblas_dcopy(2, v[i], 1, vij, 1);
						
				//Subtract vj from vi using blas, which records the value into r[i][j]
				cblas_daxpy(2, -1, v_root[j], 1, vij, 1);
							
				//Calculate Fv
				cblas_daxpy(2, coeff_v * qcoeff/rho_root[j], vij, 1, Fv[i], 1);
							
			}
						
		}
				//Calculate Gravitational Force
				
				Fg[i][1] = -rho_root[i]*g;
				
				
				//Calculate Acceleration
				
				//Calculate Fp + Fv first (Value is logged into Fv)
				cblas_daxpy(2, 1, Fp[i], 1, Fv[i], 1);

				//Calculate Fp + Fv + Fg (Value is logged into Fv)
				cblas_daxpy(2, 1, Fg[i], 1, Fv[i], 1);

				//Copy value from Fv to a
				cblas_dcopy(2, Fv[i], 1, a[i], 1);

				//Scale value of a by 1/rho[i]
				cblas_dscal(2, 1.0/rho_root[i], a[i], 1);
	}
}
					
			
//Enforce boundary condition

void SPH::calcBC(){

	for (int i = start; i < finish; ++i){
					
		//Right bc
		if (x[i][0] > (1.0 - h)){
			
			x[i][0] = 1.0 - h;
			v[i][0] = -e * v[i][0];
		}

		//Left bc
		else if (x[i][0] < (0.0 + h)){
						
			x[i][0] = 0.0 + h;
			v[i][0] = -e * v[i][0];
		}
					
		//Top bc
		if (x[i][1] > (1.0 - h)){
						
			x[i][1] = 1.0 - h;
			v[i][1] = -e * v[i][1];
		}

		//Bottom bc
		else if (x[i][1] < (0.0 + h)){
						
				x[i][1] = 0.0 + h;
				v[i][1] = -e * v[i][1];
		}
			
	}
		
}
			
			//Calculate Energy
void SPH::calcE(){
				
	//Reset Energy values everytime they are calculated
	Ek = 0;
	Ep = 0;
	Et = 0;
				
	for (int i = start; i < finish; ++i){
					
		Ek += pow(cblas_dnrm2(2, v[i], 1), 2);
		
		Ep += x[i][1];
					
	}
				
	Ek *= 0.5*m;
	Ep *= m*g;
	Et  = Ek + Ep;
}
			
//Calculate a
void SPH::calcA(){

	for (int i = start; i < finish; ++i){

		//Calculate Fp + Fv first (Value is logged into Fv)
		cblas_daxpy(2, 1, Fp[i], 1, Fv[i], 1);

		//Calculate Fp + Fv + Fg (Value is logged into Fv)
		cblas_daxpy(2, 1, Fg[i], 1, Fv[i], 1);

		//Copy value from Fv to a
		cblas_dcopy(2, Fv[i], 1, a[i], 1);

		//Scale value of a by 1/rho[i]
		cblas_dscal(2, 1.0/rho_root[i], a[i], 1);
					
	}
}

//Check if No. of proccsses exceeds no. of particles.
		
bool SPH::checksize(){

	if(size > N){

		if (rank == 0){

			cout << "Number of proccesses cannot exceed " << N << " for "<< particles << " case." << endl;

		}

		return false;
					
	}
					
	else{

		return true;

	}
}

/////
//Main Solver
/////			
			
void SPH::solver(){
				
	//Distribute positions to x arrays in individual ranks 
	initX();
				
	
	//Open energy file for output at every time interval
	
	ofstream EnergyOut("energy.txt", ios::out | ios::trunc);
				
		if(!EnergyOut.good()){
					
			cout << "Energy file could not be opened." << endl;
				
		}

		else{

			EnergyOut.precision(5);
			EnergyOut.width(15);
			EnergyOut << "Time" << " " << "KE" << " " << "PE" << " " << "TE" << endl;
	
	}
				
	//Initialize time loop
	double t = 0;

	cout << "Evaluating..."<< endl;

	while (t < T){
					
					
		//cout << "Step: " << t << endl;
		
//		if (rank == 0){
//		cout << "Initial Positions: " << endl;
//		printX();
//					
//		}
					
		//cout << "Velocities: " << endl;
		//printV();
					
		//Calculate q array in each process
			
		calcQRVIJ();
					
		//Reduce q array via MPI_SUM to q_root matrix in all processes for use in calculating density and forces
		
//		MPI_Reduce(q, q_root, N*N, MPI_DOUBLE, MPI_SUM, 0, comm);
//		MPI_Bcast(q_root, N*N, MPI_DOUBLE, 0, comm);
		MPI_Allreduce(q, q_root, N*N, MPI_DOUBLE, MPI_SUM, comm);
					
													
		//If First step, scale mass and recalculate density and force coefficients
		
		if (t == 0){
			
			//Calculate Density
			calcRho();

			//Reduce individual density matrices in each process to rho_root array for use in force calculations
			
//			MPI_Reduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, 0, comm);
//			MPI_Bcast(rho_root, N, MPI_DOUBLE, 0, comm);
			MPI_Allreduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, comm);
						
						
			//calcP(); //REMEMBER TO REMOVE
						
			//MPI_Allreduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, comm); //RMB TO REMOVE
						
						
			//Scale Mass			
			m = N * rho_0/cblas_dasum(N, rho_root, 1);
						
			//Reset Rho and recalculate coefficient before calculating again
						
			coeff_rho = m*4.0/(M_PI*h*h);
						
			cblas_dscal(N, 0.0, rho, 1);
						
			
			//Recalculate Rho and P and reduce to root arrays
			calcRho();
//			MPI_Reduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, 0, comm);
//			MPI_Bcast(rho_root, N, MPI_DOUBLE, 0, comm);						
			MPI_Allreduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, comm);
						
			calcP();  //REMEMBER TO UNCOMMENT
			
//			MPI_Reduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, 0, comm);
//			MPI_Bcast(p_root, N, MPI_DOUBLE, 0, comm);
			MPI_Allreduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, comm); //RMB TO uncomment
						
			//Rescale calculate coefficients after mass has been calculated
						
			coeff_p   = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //scale the Fp value at the end

			coeff_v   = -40.0*mu*m/(M_PI*h*h*h*h);
						     
						
		}

		else{
			
			//Calculate Density
			
			calcRho();
//			MPI_Reduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, 0, comm);
//			MPI_Bcast(rho_root, N, MPI_DOUBLE, 0, comm);
			MPI_Allreduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, comm);

														
			//Calculate Pressure
			calcP();
//			MPI_Reduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, 0, comm);
//			MPI_Bcast(p_root, N, MPI_DOUBLE, 0, comm);
			MPI_Allreduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, comm);

						
		}
		//printRho();
//		if (rank == 0){
//			
//		cout <<"Gathered Q: " << endl;	
//		printQ();
//		printP();
//					
//		}
		//cout << "Pressure: "<< endl;
		//printP();
		//cout << "Velocity Difference: " << endl;
		//printVIJ();
		//cout << "Mass: " << endl;
//		printMass();
		
		//Caculate Forces
					
		calcFA();

		//Calculate acceleration

//		calcA();
					
//		if (rank == 0){
//		cout << "Acceleration: " << endl;
//		printA();
//		}

		//Update values of x and t
	
		if (t == 0){

			for (int i = start; i < finish; ++i){

				//For first time step
				cblas_daxpy(2, dt/2.0, a[i], 1, v[i], 1); // v
				cblas_daxpy(2, dt, v[i], 1, x[i], 1);  // x
							
			}
		}
	
		else {

			//For subsequent time steps

			for (int i = start; i < finish; ++i){
				cblas_daxpy(2, dt, a[i], 1, v[i], 1); // v
				cblas_daxpy(2, dt, v[i], 1, x[i], 1); // x
			}
		}
					
					
//					
//		if (rank == 0){
//		cout << "Positions and velocitybefore bc: " << endl;
//		printX();
//		printV();
//		}
					
		//Validate Boundary conditions
		calcBC();
					
		//Calculate energy
		calcE();
					
//		cout << "Rank EK:" << rank << " " << Ek << endl;
	
		//Update x_root and v_root by using reduce and bcast so that their values can be used in calculation for next iteration
		
//		MPI_Reduce(&x[0][0], &x_root[0][0], N*2, MPI_DOUBLE, MPI_SUM, 0, comm);
//		MPI_Bcast(&x_root[0][0], N*2, MPI_DOUBLE, 0, comm);
//		MPI_Reduce(&v[0][0], &v_root[0][0], N*2, MPI_DOUBLE, MPI_SUM, 0, comm);
//		MPI_Bcast(&v_root[0][0], N*2, MPI_DOUBLE, 0, comm);
		
		MPI_Allreduce(&x[0][0], &x_root[0][0], N*2, MPI_DOUBLE, MPI_SUM, comm);
		MPI_Allreduce(&v[0][0], &v_root[0][0], N*2, MPI_DOUBLE, MPI_SUM, comm);

		//Sum energy contributions from each process onto root process
	
		MPI_Reduce(&Ek, &Ek_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
		MPI_Reduce(&Ep, &Ep_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
		MPI_Reduce(&Et, &Et_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

//		cout << "Updated Positions and velocities: " << endl;
					
					
		//Display(Print energy values)
		if (rank == 0){
//			printX();
//			printV();
					
//			cout << "Energies: " << endl;
			cout << "KE:" << Ek_root << endl;
			cout << "PE:" << Ep_root << endl;
			cout << "TE:" << Et_root << endl;
					
			//Output Energy values into file
			EnergyOut << t << " " << Ek_root << " " << Ep_root << " " << Et_root << endl;
	
		}

		//Reset Values for next iteration			
					
					
		cblas_dscal(N, 0.0, rho, 1);
	
//		cblas_dscal(2,0.0,r,1);
	
//		cblas_dscal(2, 0.0, vij, 1);
	
//		cblas_dscal(N*N, 0.0, q, 1);
	
		for (int i = start; i < finish; ++i){
			cblas_dscal(2, 0, Fv[i], 1);
			cblas_dscal(2, 0, Fp[i], 1);
//			cblas_dscal(2, 0, Fg[i], 1);
		}
					

		//Increment Time step
					
		t += dt;
	
				
					
	}
				
				
				
	if (rank == 0){

		//Output particle positions into file only in rank 0
		
		ofstream PositionOut("output.txt", ios::out | ios::trunc);
				
		if(!PositionOut.good()){
					
			cout << "output.txt file could not be opened." << endl;
				
		}

		else {
						
			PositionOut.precision(5);
			PositionOut.width(15);
				
			PositionOut<< "x coordinate" << " " <<"y coordinate"<< endl; 

			for (int i = 0; i < N; ++i){

				PositionOut << x_root[i][0] << " " << x_root[i][1]<< endl;

			}
		}
				
				
				
		//Close files

		EnergyOut.close(); 

		PositionOut.close();
				
				
		cout << "Completed" << endl;
	}
			
}