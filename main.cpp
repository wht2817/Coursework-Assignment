#include <iostream>
#include <cmath>
#include "cblas.h"
#include <ctime>
#include <random>
#include "CommandLineOptions.h"
#include <fstream>
#include "mpi.h"
using namespace std;


//Forming the class first here (Easier to test)
class SPH {

	public:
			//Define Constructor 
			SPH(string configuration, double pdt, double pT, double ph, MPI_Comm pcomm, int prank, int psize){
				
				//Assign values to private variable
				dt = pdt;
				T  = pT;
				h = ph;
				particles = configuration;
				
				this -> comm = pcomm;
				this -> rank = prank;
				this -> size = psize;
				//comm = pcomm;
				//rank = prank;
				//size = psize;
				
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
				
				//Only allocate memory if rank of proccess is less than N

				//Allocate memory for root variables in all proccesses
				//if (rank == 0){
					x_root = new double*[N];
					x_rootpool = new double[N*2];
					v_root = new double*[N];
					v_rootpool = new double[N*2];
					q_root = new double[N*N];
					rho_root = new double[N];
					p_root = new double[N];
					
					for (int i = 0; i < N; ++i, x_rootpool += 2, v_rootpool += 2){
						x_root[i] = x_rootpool;
						v_root[i] = v_rootpool;
					}
				//}
				
				//Allocating Memory to variable arrays within constructor
				
				x = new double*[N];					//Position of particles
				xpool = new double[N*2];
				//r = new double**[N];                //Array of difference in displacement
			
				v = new double*[N];               	//Velocity of particle
				vpool = new double[N*2];
				//vij = new double**[N];              //Array of difference in velocities
			
				q = new double[N*N]; 				//q array
				
				//qr = new double*[N];
			
				rho = new double[N];                //Density
			
				p   = new double[N];                //Pressure
			
				//Forces
				Fp = new double*[N];          	     //Pressure Force
			
				Fv = new double*[N];         	     //Viscous Force
			
				Fg = new double*[N];             	 //Gravitational Force
			
				a  = new double*[N];              	 //Acceleration
				
				//Try contiguous thing
				for (int i = 0; i < N; ++i, xpool += 2, vpool += 2){
					x[i] = xpool;
					v[i] = vpool;
				}
				
				for (int i = 0; i < N; i++){
					//x[i] = new double[2];
					
					//v[i] = new double[2];
					
					Fp[i] = new double[2];
					
					Fv[i] = new double[2];
					
					Fg[i] = new double[2];
					
					a[i]  = new double[2];
				}
				
				//Calculate start and end points of loops for each rank .
				
				int r = N % size;
				int k = (N - r) / size;
				
				if (rank < r) {
					k++;
					start = k * rank;
					end   = k * (rank + 1);
				}
				else{
					start = (k + 1) * r + k * (rank -r);
					end   = (k + 1) * r + k * (rank -r + 1);
				}
				
				//Initialize Particles within constructor
				if (particles == "ic-one-particle"){
					x_root[0][0] = 0.5;
					x_root[0][1] = 0.5;
//					if (rank == 0){
//						x[rank*N/size][0] = 0.5;
//						x[rank*N/size][0] = 0.5;
//					}
				}
				else if (particles == "ic-two-particles"){ //REMEMBER TO CHANGE BACK TO ORIGINAL TEST CASE
					x_root[0][0] = 0.5;
					x_root[0][1] = 0.5;
					x_root[1][0] = 0.509;
					x_root[1][1] = 0.5;
					
//					if(rank == 1){
//						x[start][0] = 0.509;
//						x[start][1] = 0.5;
//					}
//					else if (rank == 0){
//						x[start][0] = 0.5;
//						x[start][1] = 0.5;
//					}
					
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
					//Assign data for root
					for (int i = 0; i < 5; ++i){
						for (int j = 0; j < 5; ++j){
							//srand(time(0));
							//srand(time(0));
							x_root[i*5+j][0] = i*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise 
							x_root[i*5+j][1] = j*0.05 + rand()/(RAND_MAX*100.0);
						}
					}
					
				}
				else if (particles == "ic-block-drop"){
					for (int i = 0; i < 5; ++i){
						for (int j = 0; j < 7; ++j){
							//srand(time(0));
							//srand(time(0));
							x_root[i*7+j][0] = 0.1 + i*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise
							x_root[i*7+j][1] = 0.3 + j*0.05 + rand()/(RAND_MAX*100.0);
							
							//x[i*7+j][0] = 0.1 + i*(0.2/7.0) + rand()/(RAND_MAX*100.0);
							//x[i*7+j][1] = 0.3 + j*(0.3/7.0) + rand()/(RAND_MAX*100.0);
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
			
			
			
			
			//Define Destructor
			~SPH(){
				
				//deallocate arrays
				

				
			}
			
			//Member Functions
			
			
			////////////////////
			//Print Functions
			////////////////////
			
			
			//Print x
			void printX(){
				for (int i = 0; i < N; ++i){
					cout << "Particle "<< i+1 << endl;
					cout << x_root[i][0] << " " << x_root[i][1] << endl;
				}
			}
			
			void printV(){
				for (int i = 0; i < N; ++i){
					cout << "Particle "<< i+1 << endl;
					cout << v_root[i][0] << " " << v_root[i][1] << endl;
				}
			}
			
			//Print q
			void printQ(){
				
				cout << "Q Array: " << endl;
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
						cout << q_root[i*N + j] << " ";
					}
					cout << endl;
				}
			}
			
			//print r and vij
			void printVIJ(){
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
						cout << vij[0] << " " << vij[1];
					}
					cout << endl;
				}
				cout << endl;
			}
			//Print Rho
			void printRho(){
				cout << "Density Array:" << endl;
				
				for (int i = 0; i < N; ++i){
					cout << rho[i] << " ";
				}
				
				cout << endl;
			}
			
			void printP(){
				cout << "Pressure Array: " << endl;
				
				for (int i = 0; i < N; ++i){
					cout << p_root[i] << " ";
				}
				cout << endl;
			}
			
			void printMass(){
				
				cout << "Mass: " << m << endl;
			}
			
			void printForce(){
				
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
			
			void printA(){
				cout << "Acceleration:" << endl;
				for (int i = 0; i < N; ++i){
					cout << a[i][0] << " " << a[i][1] << endl;
				}
				cout << endl;
			}
			
			
			////////////////////////
			//Iteration Functions
			////////////////////////
			
			//Initialize by copying appropriate values from x_root to local x in each rank.
			void initX(){
				for (int i = start; i < end; ++i){
					cblas_dcopy(2, x_root[i], 1, x[i], 1);
				}
			}
			
			
			//Calculate r array, q array and vij array
			void calcQRVIJ(){
				
				for (int i = start; i < end; ++i){
					for(int j = 0; j < N; ++j){
						///////
						//Calc r and q
						///////
						
						//copy value from x[i] to r[i][j]
						cblas_dcopy(2, x[i], 1, r, 1);
						
						//Subtract xj from xi using blas, which records the value into r[i][j]
						cblas_daxpy(2, -1, x_root[j], 1, r, 1);
						
						//Calculate qij
						q[i*N+j] = cblas_dnrm2(2, r, 1)/h;
						
						
					}
				}
			}
			
			//Calculate Density and pressure at the same time
			void calcRho(){
				
				//Define coefficient outside of loop to save time
				double coeff = m*4.0/(M_PI*h*h);
				
				for (int i = start; i < end; ++i){
					for (int j = 0; j < N; ++j){
						
						if (q_root[i*N+j] < 1){
							
							rho[i] +=  coeff * (1 - (q_root[i*N+j]*q_root[i*N+j]))*(1 - (q_root[i*N+j]*q_root[i*N+j]))*(1 - (q_root[i*N+j]*q_root[i*N+j]));
							
							}
						}
					}
					
					//cblas_dscal(N, m*coeff , rho, 1); //Multiply matrix of densities by coefficient to save computational time
					
				}
			//Calculate Pressure
	
			void calcP(){
				
				for (int i = rank*N/size; i < (rank+1)*N/size; ++i){
					p[i] = k*(rho_root[i] - rho_0);
				}
			}
			
			//Scale Mass
			void scaleMass(){
				m = N * rho_0/cblas_dasum(N, rho_root, 1);
			}
			
			//Calculate Pressure Force
			
			void calcF(){
				
				//Calculate coefficient
				double coeff_p = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //scale the Fp value at the end
				double coeff_v;
				coeff_v = -40.0*mu*m/(M_PI*h*h*h*h);
				
				for (int i = start; i < end; ++i){
					for (int j = 0; j < N; ++j){
							
							if ((q_root[i*N+j] < 1) && (i != j)){
								//Calculate Pressure Force
								//copy value from x[i] to r[i][j]
								cblas_dcopy(2, x[i], 1, r, 1);
						
								//Subtract xj from xi using blas, which records the value into r[i][j]
								cblas_daxpy(2, -1, x_root[j], 1, r, 1);
								//Perform scaling on r vectors. Can be overwrriten as will not be using them anymore
						
								cblas_daxpy(2, coeff_p * (p_root[i]+p_root[j])*(1.0-q_root[i*N+j])*(1.0-q_root[i*N+j])/(rho_root[j]*q_root[i*N+j]),r,1, Fp[i], 1);
						
								//Calculate Viscous Force
								
								//copy value from v[i] to v[i][j]
								cblas_dcopy(2, v[i], 1, vij, 1);
						
								//Subtract vj from vi using blas, which records the value into r[i][j]
								cblas_daxpy(2, -1, v_root[j], 1, vij, 1);
							
								//Calculate Fv
								cblas_daxpy(2, coeff_v*(1-q_root[i*N+j])/rho_root[j], vij, 1, Fv[i], 1);
							
							}
							
					}
					//Calculate Gravitational Force
					Fg[i][1] = -rho_root[i]*g;
				}
			}
					
			
			//Enforce boundary condition
			void calcBC(){
				for (int i = start; i < end; ++i){
					
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
					else if (x[i][1] < (0.0 + h)){
						
						x[i][1] = 0.0 + h;
						v[i][1] = -e * v[i][1];
					}
					
				}
				
			}
			
			//Calculate Energy
			void calcE(){
				
				//Reset Energy values everytime they are calculated
				Ek = 0;
				Ep = 0;
				Et = 0;
				
				for (int i = start; i < end; ++i){
					
					Ek += pow(cblas_dnrm2(2, v[i], 1), 2);
					
					Ep += x[i][1];
					
				}
				
				Ek *= 0.5*m;
				Ep *= m*g;
				Et  = Ek + Ep;
			}
			
			//Calculate a
			void calcA(){
				for (int i = start; i < end; ++i){
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
			
			bool checksize(){
				if(rank >= N){
					return false;
					
				}
					
				else{
					return true;
				}
			}
			
			
			void solver(){
				
				//Initialize x array in indiviudal ranks. 
				initX();
				
				//if (rank == 0){
				
					ofstream EnergyOut("energy.txt", ios::out | ios::trunc);
				
					if(!EnergyOut.good()){
					
						cout << "Energy file could not be opened." << endl;
				
					}
					else{
						EnergyOut.precision(5);
						EnergyOut.width(15);
						EnergyOut << "Time" << " " << "KE" << " " << "PE" << " " << "TE" << endl;
					}
				
				//}
				//Initialize time loop
				double t = 0;
				cout << "Evaluating..."<< endl;
				while (t < T){
					
					//Initialize Positions of Particles
					//cout << "Step: " << t << endl;
					
					
					//cout << "Positions: " << endl;
					//printX();
					
					//cout << "Velocities: " << endl;
					//printV();
					
					
					calcQRVIJ();
					
					//Reducing q matrix to root and copy values from root to q
					
					
					MPI_Allreduce(q, q_root, N*N, MPI_DOUBLE, MPI_SUM, comm);
					
					
						
													
					//Calculate Pressure
					
					
					
					//Scale mass and calculate density again
					if (t == 0){
						
						calcRho();
						MPI_Allreduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, comm);
						//cblas_dcopy(N,rho_root, 1, rho, 1);
						
						calcP(); //REMEMBER TO REMOVE
						
						MPI_Allreduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, comm); //RMB TO REMOVE
						
						//cblas_dcopy(N, p_root, 1, p, 1);                     //RMB TO REMOVE
						
						scaleMass();
						
						//Reset Rho and coefficient before calculating again
												
						cblas_dscal(N, 0.0, rho, 1);
						
						calcRho();
						
						MPI_Allreduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, comm);
						//cblas_dcopy(N,rho_root, 1, rho, 1);
						
						//calcP();  //REMEMBER TO UNCOMMENT
						//MPI_Allreduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, comm); //RMB TO uncomment
						
						     
						
					}
					else{
						calcRho();
						MPI_Allreduce(rho, rho_root, N, MPI_DOUBLE, MPI_SUM, comm);
						//cblas_dcopy(N,rho_root, 1, rho, 1);
					
														
						//Calculate Pressure
						calcP();
						MPI_Allreduce(p, p_root, N, MPI_DOUBLE, MPI_SUM, comm);
						//cblas_dcopy(N, p_root, 1, p, 1);
						
					}
					//printRho();
					if (rank == 0){
						
					cout <<"Gathered Q: " << endl;	
					printQ();
					printP();
					
					}
					//cout << "Pressure: "<< endl;
					//printP();
					//cout << "Velocity Difference: " << endl;
					//printVIJ();
					//cout << "Mass: " << endl;
					//printMass();
					
					//Caculate Forces
					
					calcF();
					//calcFv();
					//calcFg();
					if (rank == 0){
					printForce();
					
					}
					//Calculate acceleration
					calcA();
					
					//cout << "Acceleration: " << endl;
					//printA();
					
					//Update values of x and t
					if (t == 0){
						for (int i = start; i < end; ++i){
							
							//For first time step
							cblas_daxpy(2, dt/2.0, a[i], 1, v[i], 1); // x
							cblas_daxpy(2, dt, v[i], 1, x[i], 1);  // v
							
						}
					}
					else {
						
						//For subsequent time steps
						for (int i = start; i < end; ++i){
						cblas_daxpy(2, dt, a[i], 1, v[i], 1); // v
						cblas_daxpy(2, dt, v[i], 1, x[i], 1); // x
						}
					}
					
					
					
					
					cout << "Positions before bc: " << endl;
					//printX();
					
					//Validate Boundary conditions
					calcBC();
					
					//Calculate energy
					calcE();
					cout << "Rank EK:" << rank << " " << Ek << endl;
					//Reduce Particle positions and velocities
					MPI_Allreduce(&x[0][0], &x_root[0][0], N*2, MPI_DOUBLE, MPI_SUM, comm);
					MPI_Allreduce(&v[0][0], &v_root[0][0], N*2, MPI_DOUBLE, MPI_SUM, comm);
					MPI_Reduce(&Ek, &Ek_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
					MPI_Reduce(&Ep, &Ep_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
					MPI_Reduce(&Et, &Et_root, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
					cout << "Updated Positions and velocities: " << endl;
					
					
					//Display(Print energy values)
					if (rank == 0){
						printX();
						printV();
					
						cout << "Energies: " << endl;
						cout << "KE:" << Ek_root << endl;
						cout << "PE:" << Ep_root << endl;
						cout << "TE:" << Et_root << endl;
					
						//Output Energy values into file
						EnergyOut << t << " " << Ek_root << " " << Ep_root << " " << Et_root << endl;
					}
					//Reset Values for next iteration			
					
					
					cblas_dscal(N, 0.0, rho, 1);
					cblas_dscal(2,0.0,r,1);
					cblas_dscal(2, 0.0, vij, 1);
					cblas_dscal(N*N, 0.0, q, 1);
					for (int i = 0; i < N; ++i){
						//for (int j = 0; j < N; ++j){
							
						//}
						//cblas_dscal(2, 0, q[i], 1);
						cblas_dscal(2, 0, Fv[i], 1);
						cblas_dscal(2, 0, Fp[i], 1);
						cblas_dscal(2, 0, Fg[i], 1);
					}
					
					
					
					//for (int i = 0; i < N; i++){
					//	cblas_dcopy(2, x_root[i], 1, x[i], 1);
					//	cblas_dcopy(2, v_root[i], 1, v[i], 1);
					//}
					
					
					t += dt;
					cout << endl;
				
					
				}
				if (rank == 0){
					//Output particle positions into file
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
	
	
	
	
	
	
	
	
			////////////////////////////
			//Getters and Setters
			////////////////////////////
			
			
	
	
	
	
	private:
	
	        /////////////////
			//Define constants
			/////////////////
			double k     = 2000.0; // Gas constant
			double rho_0 = 1000.0; // Resting density
			double mu    = 1.0;    // Viscosity
			double g     = 9.81;   // Gravitational acceleration (May not even need to define this)
			double h     = 0.01;   // Radius of influence
			double e     = 0.5;    // Coefficient of restitution
			
			/////////////////
			//Define Particle parameters
			/////////////////
			double m  = 1.0;           					// Initial Mass of particle
			
			double dt = 0;                                  //time step
			
			double T  = 0;                                  // Total Time
			
			int    N;                  					//No. of particles
			
			double **x = nullptr;
			double *xpool = nullptr;
			
			double *r = new double [2];
			
			double **v = nullptr;
			double *vpool = nullptr;
			//double ***vij = nullptr;
			double *vij = new double [2];
			
			double *q = nullptr;
			
			//double **qr = nullptr;
			
			double *rho = nullptr;
			
			double *p = nullptr;
			
			double **Fp = nullptr;
			
			double **Fv = nullptr;
			
			double **Fg = nullptr;
			
			double **a = nullptr;
			
			//Root Parameters
			double **x_root = nullptr;
			double *x_rootpool = nullptr;
			double **v_root = nullptr;
			double *v_rootpool = nullptr;
			double *q_root = nullptr;
			double Ek_root;
			double Ep_root;
			double Et_root;
			double *rho_root = nullptr;
			double *p_root   = nullptr;
			
			//Energy
			double Ek;
			double Ep;
			double Et;
			
			//Particle Initialization String
			
			string particles;
			
			//Initialize Communicator
			MPI_Comm comm;
			int rank;
			int size;
			
			//Initialize Loop Variables
			//Start of loop, end of loop
			int start;
			int end;
			
			
	
	
	
	
};

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
	
	//SPH test(InitialCondition, dt, T, h);
	
	//END OF BOOST STUFF
	
	MPI_Init(&argc, &argv);
	int rank, size;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	
	SPH test("ic-two-particles", 0.0001, 0.0002, 0.01, MPI_COMM_WORLD, rank, size);
	
	if (test.checksize()){
		test.solver();
		
	}
	
	
	MPI_Finalize();
	return 0;
}
