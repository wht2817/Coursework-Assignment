#include <iostream>
#include <cmath>
#include "cblas.h"
#include <ctime>
#include <random>
#include "CommandLineOptions.h"
#include <fstream>
//#include "mpi.h"
using namespace std;


//Forming the class first here (Easier to test)
class SPH {

	public:
			//Define Constructor 
			SPH(string configuration, double pdt, double pT, double ph){
				
				//Assign values to private variable
				dt = pdt;
				T  = pT;
				h = ph;
				hprime = 1.0/h;
				particles = configuration;
				
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
					N = 400;
				}
				else if (particles == "ic-block-drop"){
					N = 35;
				}
				else if (particles == "ic-droplet"){
					N = 51;
				}
				
				//Allocating Memory to variable arrays within constructor
				
				x = new double*[N];					  //Position of particles
				xpool = new double[N*2]();
				//r = new double**[N];                //Array of difference in displacement
			  
				v = new double*[N];                   //Velocity of particle
				vpool = new double[N*2]();
				//vij = new double**[N];              //Array of difference in velocities
			
				q = new double[N*N](); 				  //q array
			
				rho = new double[N]();                //Density
			
				p   = new double[N]();                //Pressure
			
				//Forces
				Fp = new double*[N];          	     //Pressure Force
				
				Fppool = new double[N*2]();            //
			
				Fv = new double*[N];         	     //Viscous Force
				
				Fvpool = new double[N*2]();
			
				Fg = new double*[N];             	 //Gravitational Force
				
				Fgpool = new double[N*2]();
			
				a  = new double*[N];              	 //Acceleration
				apool = new double[N*2]();
				
				for (int i = 0; i < N; ++i, xpool += 2, vpool += 2, apool += 2, Fppool += 2, Fvpool += 2, Fgpool += 2){
					x[i]  = xpool;
					v[i]  = vpool;
					a[i]  = apool;
					Fp[i] = Fppool;
					Fv[i] = Fvpool;
					Fg[i] = Fgpool;
				}
				
				
				for (int i = 0; i < N; i++){
//					x[i] = new double[2]();
					
//					v[i] = new double[2]();
					
//					Fp[i] = new double[2]();
//					
//					Fv[i] = new double[2]();
//					
//					Fg[i] = new double[2]();
					
//					a[i]  = new double[2]();
					
					//q[i] = new double[N]();
					
					//r[i] = new double*[N]();
					
					//vij[i] = new double*[N];
					
					for (int j = 0; j < N; j++){
						//r[i][j] = new double[2]();
						//vij[i][j] = new double[2];   // 2D array of coordinates (array of array of pointers)
					}
				}
				
				//Initialize Particles within constructor
				if (particles == "ic-one-particle"){
					x[0][0] = 0.5;
					x[0][1] = 0.5;
				}
				else if (particles == "ic-two-particles"){ //REMEMBER TO CHANGE BACK TO ORIGINAL TEST CASE
					x[0][0] = 0.5;
					x[0][1] = 0.5;
					x[1][0] = 0.509;
					x[1][1] = 0.5;
				}
				else if (particles == "ic-four-particles"){
					x[0][0] = 0.505;
					x[0][1] = 0.5;
					x[1][0] = 0.515;
					x[1][1] = 0.5;
					x[2][0] = 0.51;
					x[2][1] = 0.45;
					x[3][0] = 0.5;
					x[3][1] = 0.45;
				}
				else if (particles == "ic-dam-break"){
					for (int i = 0; i < 20; ++i){
						for (int j = 0; j < 20; ++j){
							//srand(time(0));
							//srand(time(0));
							x[i*20+j][0] = 0.01 + i*0.01 + rand()/(RAND_MAX*1000.0); //Random values for noise 
							x[i*20+j][1] = 0.01 + j*0.01 + rand()/(RAND_MAX*1000.0);
						}
					}
				}
				else if (particles == "ic-block-drop"){
					for (int i = 0; i < 5; ++i){
						for (int j = 0; j < 7; ++j){
							//srand(time(0));
							//srand(time(0));
							x[i*7+j][0] = 0.1 + i*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise
							x[i*7+j][1] = 0.3 + j*0.05 + rand()/(RAND_MAX*100.0);
							
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
					
					x[0][0] = 0.5;
					x[0][1] = 0.7;
					
					
					for (int i = 1; i < 5; ++i){
						for (int j = Npoints[i-1]+loopdummy; j < Npoints[i] + Npoints[i-1] + loopdummy; j++){
							
							x[j][0] = 0.5 + R[i]*cos(2.0 * (j-angledummy)*M_PI/Npoints[i]);
							x[j][1] = 0.7 + R[i]*sin(2.0 * (j-angledummy)*M_PI/Npoints[i]);
						}
						loopdummy += Npoints[i-1];
						angledummy += Npoints[i];
						
					}
					
					
				}
				
			}
			
			
			
			
			//Define Destructor
			~SPH(){
				
				//deallocate arrays
				delete[] q;
				delete[] p;
				delete[] rho;
				
				for (int i = 0; i < N; ++i){
					
//					delete[] x[i];
//					delete[] v[i];
//					delete[] Fg[i];
//					delete[] Fv[i];
//					delete[] Fp[i];
//					delete[] a[i];
				}
			delete[] x[0];
			delete[] v[0];
			delete[] a[0];
			delete[] x;
			delete[] v;
			delete[] Fg[0];
			delete[] Fv[0];
			delete[] Fp[0];
			delete[] Fg;
			delete[] Fv;
			delete[] Fp;
			delete[] a;

				
			}
			
			//Member Functions
			
			
			////////////////////
			//Print Functions
			////////////////////
			
			
			//Print x
			void printX(){
				for (int i = 0; i < N; ++i){
					cout << "Particle "<< i+1 << endl;
					cout << x[i][0] << " " << x[i][1] << endl;
				}
			}
			
			void printV(){
				for (int i = 0; i < N; ++i){
					cout << "Particle "<< i+1 << endl;
					cout << v[i][0] << " " << v[i][1] << endl;
				}
			}
			
			//Print q
//			void printQ(){
//				
//				cout << "Q Array: " << endl;
//				for (int i = 0; i < N; ++i){
//					for (int j = 0; j < N; ++j){
//						cout << q[i][j] << " ";
//					}
//					cout << endl;
//				}
//			}
			
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
					cout << p[i] << " ";
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
			
			
			//Calculate r array, q array and vij array
			void calcQRVIJ(){
				
				for (int i = 0; i < N; ++i){
					for(int j = 0; j < N; ++j){
						
						
						//Calculate r
//						r[0] = (x[i][0] - x[j][0]);
//						r[1] = x[i][1] - x[j][1];
//						
						//Calculate qij 
						q[i*N+j] = sqrt((x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) + (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]));
						
						
					}
				}
				
				cblas_dscal(N*N, hprime,q, 1);

			}
			
			//Calculate Density and pressure at the same time
			void calcRho(){
				
				//Define q_coeff
				
				double q_coeff;
				
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
						
						if (q[i* N + j] < 1){
							
							q_coeff = 1.0 - (q[i* N + j]*q[i* N + j]); //Precalculating q_coeff so this operation only needs to be done once instead of thrice
							
							rho[i] += coeff_rho * q_coeff * q_coeff * q_coeff;
							
							}
						}
					}
					
				}
			//Calculate Pressure
	
			void calcP(){
				
				for (int i = 0; i < N; ++i){
					p[i] = k*(rho[i] - rho_0);
				}
			}
			
			//Scale Mass
			void scaleMass(){
				m = N * rho_0/cblas_dasum(N, rho, 1);
			}
			
			//Calculate Pressure Force
			
			void calcFA(){
				
				//Initialize coefficients for repeated calculations
				double qcoeff;
				double FP_calc;
				double FV_calc;
				
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
							
							if ((q[i* N + j] < 1) && (i != j)){
							
								//Precalculate repeated coefficients
								qcoeff = (1.0-q[i* N + j]);
								
								FP_calc =  (coeff_p *(p[i]+p[j]) * qcoeff * qcoeff/(rho[j]*q[i* N + j]));
								
								Fp[i][0] += FP_calc * (x[i][0] - x[j][0]);

								Fp[i][1] += FP_calc * (x[i][1] - x[j][1]);
																			
						
//								cblas_daxpy(2, coeff_p * (p[i]+p[j]) * qcoeff * qcoeff/(rho[j]*q[i* N + j]),r,1, Fp[i], 1);
						
								FV_calc = (coeff_v * qcoeff / rho[j]);
								
								Fv[i][0] += FV_calc * (v[i][0] - v[j][0]);
								Fv[i][1] += FV_calc * (v[i][1] - v[j][1]);
								
								
								//Calculate Fv
//								cblas_daxpy(2, coeff_v * qcoeff /rho[j], vij, 1, Fv[i], 1);
							
							}
							
					}
						
						
						
						//Calculate Fg
						
						Fg[i][1] = -rho[i]*g;
						
						//Calculate Acceleration

						a[i][0] = (Fp[i][0] + Fv[i][0] + Fg[i][0])/rho[i];
						a[i][1] = (Fp[i][1] + Fv[i][1] + Fg[i][1])/rho[i];
						
//						Calculate Fp + Fv first (Value is logged into Fv)
//						cblas_daxpy(2, 1, Fp[i], 1, Fv[i], 1);
//						
//						Calculate Fp + Fv + Fg (Value is logged into Fv)
//						cblas_daxpy(2, 1, Fg[i], 1, Fv[i], 1);
//						
//						Copy value from Fv to a
//						cblas_dcopy(2, Fv[i], 1, a[i], 1);
//						
//						Scale value of a by 1/rho[i]
//						cblas_dscal(2, 1.0/rho[i], a[i], 1);
				}
			}
			
			//Enforce boundary condition
			void calcBC(){
				for (int i = 0; i < N; i++){
					
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
		
				
				for (int i = 0; i < N; i++){
				
					
					Ek += (v[i][0] * v[i][0] + v[i][1] * v[i][1]);//cblas_dnrm2(2, v[i], 1) * cblas_dnrm2(2, v[i], 1);
					
					Ep += x[i][1];
					
				}
				
				Ek *= 0.5*m;
				Ep *= m*g;
				Et  = Ek + Ep;
			}
			
			
			//Reset Value
			
			
			void solver(){
				//Open Files for output
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
					
					
					//Calculate Density and scale mass if first step
					calcQRVIJ();
					
					
					
					
					
					
					//Scale mass and calculate density again
					if (t == 0){
						
						calcRho();
						
						//calcP();   //Remember to delete 
						
						m = N * rho_0/cblas_dasum(N, rho, 1);
						
						
						//Reset Rho and coefficient before calculating again
						
						coeff_rho = m*4.0/(M_PI*h*h);
						
												
						cblas_dscal(N, 0.0, rho, 1);
						
						calcRho();
						
						calcP(); 
						
						coeff_p = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //scale the Fp value at the end
						coeff_v = -40.0*mu*m/(M_PI*h*h*h*h);
					}
					
					
					else{
						
						// Calculate Rho
						
						calcRho();
					
					
														
						//Calculate Pressure
						calcP();
						
					}

					
					//Caculate Forces and acceleration
					calcFA();

//					printForce();
					
					//Update values of x and v (inlined function)
					
					if (t == 0){
//						for (int i = 0; i < N ; ++i){
							
							//For first time step
							cblas_daxpy(N*2, dt/2.0, &a[0][0], 1, &v[0][0], 1); // x
							cblas_daxpy(N*2, dt, &v[0][0], 1, &x[0][0], 1);  // v
							
//						}
					}
					else {
						
						//For subsequent time steps
//						for (int i = 0; i < N; ++i){
							
						cblas_daxpy(N*2, dt, &a[0][0], 1, &v[0][0], 1); // v
						cblas_daxpy(N*2, dt, &v[0][0], 1, &x[0][0], 1); // x
						
//						}
					}
					
					//Validate Boundary conditions
					calcBC();
//					printX();
					//Calculate energy
					calcE();
					
//					cout << "Updated Positions and velocities: " << endl;
					
					//Display(Print energy values)
					//cout << "Energies: " << endl;
//					printMass();
					cout << "KE:" << Ek << endl;
					cout << "PE:" << Ep << endl;
					cout << "TE:" << Et << endl;
					
					//Output Energy values into file
					EnergyOut << t << " " << Ek << " " << Ep << " " << Et << endl;
					
					//Reset Values for next iteration			
					
					
					cblas_dscal(N, 0.0, rho, 1);
					
//					for (int i = 0; i < N; ++i){
						//for (int j = 0; j < N; ++j){
							
						//}
						//cblas_dscal(2, 0, q[i], 1);
						cblas_dscal(2 * N, 0, &Fv[0][0], 1);
						cblas_dscal(2 * N, 0, &Fp[0][0], 1);
						
//					}
					
					
					
					t += dt;
//					cout << endl;
				
					
				}
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
						PositionOut << x[i][0] << " " << x[i][1]<< endl;
					}
				}
				
				
				
				//Close files
				EnergyOut.close(); 
				PositionOut.close();
				
				
				cout << "Completed" << endl;
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
			
			double *r = new double[2]();
			
			double **v = nullptr;
			double *vpool = nullptr;
			//double ***vij = nullptr;
			double *vij = new double [2]();
			
			double *q = nullptr;
			
			double *rho = nullptr;
			
			double *p = nullptr;
			
			double **Fp = nullptr;
			
			double *Fppool  = nullptr;
			
			double **Fv = nullptr;
			
			double *Fvpool = nullptr;
			
			double **Fg = nullptr;
			
			double *Fgpool = nullptr;
	
			double **a = nullptr;
			double *apool = nullptr;
			
			//Energy
			double Ek;
			double Ep;
			double Et;
			
			//Particle Initialization String
			
			string particles;
			
			//Initialize calculation Coefficients
			
			
			double coeff_p = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //scale the Fp value at the end
			
			double coeff_v = -40.0*mu*m/(M_PI*h*h*h*h);
			
			double coeff_rho = m*4.0/(M_PI*h*h);
			
			double hprime;
			
	
	
	
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
	
	SPH test("ic-dam-break", 0.0001, 1, 0.01);
	test.solver();
	
	
	return 0;
}
