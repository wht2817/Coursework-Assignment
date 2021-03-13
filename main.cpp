#include <iostream>
#include <cmath>
#include "cblas.h"
#include <ctime>
#include <random>
#include "CommandLineOptions.h"

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
					N = 25;
				}
				else if (particles == "ic-block-drop"){
					N = 35;
				}
				else if (particles == "ic-droplet"){
					N = 51;
				}
				
				//Allocating Memory to variable arrays within constructor
				
				x = new double*[N];					//Position of particles
			
				r = new double**[N];                //Array of difference in displacement
			
				v = new double*[N];               	//Velocity of particle
				
				vij = new double**[N];              //Array of difference in velocities
			
				q = new double*[N]; 				//q array
			
				rho = new double[N];                //Density
			
				p   = new double[N];                //Pressure
			
				//Forces
				Fp = new double*[N];          	     //Pressure Force
			
				Fv = new double*[N];         	     //Viscous Force
			
				Fg = new double*[N];             	 //Gravitational Force
			
				a  = new double*[N];              	 //Acceleration
				
				for (int i = 0; i < N; i++){
					x[i] = new double[2];
					
					v[i] = new double[2];
					
					Fp[i] = new double[2];
					
					Fv[i] = new double[2];
					
					Fg[i] = new double[2];
					
					a[i]  = new double[2];
					
					q[i] = new double[N];
					
					r[i] = new double*[N];
					
					vij[i] = new double*[N];
					
					for (int j = 0; j < N; j++){
						r[i][j] = new double[2];
						vij[i][j] = new double[2];   // 2D array of coordinates (array of array of pointers)
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
					for (int i = 0; i < 5; ++i){
						for (int j = 0; j < 5; ++j){
							//srand(time(0));
							//srand(time(0));
							x[i*5+j][0] = i*0.05 + rand()/(RAND_MAX*100.0); //Random values for noise 
							x[i*5+j][1] = j*0.05 + rand()/(RAND_MAX*100.0);
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
			void printQ(){
				
				cout << "Q Array: " << endl;
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
						cout << q[i][j] << " ";
					}
					cout << endl;
				}
			}
			
			//print r and vij
			void printVIJ(){
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
						cout << vij[i][j][0] << " " << vij[i][j][1];
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
			void calcQRVIJ(int i){
				
					for(int j = 0; j < N; ++j){
						///////
						//Calc r and q
						///////
						
						//copy value from x[i] to r[i][j]
						cblas_dcopy(2, x[i], 1, r[i][j], 1);
						
						//Subtract xj from xi using blas, which records the value into r[i][j]
						cblas_daxpy(2, -1, x[j], 1, r[i][j], 1);
						
						//Calculate qij
						q[i][j] = cblas_dnrm2(2, r[i][j], 1)/h;
						
						//////
						//calc vij
						//////
						
						//copy value from v[i] to v[i][j]
						cblas_dcopy(2, v[i], 1, vij[i][j], 1);
						
						//Subtract vj from vi using blas, which records the value into r[i][j]
						cblas_daxpy(2, -1, v[j], 1, vij[i][j], 1);
						
						
					}
			}
			
			//Calculate Density and pressure at the same time
			void calcRho(int i){
				
				//Define coefficient outside of loop to save time
				//double coeff = 4.0/(M_PI*h*h);
				
					for (int j = 0; j < N; ++j){
						
						if (q[i][j] < 1){
							
							rho[i] +=  (1 - (q[i][j]*q[i][j]))*(1 - (q[i][j]*q[i][j]))*(1 - (q[i][j]*q[i][j]));
							
							}
						
						else{
							
							rho[i] += 0;
							}
						}
					
					
					//cblas_dscal(N, m*coeff_rho , rho, 1); //Multiply matrix of densities by coefficient to save computational time
					
				}
			//Calculate Pressure
	
			void calcP(){
				
				for (int i = 0; i < N; ++i){
					p[i] = k*(rho[i] - rho_0);
				}
			}
			
			//Scale Mass for initial condition
			void scaleMass(){
				m = N * rho_0/cblas_dasum(N, rho, 1);
			}
			
			//Calculate Pressure Force
			
			void calcFp(){
				
				//Calculate coefficient
				double coeff = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //scale the Fp value at the end
				
				for (int i = 0; i < N; ++i){
					for (int j = 0; j < N; ++j){
							
							if ((q[i][j] < 1) && (i != j)){
							
								//Perform scaling on r vectors. Can be overwrriten as will not be using them anymore
						
								cblas_daxpy(2, coeff * (p[i]+p[j])*(1.0-q[i][j])*(1.0-q[i][j])/(rho[j]*q[i][j]),r[i][j],1, Fp[i], 1);
						
								//Add value to appropriate element in Fp matrix
						
								//cblas_daxpy(2, 1, r[i][j], 1, Fp[i], 1);
							
							}
							//Do nothing if the above if condition is not satisfied
							//else{
								
								
								//cblas_daxpy(2, 0,r[i][j],1, Fp[i], 1);	
								//cblas_daxpy(2, 1, r[i][j], 1, Fp[i], 1);
							//}
					}
						//Scale Fp matrix by coefficient to get the correct value
						//cblas_dscal(2, coeff, Fp[i], 1);
				}
			}
			
			//Calculate Viscous Force
			void calcFv(){
				
				//Calculate Coefficient
				double coeff;
				coeff = -40.0*mu*m/(M_PI*h*h*h*h);
				
				for (int i = 0; i < N; i++){
					for (int j = 0; j < N; j++){
						
						if ((q[i][j] < 1) && (i != j)){
							
							//Calculate Fv
							cblas_daxpy(2, coeff*(1-q[i][j])/rho[j], vij[i][j], 1, Fv[i], 1);
						}
					}
					//cblas_dscal(2, coeff, Fv[i], 1); //Maybe can just pop this into the above calculation? cos its literally inlining it...
				}
				
			}
			
			//Calculate gravitational Force
			void calcFg(){
				for (int i = 0; i < N; i++){
					Fg[i][1] = -rho[i]*g;
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
					
					Ek += pow(cblas_dnrm2(2, v[i], 1), 2);
					
					Ep += x[i][1];
					
				}
				
				Ek *= 0.5*m;
				Ep *= m*g;
				Et  = Ek + Ep;
			}
			
			//Calculate a
			void calcA(){
				for (int i = 0; i < N; i++){
					//Calculate Fp + Fv first (Value is logged into Fv)
					cblas_daxpy(2, 1, Fp[i], 1, Fv[i], 1);
					//Calculate Fp + Fv + Fg (Value is logged into Fv)
					cblas_daxpy(2, 1, Fg[i], 1, Fv[i], 1);
					//Copy value from Fv to a
					cblas_dcopy(2, Fv[i], 1, a[i], 1);
					//Scale value of a by 1/rho[i]
					cblas_dscal(2, 1.0/rho[i], a[i], 1);
					
				}
			}
			
			//Reset Value
			
			
			void solver(){
				//initArray();
				
				
				//Initialize time loop
				double t = 0;
				while (t < T){
				
					
					//Initialize Positions of Particles
					cout << "Step: " << t << endl;
					
					
					cout << "Positions: " << endl;
					printX();
					
					cout << "Velocities: " << endl;
					printV();
					//Calculate Density
					for (int i = 0; i < N; ++i ){
						calcQRVIJ(i);
						calcRho(i);
					}
					
					//Scale Density (Inline)
					cblas_dscal(N, m*coeff_rho , rho, 1);
					
					//Calculate Pressure
					calcP();
					
					
					//Scale mass and calculate density again
					if (t == 0){
						scaleMass();
						
						//Reset Rho and coefficient before calculating again
						coeff_rho = 4.0/(M_PI*h*h);
						
						cblas_dscal(N, 0.0, rho, 1);
						for (int i = 0; i < N; ++i){
							calcRho(i);
						
						}
						cblas_dscal(N, m*coeff_rho , rho, 1);
						//calcP(); //DO I NEED TO DO THIS????? If need, just put it after the scaling (REMEMBER TO PUT THIS BACK!)
						
						
					}
					printRho();
					printQ();
					cout << "Pressure: "<< endl;
					printP();
					cout << "Velocity Difference: " << endl;
					printVIJ();
					cout << "Mass: " << endl;
					printMass();
					
					//Caculate Forces
					calcFp();
					calcFv();
					calcFg();
					
					printForce();
					//Calculate acceleration
					calcA();
					
					cout << "Acceleration: " << endl;
					printA();
					
					//Update values of x and t
					if (t == 0){
						for (int i = 0; i < N ; ++i){
							
							//For first time step
							cblas_daxpy(2, dt/2.0, a[i], 1, v[i], 1); // x
							cblas_daxpy(2, dt, v[i], 1, x[i], 1);  // v
							
						}
					}
					else {
						
						//For subsequent time steps
						for (int i = 0; i < N; ++i){
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
					
					cout << "Updated Positions and velocities: " << endl;
					
					printX();
					printV();
					
					//Display(Print energy values)
					cout << "Energies: " << endl;
					cout << "KE:" << Ek << endl;
					cout << "PE:" << Ep << endl;
					cout << "TE:" << Et << endl;
					
					//Reset Values for next iteration
					
					cblas_dscal(N, 0.0, rho, 1);
					for (int i = 0; i < N; ++i){
						for (int j = 0; j < N; ++j){
							cblas_dscal(2,0.0,r[i][j],1);
							cblas_dscal(2, 0.0, vij[i][j], 1);
						}
						cblas_dscal(2, 0, q[i], 1);
						cblas_dscal(2, 0, Fv[i], 1);
						cblas_dscal(2, 0, Fp[i], 1);
						cblas_dscal(2, 0, Fg[i], 1);
					}
					
					
					
					t += dt;
					cout << endl;
				
					
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
			
			double ***r = nullptr;
			
			double **v = nullptr;
			
			double ***vij = nullptr;
			
			double **q = nullptr;
			
			double *rho = nullptr;
			
			double *p = nullptr;
			
			double **Fp = nullptr;
			
			double **Fv = nullptr;
			
			double **Fg = nullptr;
			
			double **a = nullptr;
			
			
			//Energy
			double Ek;
			double Ep;
			double Et;
			
			//Particle Initialization String
			
			string particles;
			
			//Initialize calculation Coefficients
			double coeff_rho = 4.0/(M_PI*h*h);
			
			
	
	
	
	
};

int main(int argc, char **argv)
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
	
	SPH test("ic-two-particles", 0.0001, 0.0002, 0.01);
	test.solver();
	
	
	return 0;
}
