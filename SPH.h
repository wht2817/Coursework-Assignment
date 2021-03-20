#ifndef SPH_H
#define SPH_H
#include "mpi.h"

using namespace std;


class SPH
{
public:

	//Constructor which initializes initial condtion and allocates data based on rank in MPI.
	SPH(string configuration, double pdt, double pT, double ph, MPI_Comm pcomm, int prank, int psize);
	
	//Destructor
	~SPH();
	
	/////
	//Print Functions
	/////
		
	void printX();
		
	//Print Velocity
	void printV();
	
	//Print q matrix
	void printQ(); 
	
	//Print Density
	void printRho();
	
	//Print Pressure
	void printP();
	
	//Print Mass
	void printMass();
	
	//Print Forces
	void printForce();
		
	//Print Acceleration
	void printA();
	
	
	/////
	//Iteration Functions
	/////
	
	//Initialize parameters to aid in distributing work among proccesses
	void initMPIVariables();
	
	//Initialize Particles
	
	void initRootParticles();
	
	//Initialize Local Particles
	void initLocParticles();
	
	//Initialize values in local x matrix from x_root matrix
	
	void initX();
			
	//Precalculate array of q values
	void calcQRVIJ();
	
	//Calculate Density array
	void calcRho();
	
	//Scale Mass
	void scaleMass();
	
	//Calculate Pressure
	void calcP();

	//Calculate All Forces		
	void calcFA();
			
	//Enforce boundary condition
	void calcBC();
			
	//Calculate Energy
	void calcE();
	
	//Calculate Acceleration
	void calcA();
	
	
	
	/////
	//Main Solver
	/////
	
	//Solves SPH algorithm
	void solver();
	
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
	double m  = 1.0;           					    // Initial Mass of particle
			
	double dt = 0.0001;                             //Default timestep
			
	double T  = 5;                                  // Default total time
			
	int    N;                  						//No. of particles
	
	/////
	//Process Local arrays
	/////
	
	//Pointers for contiguous 2D x array		
	double **x = nullptr;
	
	double *xpool = nullptr;
	
	//Pointer for array rij = xi - xj
	double *r = new double [2]();
	
	//Pointers for contiguous 2D v array
	double **v = nullptr;
	
	double *vpool = nullptr;
	
	//Pointer for array vij = vi - vj;
	double *vij = new double [2]();
	
	//Pointer for array storing q values
	double *q = nullptr;
	
	//Pointer for array storing density values
	double *rho = nullptr;
	
	//Pointer for array storing pressure values
	double *p = nullptr;
	
	//Pointer for array storing Forces
	double **Fp = nullptr;
	
	double **Fv = nullptr;
	
	double **Fg = nullptr;
	
	//Pointer for array storing Acceleration
	double **a = nullptr;
	
	double *apool = nullptr;
	
	//Energy
	double Ek;
	
	double Ep;
	
	double Et;
	
	/////
	//Root parameters
//	Contain all the values of x, v, q, rho and p and total Ek, Ep and Et.
//	Only the energies are stored exclusively on the root, the rest of the parameters can be accessed by all ranks
//	so that calculation of forces, density and pressure can take place within each rank.
	/////
	
	double **x_root = nullptr;
	
	double *x_rootpool = nullptr;
	
	double **v_root = nullptr;
	
	double *v_rootpool = nullptr;
	
	double Ek_root;
	
	double Ep_root;
	
	double Et_root;
	
	double *rho_root = nullptr;
	
	double *p_root   = nullptr;
			
	
			
	/////
	//Coefficients for calculating density and forces
	/////
	
//	Only calculated twice per rank in whole code to save time
			
	double coeff_rho = m*4.0/(M_PI*h*h);            //Density
	
	double coeff_p = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //Pressure force
	
	double coeff_v = -40.0*mu*m/(M_PI*h*h*h*h);     //Viscous force
	
	double hprime;                                  // 1.0/h
			
	/////		
	//Particle Initialization String
	/////
	
	string particles;
	
	/////
	//MPI communicator variables
	/////
	
	MPI_Comm comm;
	
	int rank;
	
	int size;
	
	//Loop variables to evenly distribute work amongst processes
	
	int start;
	
	int finish;
	
	int lengthloc;
	
	int *sendsizex;
	
	int *stridex;
	
	int *sendsize;
	
	int *stride;
	
	
};

#endif // SPH_H