#ifndef SPH_H
#define SPH_H
#include "mpi.h"

using namespace std;


class SPH
{
public:
	
	//Delete default constructor (Do not allow object to be created without parameters)
	SPH() = delete;
	
	//Constructor which initializes initial condtion and allocates data based on rank in MPI.
	SPH(string configuration, double pdt, double pT, double ph, MPI_Comm pcomm, int prank, int psize);
	
	//Destructor
	~SPH();
		
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
	//Define private constants
	/////////////////
	double k     = 2000.0; // Gas constant
	double rho_0 = 1000.0; // Resting density
	double mu    = 1.0;    // Viscosity
	double g     = 9.81;   // Gravitational acceleration (May not even need to define this)
	double e     = 0.5;    // Coefficient of restitution
			
	/////////////////
	//Define Particle parameters
	/////////////////
	double m  = 1.0;           					    // Initial Mass of particle
			
	double dt = 0.0001;                             // Default timestep
			
	double T  = 5;                                  // Default total time
			
	int    N  = 0;                  				// No. of particles
	
	double h  = 0.01;   							// Radius of influence
	
	/////
	//Process Local arrays
	/////
	
	/* Local contributions of each process for each variable
	 * These represent x_i, v_i etc
	 */
	
			
	double **x = nullptr;							// Array of pointers to xpool
	
	double *xpool = nullptr;						// Pointer to particle local positions
	
	double *r = new double [2]();					// Pointer to array r = xi - xj
	
	double **v = nullptr;							// Array of pointers to vpool
	
	double *vpool = nullptr;						// Pointer to particle local velocities
	
	double *q = nullptr;							// Pointer to array of local q values
	
	double *rho = nullptr;							// Pointer to array of local densities

	double *p = nullptr;							// Pointer to array of local pressures
	
	double **Fp = nullptr;							// Array of pointers to Fppool
	
	double *Fppool = nullptr;						// Pointer to local particle Fp
	
	double **Fv = nullptr;							// Array of pointers to Fvpool
	
	double *Fvpool = nullptr;						// Pointer to local particle Fv
	
	double **Fg = nullptr;							// Array of pointers to Fgpool
	
	double *Fgpool = nullptr;						// Pointer to local particle Fg

	double **a = nullptr;							// Array of pointers to apool
	
	double *apool = nullptr;						// Pointer to local particle acceleration
	
	double Ek;										// Local contribution of kinetic energy
	
	double Ep;										// Local contribution of potential energy
	
	double Et;										// Local contribution of total energy
	
	/////
	//Root parameters and arrays
	/////
 
	/*
	*Contains all the values of x, v, q, rho and p and total Ek, Ep and Et.
	*Only the energies are stored exclusively on the root, the rest of the parameters can be accessed by all ranks
	*so that calculation of forces, density and pressure can take place within each rank.
	*These represent x_j, v_j etc...
	*/
	
	double **x_root = nullptr;						// Array of pointers to x_rootpool
	
	double *x_rootpool = nullptr;					// Pointer to all particle positions
	
	double **v_root = nullptr;						// Array of pointers to v_rootpool
	
	double *v_rootpool = nullptr;					// Pointer to all particle velocities
	
	double Ek_root;									// Total kinetic energy
	
	double Ep_root;									// Total potential energy
	
	double Et_root;									// Total total energy.
	
	double *rho_root = nullptr;						// Pointer to all particle densities
	
	double *p_root   = nullptr;						// Pointer to all particle pressures
			
	
			
	/////
	//Coefficients for calculating density and forces
	/////
	
	/*Only calculated twice per process in whole code to save time
	 *Once when calling constructor and once after scaling mass
	 */
			
	double coeff_rho = m*4.0/(M_PI*h*h);            //Density
	
	double coeff_p = (-m/2.0)*(-30.0/(M_PI*h*h*h)); //Pressure force
	
	double coeff_v = -40.0*mu*m/(M_PI*h*h*h*h);     //Viscous force
	
	double hprime;                                  // 1.0/h
			
	/////		
	//Particle Initialization String
	/////
	
	string particles;								// String to store initial condition
	
	/////
	//MPI communicator variables
	/////
	
	MPI_Comm comm;									// MPI_COMM_WORLD
	
	int rank;										// Process rank
	
	int size;										// Process size
	
	//Variables to evenly distribute work amongst processes and to aid in scattering and gathering
	
	int start;										// Start of send array
	
	int finish;										// End of send array
	
	int lengthloc;									// No. of particles distributed to per process
	
	int *sendsizex;									// Array of sendcounts for x where sendsizex[i] corresponds to sendcount for ith rank
	
	int *stridex;									// Array of displ for x where stridex[i] corresponds to displ for ith rank
	
	int *sendsize;									// Array of sendcounts for p & rho where sendsize[i] corresponds to sendcount for ith rank
	
	int *stride;									// Array of displ for p & rho where stride[i] corresponds to displ for ith rank
	
	
};

#endif // SPH_H