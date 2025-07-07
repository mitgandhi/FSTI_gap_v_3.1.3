# ifndef __inputdata__
# define __inputdata__

# include <vector>
# include <sstream>
# include "./inputfile.h"

// ------------------------------------------------------------------------- //
struct material
{
	std::string name;
	double E;
	double nu;
	double rho;
	double lambda;
	double alpha;
};
// ------------------------------------------------------------------------- //
struct dirichlet
{
	std::string setname;
	double Tp;
};
// ------------------------------------------------------------------------- //
struct mixed
{
	std::string setname;
	double h;
	double Tinf;
};
// ------------------------------------------------------------------------- //
struct neumann
{
	std::string setname;
	std::vector<double> q;
};
// ------------------------------------------------------------------------- //
struct constraint 
{
	std::string setname;
	bool x;
	bool y;
	bool z;
};
// ------------------------------------------------------------------------- //
struct input_data
{

public:

	// --------------------------- pump geometry ----------------------------- //
	
	struct
	{
		// Pitch Diameter Cylinder Block [mm]
		double dB;
		// Piston Diameter - [mm]
		double dK;
		// Bushing Diameter - [mm]
		double dZ;
		// Piston Length - [mm]  
		double lK;
		// Bushing Length - [mm]
		double lF;
		// Bushing Beginning Position	- [mm]	
		double le;
		// Displacement Chamber Length with Piston at ODP - [mm]		
		double lZ0;
		// Diameter Orifice Piston Head - [mm]
		double dDK;
		// Length Orifice Piston Head - [mm]		
		double lDK;
		// Length of Piston Gap surface - [mm]	
		double lKG;
		// Length of Piston Chamfer - [mm]	
		double lch;
		// Reduction Radius for Spherical Piston - [microns]
		double rK_red;
		// Length Cylindrical Section of the Piston (only for half spherical piston) - [mm]
		double lK_hs;
		// Piston Relative Rotational Speed [-]
		double speedK;
		// Slipper Outer Diameter - [mm]		
		double doutG;
		// Slipper Inner Diameter - [mm]
		double dinG;
		// Slipper Orifice Diameter - [mm]
		double dDG;
		// Slipper Orifice Length - [mm]	
		double lDG;
		// Slipper Socket Area (projected normal to the swashplate) - [mm^2]	
		double aSock;
		// Slipper Socket Radius - [mm]	
		double r_socket;
		// Distance to the Center of Mass of the Piston/Slipper Assembly - [mm]
		double lSK;
		// Distance from Piston Head to Slipper Center of Mass - [mm]	
		double lSG;
		// Distance from Piston Head to Gap Area - [mm]
		double lG;
		// Mass of Piston/Slipper Assembly - [g]	
		double mK;
		// Volume of slipper pocket - [m^3]
		double vPocket;
		
		// Slipper Mass - [g]					
		double mG;

		// gap inner diameter
		double d_gap_in;
		// opening inner diameter
		double d_ope_in;
		// opening outer diameter
		double d_ope_out;
		// groove separing the sealing land to the outer bearing, inner radius
		// (use 0 if not present)
		double d_groove_in;
		// groove separing the sealing land to the outer bearing, outer radius
		// (use 0 if not present)
		double d_groove_out;
		// gap outer diameter
		double d_gap_out;
		// diameter of the spherical gap (0 for flat gap)
		double d_spherical;
		// outer radius cylinder block
		double dBa;
		// Cylinder Block Total Length - [mm]		
		double lengthB;
		// Cylinder Block Canal Length (including the taper length if present) - [mm]	
		double lengthcanalB;
		// distance from the ideal shaft's crowing point 
		double delta_z0;
		// displacement chamber area (projection in the axial direction)
		double ADC;
		// file path where the surface mesh of the displacement chamber is stored
		std::string DC_mesh;
		// mass of the block [kg]
		double mB;
		// Moment of inertia of the block [kg*m^2] 
		// the reference is the spline center, 
		// the moment of inertia referes to any radial position (Ixx = Iyy)
		double IMB;
		// spring force cylinder block [N]
		double Fblock;
		// area where the spring force (pushing down the block) is acting on
		double spring_area;
		// Swashplate Cross Angle [deg]
		double gamma;
		// offset of the swashplate rotation axis from the pump origin in the y-axis. [m]
		double offset_J;
		// offset of the swashplate rotation axis from the pump origin in the z-axis. [m]
		double offset_K;

	} geometry;

	// ------------------------------ materials ------------------------------ //

	std::vector<material> materials;

	// --------------------------- oil properties ---------------------------- //

	struct
	{
		// general
		struct
		{
			int oiltype;					// oil type
			double oillambda;			// oil thermal conductivity
			double oilC;					// oil heat capacity
		} general;
		
		// oil with constant properties
		struct
		{
			double oilviscosity;	// oil dynamic viscosity
			double oildensity;		// oil density 
			double oilK;
		} constant_properties;
		
		// user defined oil
		struct
		{
			// range of validity
			double pmin;
			double pmax;
			double Tmin;
			double Tmax;
			
			// kinematic viscosity
			double w;			// weighting factor
			double nupf;	// first coeff for fp
			double nup1;	// second coeff for fp
			double nup2;  // first coeff for fT
			double nuT1;	// first coeff for fT
			double nuT2;	// second coeff for fT
			
			// density and bulk modulus
			double rho0;
			double rhop1;	// first coeff for density
			double rhop2;	// second coeff for density
			double rhoT;	// coeff for temperature
			double rhopT;	// coeff for temperature*pressure
			double rhop2T;	// coeff for temperature*pressure

			// Roealands Eq coefficients
			double s_0;	
			double g_0;  
			double c_2;	
			double d_2;	
		} user_defined;
	} oil;

	// ------------------------ operating conditions ------------------------- //

	struct 
	{
		// working mode 1:pumping 2:motoring
		unsigned int mode;
		// Number of Pistons
		unsigned int npistons;
		// Rotational speed [rpm]
		double speed;
		// Swash Plate Angle [deg]
		double beta;
		// Max Swash Plate Angle [deg]
		double betamax;

		// name of the pressure module file
		std::string pModuleFile;
		
		// High Pressure [bar]
		double HP;
		// Low pressure [bar]
		double LP;
		// Case pressure [bar]
		double pCase;

		// Temperature High Pressure [C]
		double T_HP;
		// Temperature Low Pressure [C]
		double T_LP;
		// Temperature Leakage [C]
		double T_Leak;

	} operating_conditions;

	// ------------------------ lubrication module ------------------------- //

	struct 
	{
		// enable the piston interface - 0: no, 1: yes
		unsigned int solve_piston;
		// enable the block interface - 0: no, 1: yes
		unsigned int solve_block;
		// enable the slipper interface - 0: no, 1: yes
		unsigned int solve_slipper;

		//how many revolutions should the lubrication module run
		int n_lubrication_revolutions;

		// enable the piston interface - 0: no, 1: yes
		unsigned int wobblepump;

		bool temperature_convergence;

	} lubrication_module;

	// ------------------------------ p module ------------------------------- //
	
	struct 
	{
		// area file
		std::string	AreaFile;
		// Number of Pressure Module revolutions []
		int nrevolutions;
		// Displacement chamber dead volume - [m^3]
		double Vdead;
		// Low pressure line flow coefficient - [-]
		double alphaD_LP;	
		// High pressure line flow coefficient - [-]
		double alphaD_HP;
		// Low pressure line volume - [m^3]
		double V_LP;
		// High pressure line volume - [m^3]
		double V_HP;
		// Low pressure cross section area - [m^2]
		double AD_LP;
		// High pressure cross section area - [m^2]
		double AD_HP;
		// Pumping Mode: upstream pressure low pressure throttle - [bar]
		// Motoring Mode: downstream pressure low pressure throttle - [bar]
		double P1;
		// Pumping mode: downstream pressure high pressure throttle - [bar]
		// Motoring mode: upstream pressure high pressure throttle - [bar]
		double P2;
		// Leakage Options
		// 0: use input file
		// 1: constant leakage value
		int leakageoption;
		// Constant Leakage Value - [l/min]
		double Q_Leak;
		//Automatic area control
		int Auto_Area;
		// Area Cover Delivery Volume [m^2](port plate side)		
		double HPvp_area;   
		// Area Cover Delivery Volume [m^2](HP port) - AD_HPtoLine_pump
		double AD_HPtoLine; 
		// Flow Coefficient Outlet port - alphaD_HPtoLine
		double alphaD_HPtoLine;
		// Volume of line [m^3]
		double V_line;
		// Beginning of ODC HP groove Angular position
		double phi0;
		// Ending of ODC HP groove Angular position
		double phi2;
		// Beginning of IDC LP groove Angular position
		double phi3;
		// Ending of IDC LP groove Angular position
		double phi4;
		//pre-compression volume [m3]
		double vpv;
		//de-compression volume [m3]
		double vdv;
		//Momentum equation
		//0: Disable
		//1: Momentum in HP port
		//2: Momentum in HP groove
		//3: Momentum in LP groove 
		//4: Momentum in HP and LP groove
		//5: Momentum in HP port and in HP groove
		//6: Momentum in HP port and in LP groove
		//7: Momentum in HP port and in HP and LP groove
		int Momentum;
		// Integral HP file
		std::string	IntegralHPfile;
		// Integral LP file
        std::string	IntegralLPfile;
		//FV-filter
		//0: Disable
		//1: PCFV filter
		//2: DCFV filter
		//3: PCFV and DCFV filter
		int FV;
        // FV area file
        std::string	FVAreaFile;
		//Air release port
		//0: Disable
		//1: Enable
		int Air;
        // Air area file
        std::string	AirAreaFile;
		//Solver
		//0: Non-stiff
		//1: Stiff
		int Solver;
	} p_module;


	// -------------------------- piston options ----------------------------- //

	struct
	{
	
		struct 
		{

			// Use simulated displacement chamber pressure
			bool ReadpFile; 
			// Solve Reynolds equation using MultiGrid
			bool ReynoldsMultiGrid;
			// Solve energy equation in the gap
			bool EnergyEquation;			
			// Solve temperature distribution in piston and cylinder
			bool HeatTransfer;
			// Solve elastic deformation due to pressure
			bool PressureDeformation;
			// Path to Piston IM
			std::string IM_piston_path;
			// Path to Bushing IM
			std::string IM_bushing_path;
			// Solve parallelized elastic deformation due to pressure
			int PressureDeformationOMP;
			// Solve elastic deformation due to temperature
			bool ThermalDeformation;
			// Simulate EHD test rig	
			bool EHDTestRig;
			// Simulate Tribo test rig	
			bool TriboTestRig;
			// Piston Macrogeometry
			int McrK;
			// Bushing Macrogeometry
			int McrB;
			// Piston Macrogeometry Path
			std::string McrK_file;
			// Bushing Macrogeometry Path
			std::string McrB_file;

			// to turn of friction 
			bool friction_less;
			// inclined piston option turn on
			bool inclined_piston;
			//piston inclination angle
			double zeta;
			
		} general;

		struct
		{
			// Piston initial positions: xA, yA, xB, yB - [m]
			// A: reference point closer to the displacement chamber
			// B: reference point closer to the swash plate
			double xA, yA, xB, yB;
		} position;

		struct
		{

			// -------------------------- force balance -------------------------- //

			// Newton Zero Finding Tolerance - [-]
			double epsilonK;
			// Newton Inner Loop Max Count - [-]
			unsigned int jmax;
			// Newton Outer Loop Max Count - [-]	
			unsigned int kmax;
			// Newton Delta Squeeze Velocity - [m/s]
			double delta_v;

			// ---------------- Outer Convergence Loop Parameters ---------------- //  
			
			// Damping Coefficient Pressure - [-]
			double AlphaP;
			// Damping Coefficient Deformation - [-]
			double AlphaDef;
			// Damping Coefficient Viscosity - [-]
			double AlphaMu;
			// Damping Coefficient Thermal - [-]
			double AlphaTh;
			// Minimum Fluid Residual - [-]
			double Rmin_p;
			// Minimum Structure Residual - [-]
			double Rmin_h;
			// Maximum Iterations Number [-]
			unsigned int nmax;

			// ------ Inner Convergence Loop Parameters ------//

			// Minimum Reynolds Equation Residual - [-]
			double Rmin_R;
			// Minimum Energy Equation Residual - [-]
			double Rmin_E;

			// ------ ODE Solver Parameters ------//

			// Simulation Maximum Step Size - [deg]
			double Simalphastep;
			// Simulation Maxumum Plot Step Size - [deg]
			double Simalphaplot;

			// ------ Antipenetration Material Properties ------ //

			// Piston elastic modulus - [Pa]
			double EmodK;
			// Piston Poisson ratio - [-]
			double vK;
			// Cylinder elastic modulus - [Pa]
			double EmodB;
			// Cylinder Poisson ratio - [-]
			double vB;

			// ------ Antipenetration Parameters ------ //

			// Minimum Gap Height - [microns]
			double hmin;

			// ------ Thermal Parameters ------ //

			// Estimate Maximum Gap Temperature - [°C]		
			double Tmax;
			//Convection Coefficient Displacement Chamber - [W/m2°C]
			double AlphaDC;
			//Convection Coefficient Case - [W/m2°C]
			double AlphaCase;

			// ------ Groove Profile before 03062015------ //

			// Groove on piston or bushing? 1: piston linear; 2: bushing linear; 3:piston constant; 4:bushing constant
			/*double pobgv;
			// Distance between the first groove and the starting point of the end of the piston closer to DC - [m]		
			double stgv;
			// Width of the groove - [m]
			double wgv;
			// number of the grooves - [N/A]
			int ngv;
			// spacing between two gooves - [m]
			double spgv;*/

			// ------ Groove Profile after 03062015------ //

			// Groove on piston or bushing? 1: piston linear; 2: bushing linear; 3:piston constant; 4:bushing constant
			std::vector<int> cgv;
			// Distance between the first groove and the starting point of the end of the piston closer to DC - [m]		
			std::vector<double> pgv;
			// Width of the groove - [m]
			std::vector<double> wgv;

			// ------ Step Profile ------ //

			// Distance from gap end (near Case) to step - [m]		
			double stploc;
			// depth of the step - [m]
			double stpdep;

		} numeric;

		struct
		{
			// Gauss-Seidel fluid grid definition
			struct
			{
				// Volumes in fluid film circumference
				unsigned int N;
				// Volumes in fluid film length
				unsigned int M;
				// Volumes in fluid film height
				unsigned int Q;
			} GS;

			// Multigrid fluid grid definition
			struct
			{
				// number of levels
				unsigned int nL;
				// Multigrid Mesh Definition [N M] Level 0 to nL 
				// (N volumes in fluid film circumference, M volumes in fluid film length)
				std::vector<unsigned int> MG_M;
				std::vector<unsigned int> MG_N;

				// Multigrid Mesh Definition [Q] All Levels (Q volumes in fluid film height)
				unsigned int Q;

				// Type of Cycle: 0 for V-Cycle; 1 for W-Cycle
				unsigned int VW;
				// Type of Prolongation: 0 for Standard; 1 for Bilinear
				unsigned int MGInt;
				// Number of GS Sweeps Down-Leg
				unsigned int v1;
				// Number of GS Sweeps Up-Leg
				unsigned int v2;

			} MG;

		} fluid_grid;		


	} options_piston;

	// -------------------------- slipper options ---------------------------- //

	struct
	{

		// ------------------------- general settings -------------------------- //
		struct
		{
			// Solve Reynolds equation using viscosity at (1) cell center or (2) cell faces
			// Also known as (1) average mu or (2) full mu
			int reynolds_mu;
			// Flow Type for Piston/Slipper
			// 0: turbulent flow through piston and slipper orifices
			// 1: laminar flow through piston and slipper orifices
			// 2: turbulent flow through piston orifice
			// 3: laminar flow through piston orifice						
			unsigned int flowtype;
			// Orifice Coefficient Piston/Slipper Assembly - [-]	
			double alphaD_KG;
			// Slipper macro geometry
			bool EnableSlipperMacro;
			// path to the slipper macro geometry file
			std::string SlipperMacroFile;
			// Consider slipper surface roughness
			bool EnableRoughness;
			//  roughness Rq [um]
			double RoughnessRq;
			// Solve energy equation in the gap
			bool CalcEnergy;
			// Solve elastic deformation due to pressure on the slipper
			bool SlipperPressureDeformation;
			// Consider the thermo elastic problem for the slipper
			bool SlipperThermoElastic;
			// Solve elastic deformation due to pressure on the swashplate
			bool SwashplatePressureDeformation;
			// Consider the thermo elastic problem for the swashplate
			bool SwashplateThermoElastic;
			// 0: Disable EHD squeeze, 1: Enable EHD squeeze
			bool EHDsqueeze;
			// 0: Use newton method. 1: Use hybrid force balance
			bool HybridForceBalance;
			// 0: Don't use a pre-FSI force balance, 1: Use a pre-FSI force balance
			bool preFSIforceBalance;
			// 0: Implicit ODE, 1: Explicit ODE
			bool Explicit;
			// This will enable the full picard iteration for all  calculations 
			// (Not normally recommended)
			int ComplexPicard;
			// Debug mode enables extra logging
			int DebugMode;
			// Dense mode will only solve for a single time step but with dense FSI results
			int DenseMode;
			// Absolute or relative path to the slipper and/or swashplate IM folder
			std::string IMpath;
		
		} general;

		struct
		{
			// Slipper initial positions: hG1, hG2, hG3 - [m]	
			double hG1, hG2, hG3;
		} position;

		struct
		{

			// ------------------------- Slipper Parameters ---------------------- //
			
			// over-relaxation coeff for Gauss-Seidel solution for Reynolds
			double AlphaReynolds;
			// over-relaxation coeff for Gauss-Seidel solution for Energy
			double AlphaEnergy;
			// under relaxation parameter for the poket pressure calculation
			double AlphapG;
			// under relaxation parameter for the TEHD calculation
			double AlphaTEHD;
			// under relaxation parameter for the TEHD calculation
			double AlphaContact;
			// under relaxation parameter for the TEHD calculation
			double contact_def_lim;
			// under relaxation parameter for the TEHD calculation
			double contact_p_lim;
			// parameters used to control the newton method
			double Newton1;
			double Newton2;
			double Newton3;
			// FSI residual tolerance
			double FSIresidTol;
			//Simulation step size
			double phistep_deg;
		
		} numeric;

		struct
		{

			// Volumes in circumferential direction
			int N;
			// Volumes in radial direction
			int M;
			// Volumes in height direction
			int Q;
			// Which ring of cells will be used to calculate slipper leakage
			int SealingLand;

			// The radial fluid grid cell id of a groove connected to pocket pressure
			int Groove1Location;
			// The center radius of Groove1
			double Groove1r;
			// The total dr of Groove1
			double Groove1dr;

			// The radial fluid grid cell id of a groove connected to case pressure
			int Groove2Location;
			// The center radius of Groove2
			double Groove2r;
			// The total dr of Groove2
			double Groove2dr;


		} fluid_grid;

		struct
		{
			// In the piston/slipper ball joint, calculate (0) no friction, (1) constant coefficients, (2) pressure-speed constant coefficients, (3) arctan, (4) Stribeck curve	
			unsigned int ball_joint_friction;

			// Friction reduction factor in the ball joint. 
			// Multiplied by F_friction calculated from any of the following methods. 
			double C_joint;

			//Option (1):  F_friction = mu * F_SK
			// Coefficient of static friction 
			//Literature value = 0.16, but stability issues are worsened when static and dynamic coeffs are not equal.
			double mu_coef_static;
			// Coefficient of dynamic friction 
			double mu_coef_dyn;


			//Option (2):  F_friction = (mu_p * p [bar] + mu_v * v [rps]) * F_SK
			// Coefficient of pressure-dependent friction
			double mu_coef_pressure;
			// Coefficient of tilting speed-dependent friction
			double mu_coef_speed;


			//Option (3):  F_friction = arctan(tilting speed) * F_SK


			//Option (4):  F_friction = fct(viscosity, pressure, tilting speed) * F_SK


		} friction;
		
		struct
		{
			// Maximal Height of the Slipper Holder - [microns]	
			// (only used for fixed clearance holder, positive spring holder: -1000) 	
			double hmaxG;
			// Slipper Holder Stiffness Coefficient [N/m] (if fixed clearance holder is chosen)	
			// Slipper Total Spring Force [N] (if positive slipper spring holder)	
			double Fslipper;
			// Slipper Spring Force for applying Hold down effect. 
			int SlipHD;
			// Slipper Hold down Stiffnes;
			double HoldDownStiffness;

			// spring constant for the spring effect from the piston
			double Fspringconstant;

			// spring length from BDC
			double LspringBDC;

		} springforce;

	} options_slipper;

	// --------------------------- block options ----------------------------- //

	struct
	{
		// ------------------------- general settings -------------------------- //
		struct
		{
			// simulation step angle
			double step_angle;
			// Activate/Deactivate the EHD on the cylinder block side
			bool EHD_CB;
			// specify the folder containing the influence matrices
			std::string IM_CB;
			// Activate/Deactivate the EHD on the valve plate side 
			bool EHD_VP;
			// specify the folder containing the influence matrices
			std::string IM_VP;
			// Activate/Deactivate the Thermal analysis on the cylinder block side
			bool Thermal_CB;
			// Activate/Deactivate the Thermal analysis on the valve plate side
			bool Thermal_VP;
			// Start the simulation with the thermal analysis
			bool StartWithTH;

			// enable macro geometry on block: (0 disabled 1,2,3,4... different macro types)
			int macro_CB;
			// enable macro geometry on valve plate: (0 disabled 1,2,3,4... different macro types)
			int macro_VP;
			// print residual vs iterations for each EHD loop
			bool EHD_loop_debug;
			// write a VTK file (light version) for each time step
			bool dense_vtk_out_light;
			// write a VTK file (full version) for each time step
			bool dense_vtk_out;
			
			//bool thermal_file_out; 

		} general;

		struct
		{
			// block position over the three control points
			// P1 is @ 0 deg, P2 is @ 120 deg, P3 is @ 240 deg
			double hB1, hB2, hB3;
		} position;

		// ------------------------- numeric settings -------------------------- //
		struct
		{
			
			// minimum film thickness [um]
			double hmin;
			// contact option: 
			// 0 -> do nothing
			// 1 -> integrate just positive squeeze in contact conditions
			// 2 -> use squeeze velocity correction
			int contact;

			// minimum heat flux
			double q_min_limit;
			// maximum heat flux
			double q_max_limit;

			// use the hydrostatic squeeze term
			bool use_sqz_hd;
			// use the hydrodynamic squeeze term
			bool use_sqz_hs;

			// ---------- force balance ---------- //

			// Newton Zero Finding Tolerance - [-]
			double epsilonB;
			// Newton Delta Squeeze Velocity - [m/s]
			double delta_v;
			// use FSI term when solving Reynolds equation in the force balance
			bool use_fsi_fb;
			// External forces calculation method: 0 simplified, 1 advanced
			int Fext_method;

			// ----------- finite element thermal solver settings ---------------- //

			// max number of iterations 
			unsigned int FEM_maxIters;
			// solver tolerance
			unsigned int FEM_tolerance;
			// relaxation factor used for the cylinder block thermal analysis
			double relax_CB;
			// relaxation factor used for the valve plate thermal analysis
			double relax_VP;


		} numeric;

		// ----------------------- fluid grid settings ------------------------- //
		struct
		{
			// Volumes in circumferential direction
			int N;
			// Volumes in radial direction
			int M;
			// Volumes in height direction
			int Q;
			// stl file used to define the cylinder block gap
			std::string stl_cb;
			// stl file used to define the valve plate gap
			std::string stl_vp;

		} fluid_grid;

	} options_block;

	// -------------------------- thermal options ---------------------------- //

	struct
	{
		// piston thermal boundaries
		struct
		{
			// thermal mesh file
			std::string meshfile;
			// inertia relief option
			bool IR;
			struct
			{
				// linear solver convergence tolerance
				double tol;
				// linear solver maximum iteration number
				int maxit;
			} solver;
			// material association: list of set -> associated material
			std::vector<std::pair<std::string, std::string>> materials;
			// boundaries and constraints
			std::vector<constraint> constraints;
			std::vector<dirichlet> dirichlet_bc;
			std::vector<mixed> mixed_bc;
			std::vector<neumann> neumann_bc;
		} piston;
		// slipper thermal boundaries
		struct
		{
			// thermal mesh file
			std::string meshfile;
			// inertia relief option
			bool IR;
			struct
			{
				// linear solver convergence tolerance
				double tol;
				// linear solver maximum iteration number
				int maxit;
			} solver;
			// material association: list of set -> associated material
			std::vector<std::pair<std::string, std::string>> materials;
			// boundaries and constraints
			std::vector<constraint> constraints;
			std::vector<dirichlet> dirichlet_bc;
			std::vector<mixed> mixed_bc;
			std::vector<neumann> neumann_bc;
		} slipper;
		// swashplate thermal boundaries
		struct
		{
			// thermal mesh file
			std::string meshfile;
			// inertia relief option
			bool IR;
			struct
			{
				// linear solver convergence tolerance
				double tol;
				// linear solver maximum iteration number
				int maxit;
			} solver;
			// material association: list of set -> associated material
			std::vector<std::pair<std::string, std::string>> materials;
			// boundaries and constraints
			std::vector<constraint> constraints;
			std::vector<dirichlet> dirichlet_bc;
			std::vector<mixed> mixed_bc;
			std::vector<neumann> neumann_bc;
		} swashplate;
		// cylinder block thermal boundaries
		struct
		{
			// thermal mesh file
			std::string meshfile;
			// inertia relief option
			bool IR;
			struct
			{
				// linear solver convergence tolerance
				double tol;
				// linear solver maximum iteration number
				int maxit;
			} solver;
			// material association: list of set -> associated material
			std::vector<std::pair<std::string, std::string>> materials;
			// boundaries and constraints
			std::vector<constraint> constraints;
			std::vector<dirichlet> dirichlet_bc;
			std::vector<mixed> mixed_bc;
			std::vector<neumann> neumann_bc;
			// volume set in the mesh used for the calculation of the
			// thermal deflection. Use ALL to run the analysis on the 
			// whole mesh
			std::vector<std::string> calc_deflect_on;
		} block;
		// valveplate thermal boundaries
		struct
		{
			// thermal mesh file
			std::string meshfile;
			// inertia relief option
			bool IR;
			struct
			{
				// linear solver convergence tolerance
				double tol;
				// linear solver maximum iteration number
				int maxit;
			} solver;
			// material association: list of set -> associated material
			std::vector<std::pair<std::string, std::string>> materials;
			// boundaries and constraints
			std::vector<constraint> constraints;
			std::vector<dirichlet> dirichlet_bc;
			std::vector<mixed> mixed_bc;
			std::vector<neumann> neumann_bc;
			// volume set in the mesh used for the calculation of the
			// thermal deflection. Use ALL to run the analysis on the 
			// whole mesh
			std::vector<std::string> calc_deflect_on;
		} valveplate;
	} thermal;

};
// ------------------------------------------------------------------------- //

# endif