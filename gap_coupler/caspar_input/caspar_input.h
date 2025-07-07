# ifndef __caspar_input__
# define __caspar_input__

# include "input.h"

class caspar_input : public input_data
{
public:

	struct Common
	{
		//Current time, phi (0-360+)
		double time;
		double phi_deg;
		double phi_rad;

		//For the current revolution (0-360), the current phi
		double phi_rev_deg;
		double phi_rev_rad;

		//what is the current revolution index
		int revolution_index;

		//The current DC pressure (Pa)
		double pDC;

		//The current HP pressure (Pa)
		double pHP;

		//The current LP pressure (Pa)
		double pLP;

		//how long is one revolution (s)
		double rev_period;

		//how many steps are in one revolution (int)
		int rev_steps;

		//how long is each timestep (s)
		double timestep;

		//how long is each timestep (deg)
		double phistep_deg;
		
		//This is the 'tolerance' in deg that will be used to "round" the phi angle at certain places
		double phi_deg_tol;

		//Flags that when true indicate 'special' timesteps
		bool first_step;
		bool last_step;
		bool first_rev_step;
		bool last_rev_step;

		//True if this is a resumed simulation
		bool resumed;

	} common;

	struct Pfile
	{
		std::vector<double> time;
		std::vector<double> pDC;
		std::vector<double> pHP;
		std::vector<double> pLP;
	} pfile;
	
	//the parent input class, if needed by methods inside the slipper
	const input * parent_input;

	//methods
	caspar_input(const input & I);
	void readpfile();

};

#endif
