# ifndef __input__
# define __input__

# include "./inputdata.h"

class input
{

	// operating conditions file
	operatingconditions_file operatingconditions;
	// operating conditions file
	lubricationmodule_file lubricationmodule;
	// geometry input file
	geometry_file geometry;
	// oil file
	oil_file oil;
	// pressure module file
	pmodule_file pmodule;
	// piston options file
	optionspiston_file optionspiston;
	// block options file
	optionsblock_file optionsblock;
	// block options file
	optionsslipper_file optionsslipper;
	// materials file
	materials_file materials;
	// thermal file for piston 
	thermal_file thermal_piston;
	// thermal file for slipper
	thermal_file thermal_slipper;
	// thermal file for swashplate
	thermal_file thermal_swashplate;
	// thermal file for cylinder block
	thermal_file thermal_block;
	// thermal file for valveplate
	thermal_file thermal_valveplate;

	//this method can be used to set input defaults throughout the input class
	//in the event they will not be read by this class for whatever reason
	void setdefaults();
	
public:

	// input data structure
	input_data data;

	// main constructor
	// p_module_only instructs the inputs class to only load:
	//		geometry.txt
	//		oil.txt
	//		operatingconditions.txt
	//		p_module.txt

	input(const bool p_module_only = false);

};

# endif