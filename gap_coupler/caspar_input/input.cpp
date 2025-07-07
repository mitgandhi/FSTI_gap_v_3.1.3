# include <iostream>
# include <fstream>
# include <string>
# include <sstream>

# include "input.h"

using namespace std;

//the global input log
logger inputlog("input_log.txt", true);

input::input(const bool p_module_only) : 
	// initialize thermal files name
	thermal_piston("./input/thermal_piston.txt"),
	thermal_slipper("./input/thermal_slipper.txt"),
	thermal_swashplate("./input/thermal_swashplate.txt"),
	thermal_block("./input/thermal_block.txt"),
	thermal_valveplate("./input/thermal_valveplate.txt")
{

	//some values need to be defaulted in the event they do not exist (or are not read)
	setdefaults();
	
	//these three files always need to be read
	geometry.read(data);
	operatingconditions.read(data);
	oil.read(data);

	//only read the pmodule file if in the pressure module
	if(p_module_only)
	{
		pmodule.read(data);
	}

	//only read the lubrication and materials file if in the gap
	if(!p_module_only)
	{
		lubricationmodule.read(data);
		materials.read(data);
	}
	
	//only read the options file if the interface is enabled
	if(data.lubrication_module.solve_piston == 1)
	{
		optionspiston.read(data);
	}
	if(data.lubrication_module.solve_block == 1)
	{
		optionsblock.read(data);
	}
	if(data.lubrication_module.solve_slipper == 1)
	{
		optionsslipper.read(data);
	}

	// thermal files
	if(!p_module_only)
	{
		thermal_piston.read(data);
		thermal_slipper.read(data);
		thermal_swashplate.read(data);
		thermal_block.read(data);
		thermal_valveplate.read(data);
	}

	// check thermal files and options
	if(!thermal_piston.is_open() && data.options_piston.general.HeatTransfer && data.lubrication_module.solve_piston == 1)
	{
		inputlog << logger::error()
				 << endl << endl 
				 << "Could not find thermal_piston.txt but the HeatTransfer option in\n"
				 << "options_piston.txt file is activated. Please check."
				 << endl << endl;
	}
	if(!thermal_slipper.is_open() && data.options_slipper.general.SlipperThermoElastic && data.lubrication_module.solve_slipper == 1)
	{
		inputlog << logger::error()
				 << endl << endl 
				 << "Could not find thermal_slipper.txt but the SlipperThermoElastic option in\n"
				 << "options_slipper.txt file is activated. Please check."
				 << endl << endl;
	}
	if(!thermal_swashplate.is_open() && data.options_slipper.general.SwashplateThermoElastic && data.lubrication_module.solve_slipper == 1)
	{
		inputlog << logger::error() 
				 << endl << endl 
				 << "Could not find thermal_swashplate.txt but the SwashplateThermoElastic option in\n"
				 << "options_slipper.txt file is activated. Please check."
				 << endl << endl;
	}
	if(!thermal_block.is_open() && data.options_block.general.Thermal_CB && data.lubrication_module.solve_block == 1)
	{
		inputlog << logger::error()
				 << endl << endl 
				 << "Could not find thermal_block.txt but the Thermal_CB option in\n"
				 << "options_block.txt file is activated. Please check."
				 << endl << endl;
	}
	if(!thermal_valveplate.is_open() && data.options_block.general.Thermal_VP && data.lubrication_module.solve_block == 1)
	{
		inputlog << logger::error()
				 << endl << endl 
				 << "Could not find thermal_valveplate.txt but the Thermal_VP option in\n"
				 << "options_block.txt file is activated. Please check."
				 << endl << endl;
	}

}

void input::setdefaults()
{
	//we need to set defaults for some of the lubrication module options
	//so the pressure module doesn't have to load the files
	data.lubrication_module.n_lubrication_revolutions = 0;
	data.lubrication_module.solve_block = 0;
	data.lubrication_module.solve_piston = 0;
	data.lubrication_module.solve_slipper = 0;

	data.options_piston.general.HeatTransfer = false;
	data.options_block.general.Thermal_CB = false;
	data.options_block.general.Thermal_VP = false;
	data.options_slipper.general.SlipperThermoElastic = false;
	data.options_slipper.general.SwashplateThermoElastic = false;

	data.options_piston.general.EHDTestRig = 0;
	data.options_piston.general.TriboTestRig = 0;

	
	data.thermal.block.calc_deflect_on.push_back("ALL");
	data.thermal.valveplate.calc_deflect_on.push_back("ALL");

	data.geometry.gamma = 0;
	data.geometry.offset_J = 0;
	data.geometry.offset_K = 0;

}
