# include "inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
const char* lubricationmodule_file::kwds[] =
{
	"solve_piston", "solve_block", "solve_slipper", "n_lubrication_revolutions","WobblePump", "Temperature_convergence"
};
// ------------------------------------------------------------------------- //
int lubricationmodule_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of lubrication_module.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 1)
	{
		inputlog << logger::error()
				 << "\n\nExpected just one section named \"lubrication_module\" "
				 << " in lubrication_module.txt input file."
				 << endl << endl;
	}

	if(sections[0]->name.compare("lubrication_module") != 0)
	{
		inputlog << logger::error()
				 << "\n\nExpected section named \"lubrication_module\" in "
				 << " lubrication_module.txt input file."
				 << endl << endl;
	}

	// ------------- check if all geometric fields are defined --------------- //
	
	int sz = sizeof(kwds)/sizeof(kwds[0]);
	sections[0]->check_keywords(kwds, sz, "lubrication_module.txt");
			
	// ---------------- copy the input information ------------------- //

	section* s = sections[0];

	for(unsigned int j=0; j<s->keywords.size(); j++)
	{
		if(s->keywords[j].compare("solve_piston") == 0)
			istringstream(s->values[j]) >> data.lubrication_module.solve_piston;
		if(s->keywords[j].compare("solve_block") == 0)
			istringstream(s->values[j]) >> data.lubrication_module.solve_block;
		if(s->keywords[j].compare("solve_slipper") == 0)
			istringstream(s->values[j]) >> data.lubrication_module.solve_slipper;
		if (s->keywords[j].compare("n_lubrication_revolutions") == 0) {
			int temp;
			istringstream(s->values[j]) >> temp;
			data.lubrication_module.n_lubrication_revolutions = (temp < 0) ? 20 : temp;
		}
		if(s->keywords[j].compare("WobblePump") == 0)
			istringstream(s->values[j]) >> data.lubrication_module.wobblepump;
		if(s->keywords[j].compare("Temperature_convergence") == 0)
			istringstream(s->values[j]) >> data.lubrication_module.temperature_convergence;
	}

	inputlog << "done!" << endl;

	return 0;
}
// ------------------------------------------------------------------------- //
