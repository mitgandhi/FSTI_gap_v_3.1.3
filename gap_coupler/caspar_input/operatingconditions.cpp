# include "inputfile.h"
# include <iostream>
# include <cmath>

using namespace std;

const double pi = 3.14159265358979323846;

// ------------------------------------------------------------------------- //
const char* operatingconditions_file::kwds[] =
{
	"mode", "npistons", "speed", "beta", "betamax", "pModuleFile",
	"HP", "LP", "pCase", "T_HP", "T_LP", "T_Leak"
};
// ------------------------------------------------------------------------- //
int operatingconditions_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of operatingconditions.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 1)
	{
		inputlog << logger::error()
				 << "\n\nExpected just one section named \"operatingconditions\" "
				 << " in operatingconditions.txt input file."
				 << endl << endl;
	}

	if(sections[0]->name.compare("operatingconditions") != 0)
	{
		inputlog << logger::error()
				 << "\n\nExpected section named \"operatingconditions\" in "
				 << " operatingconditions.txt input file."
				 << endl << endl;
	}

	// ------------- check if all geometric fields are defined --------------- //
	
	int sz = sizeof(kwds)/sizeof(kwds[0]);
	sections[0]->check_keywords(kwds, sz, "operatingconditions.txt");
			
	// ---------------- copy the geometric information ------------------- //

	section* s = sections[0];

	for(unsigned int j=0; j<s->keywords.size(); j++)
	{
		if(s->keywords[j].compare("mode") == 0)
			istringstream(s->values[j]) >> data.operating_conditions.mode;
		if(s->keywords[j].compare("npistons") == 0)
			istringstream(s->values[j]) >> data.operating_conditions.npistons;
		if(s->keywords[j].compare("speed") == 0)
		{
			istringstream(s->values[j]) >> data.operating_conditions.speed;
			data.operating_conditions.speed *= 2.0*pi/60.0; // [rad/s]
		}
		if(s->keywords[j].compare("beta") == 0)
		{
			istringstream(s->values[j]) >> data.operating_conditions.beta;
			data.operating_conditions.beta *= pi/180.0; // [rad]
		}
		if(s->keywords[j].compare("betamax") == 0)
		{
			istringstream(s->values[j]) >> data.operating_conditions.betamax;
			data.operating_conditions.betamax *= pi/180.0; // [rad]
		}
		if(s->keywords[j].compare("pModuleFile") == 0)
			istringstream(s->values[j]) >> data.operating_conditions.pModuleFile;
		if(s->keywords[j].compare("HP") == 0)
		{
			istringstream(s->values[j]) >> data.operating_conditions.HP;
			data.operating_conditions.HP *= 1e5; // [bar]
		}
		if(s->keywords[j].compare("LP") == 0)
		{
			istringstream(s->values[j]) >> data.operating_conditions.LP;
			data.operating_conditions.LP *= 1e5; // [bar]
		}
		if(s->keywords[j].compare("pCase") == 0)
		{
			istringstream(s->values[j]) >> data.operating_conditions.pCase;
			data.operating_conditions.pCase *= 1e5; // [bar]
		}
		if(s->keywords[j].compare("T_HP") == 0)
			istringstream(s->values[j]) >> data.operating_conditions.T_HP;
		if(s->keywords[j].compare("T_LP") == 0)
			istringstream(s->values[j]) >> data.operating_conditions.T_LP;
		if(s->keywords[j].compare("T_Leak") == 0)
			istringstream(s->values[j]) >> data.operating_conditions.T_Leak;
	}

	inputlog << "done!" << endl;

	return 0;
}
// ------------------------------------------------------------------------- //
