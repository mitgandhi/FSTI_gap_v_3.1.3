# include "inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
const char* materials_file::material_kwds[] =
{
	"name", "E", "nu", "rho", "lambda", "alpha" 
};

int materials_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of materials.txt input file ... ";

	// ------------- check if all geometric fields are defined --------------- //
	
	int nmaterials = 0;
	for(unsigned int i=0; i<sections.size(); i++)
	{
		if(sections[i]->name.compare("material") == 0)
		{
			int sz = sizeof(material_kwds)/sizeof(material_kwds[0]);
			sections[i]->check_keywords(material_kwds, sz, "materials.txt");
			nmaterials++;
		}
	}

	// ------------------ copy the geometric information --------------------- //

	for(unsigned int i=0; i<sections.size(); i++)
	{
		section* s = sections[i];
		// general
		if(s->name.compare("material") == 0)	
		{
			material m;
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("name") == 0)
					istringstream(s->values[j]) >> m.name;
				if(s->keywords[j].compare("E") == 0)
					istringstream(s->values[j]) >> m.E;
				if(s->keywords[j].compare("nu") == 0)
					istringstream(s->values[j]) >> m.nu;
				if(s->keywords[j].compare("rho") == 0)
					istringstream(s->values[j]) >> m.rho;
				if(s->keywords[j].compare("lambda") == 0)
					istringstream(s->values[j]) >> m.lambda;
				if(s->keywords[j].compare("alpha") == 0)
					istringstream(s->values[j]) >> m.alpha;
			}
			data.materials.push_back(m);
		}
	}

	inputlog << "done!" << endl;

	return 0;
}