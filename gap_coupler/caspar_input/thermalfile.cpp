# include "./inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
const char* thermal_file::general_kwds[] =
{
	"meshFile", "IR", "maxIters", "tolerance"
};
// ------------------------------------------------------------------------- //
const char* thermal_file::boundary_kwds[] =
{
	"type", "set"
};
// ------------------------------------------------------------------------- //
const char* thermal_file::dirichlet_kwds[] =
{
	"type", "set", "Tp"
};
// ------------------------------------------------------------------------- //
const char* thermal_file::neumann_kwds[] =
{
	"type", "set", "q"
};
// ------------------------------------------------------------------------- //
const char* thermal_file::mixed_kwds[] =
{
	"type", "set", "Tinf", "h"
};
// ------------------------------------------------------------------------- //
const char* thermal_file::constraints_kwds[] =
{
	"set", "x", "y", "z"
};
// ------------------------------------------------------------------------- //
int thermal_file::read(input_data& data)
{
	
	//if the file is not open do not check
	if(!isopen)
		return 0;

	string fname(filename);
	size_t pos = fname.find_first_not_of("./input/");
  fname = fname.substr(pos-1);

	inputlog << "Checking integrity of " << fname << " ... ";

	for(unsigned int s=0; s<sections.size(); s++)
	{
		if
		(
			sections[s]->name.compare("general") != 0 && 
			sections[s]->name.compare("materials") != 0 && 
			sections[s]->name.compare("boundary") != 0 && 
			sections[s]->name.compare("constraint") != 0 
		)
		{
			inputlog << logger::error() 
					 << "\n\nSection name " << sections[s]->name << " is invalid."
					 << "\nExpected 4 possible sections named:\n\n"
					 << "\t\"general\",\n"
					 << "\t\"materials\",\n"
					 << "\t\"boundary\",\n"
					 << "\t\"constraints\",\n"
					 << endl << "in " << filename << " input file."
					 << endl << endl;
		}	
	}

	// ------------- check if the sections are valid --------------- //
	
	int sz = 0;
	for(unsigned int i=0; i<sections.size(); i++)
	{
		if(sections[i]->name.compare("general") == 0)
		{
			sz = sizeof(general_kwds)/sizeof(general_kwds[0]);
			sections[i]->check_keywords(general_kwds, sz, filename);
		}
		else if(sections[i]->name.compare("boundary") == 0)
		{
			sz = sizeof(boundary_kwds)/sizeof(boundary_kwds[0]);
			sections[i]->check_keywords(boundary_kwds, sz, filename);
			// further check on the BC type
			istringstream iss(sections[i]->values[sections[i]->find_keyword("type")]);
			string bctype;
			iss >> bctype;
			if(bctype.compare("dirichlet") == 0)
			{
				sz = sizeof(dirichlet_kwds)/sizeof(dirichlet_kwds[0]);
				sections[i]->check_keywords(dirichlet_kwds, sz, filename);
			}
			else if(bctype.compare("mixed") == 0)
			{
				sz = sizeof(mixed_kwds)/sizeof(mixed_kwds[0]);
				sections[i]->check_keywords(mixed_kwds, sz, filename);
			}
			else if(bctype.compare("neumann") == 0)
			{
				string setname;
				istringstream(sections[i]->values[sections[i]->find_keyword("set")]) >> setname;
				sz = sizeof(neumann_kwds)/sizeof(neumann_kwds[0]);
				sections[i]->check_keywords(neumann_kwds, sz, filename);
			}
			else
			{
				inputlog << logger::error()
						 << "\n\nIn file " << filename << " (section " << sections[i]->name << ")\n" 
						 << "keyword \"type\" must have value [dirichlet/mixed/neumann], please check."
						 << endl << endl;
			}
		}
		else if(sections[i]->name.compare("constraint") == 0)
		{
			sz = sizeof(constraints_kwds)/sizeof(constraints_kwds[0]);
			sections[i]->check_keywords(constraints_kwds, sz, filename);
		}
		
	}

	inputlog << "done!" << endl;

	// ------------------ copy the geometric information --------------------- //

	// define references to limit the code repetitions
	string* _mesh = 0;
	bool* _IR = 0;
	double* _tol = 0;
	int* _maxit = 0;
	vector<pair<string,string>>* _materials;
	vector<constraint>* _constraints;
	vector<dirichlet>* _dirichlet_bc;
	vector<mixed>* _mixed_bc;
	vector<neumann>* _neumann_bc;
	vector<string>* _calc_deflect_on;

	
	if(string(filename).compare("./input/thermal_piston.txt") == 0)
	{
		_mesh = &data.thermal.piston.meshfile;
		_IR = &data.thermal.piston.IR;
		_tol = &data.thermal.piston.solver.tol;
		_maxit = &data.thermal.piston.solver.maxit;
		_constraints = &data.thermal.piston.constraints;
		_dirichlet_bc = &data.thermal.piston.dirichlet_bc;
		_mixed_bc = &data.thermal.piston.mixed_bc;
		_neumann_bc = &data.thermal.piston.neumann_bc;
		_materials = &data.thermal.piston.materials;

	}
	else if(string(filename).compare("./input/thermal_slipper.txt") == 0)
	{
		_mesh = &data.thermal.slipper.meshfile;
		_IR = &data.thermal.slipper.IR;
		_tol = &data.thermal.slipper.solver.tol;
		_maxit = &data.thermal.slipper.solver.maxit;
		_constraints = &data.thermal.slipper.constraints;
		_dirichlet_bc = &data.thermal.slipper.dirichlet_bc;
		_mixed_bc = &data.thermal.slipper.mixed_bc;
		_neumann_bc = &data.thermal.slipper.neumann_bc;
		_materials = &data.thermal.slipper.materials;
	}
	else if(string(filename).compare("./input/thermal_swashplate.txt") == 0)
	{
		_mesh = &data.thermal.swashplate.meshfile;
		_IR = &data.thermal.swashplate.IR;
		_tol = &data.thermal.swashplate.solver.tol;
		_maxit = &data.thermal.swashplate.solver.maxit;
		_constraints = &data.thermal.swashplate.constraints;
		_dirichlet_bc = &data.thermal.swashplate.dirichlet_bc;
		_mixed_bc = &data.thermal.swashplate.mixed_bc;
		_neumann_bc = &data.thermal.swashplate.neumann_bc;
		_materials = &data.thermal.swashplate.materials;
	}
	else if(string(filename).compare("./input/thermal_block.txt") == 0)
	{
		_mesh = &data.thermal.block.meshfile;
		_IR = &data.thermal.block.IR;
		_tol = &data.thermal.block.solver.tol;
		_maxit = &data.thermal.block.solver.maxit;
		_constraints = &data.thermal.block.constraints;
		_dirichlet_bc = &data.thermal.block.dirichlet_bc;
		_mixed_bc = &data.thermal.block.mixed_bc;
		_neumann_bc = &data.thermal.block.neumann_bc;
		_materials = &data.thermal.block.materials;
		_calc_deflect_on = &data.thermal.block.calc_deflect_on;
	}
	else if(string(filename).compare("./input/thermal_valveplate.txt") == 0)
	{
		_mesh = &data.thermal.valveplate.meshfile;
		_IR = &data.thermal.valveplate.IR;
		_tol = &data.thermal.valveplate.solver.tol;
		_maxit = &data.thermal.valveplate.solver.maxit;
		_constraints = &data.thermal.valveplate.constraints;
		_dirichlet_bc = &data.thermal.valveplate.dirichlet_bc;
		_mixed_bc = &data.thermal.valveplate.mixed_bc;
		_neumann_bc = &data.thermal.valveplate.neumann_bc;
		_materials = &data.thermal.valveplate.materials;
		_calc_deflect_on = &data.thermal.valveplate.calc_deflect_on;
	}

	for(unsigned int i=0; i<sections.size(); i++)
	{
		section* s = sections[i];
		// general
		if(s->name.compare("general") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("meshFile") == 0)
					istringstream(s->values[j]) >> *_mesh;
				// this is used sometime for the CB/VP interface, there is no check for 
				// this keyword, since is not always used
				if(s->keywords[j].compare("calc_deflect_on") == 0)
				{
					_calc_deflect_on->resize(0);
					istringstream iss(s->values[j]);
					while(true)
					{
						string tmp;
						iss >> tmp;
						if(tmp.size() > 0)
						{
							string val;
							istringstream(tmp) >> val;
							_calc_deflect_on->push_back(val);
						}
						else
							break;
					}
				}
				if(s->keywords[j].compare("IR") == 0)
					istringstream(s->values[j]) >> *_IR;
				if(s->keywords[j].compare("maxIters") == 0)
					istringstream(s->values[j]) >> *_maxit;
				if(s->keywords[j].compare("tolerance") == 0)
					istringstream(s->values[j]) >> *_tol;
			}
		}
		if(s->name.compare("materials") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				string thisset = s->keywords[j];
				string thismaterial;
				istringstream(s->values[j]) >> thismaterial;
				pair<string, string> entry;
				entry.first = thisset;
				entry.second = thismaterial;
				_materials->push_back(entry);
			}
		}
		if(s->name.compare("boundary") == 0)	
		{
			string bctype;
			istringstream(s->values[s->find_keyword("type")]) >> bctype;
			if(bctype.compare("dirichlet") == 0)
			{
				dirichlet bc;
				for(unsigned int j=0; j<s->keywords.size(); j++)
				{
					if(s->keywords[j].compare("set") == 0)
						istringstream(s->values[j]) >> bc.setname;
					if(s->keywords[j].compare("Tp") == 0)
						istringstream(s->values[j]) >> bc.Tp;
				}
				// add to the bc list
				_dirichlet_bc->push_back(bc);
			}
			else if(bctype.compare("mixed") == 0)
			{
				mixed bc;
				for(unsigned int j=0; j<s->keywords.size(); j++)
				{
					if(s->keywords[j].compare("set") == 0)
						istringstream(s->values[j]) >> bc.setname;
					if(s->keywords[j].compare("Tinf") == 0)
						istringstream(s->values[j]) >> bc.Tinf;
					if(s->keywords[j].compare("h") == 0)
						istringstream(s->values[j]) >> bc.h;
				}
				// add to the bc list
				_mixed_bc->push_back(bc);
			}
			else if(bctype.compare("neumann") == 0)
			{
				neumann bc;
				for(unsigned int j=0; j<s->keywords.size(); j++)
				{
					if(s->keywords[j].compare("set") == 0)
						istringstream(s->values[j]) >> bc.setname;
					if(s->keywords[j].compare("q") == 0)
					{
						double q;
						istringstream(s->values[j]) >> q;
						bc.q.push_back(q);
					}
				}
				// set to zero the heat flux value if is the gap surface
				if(bc.setname.compare("gap") == 0)
					bc.q.resize(0);
				// add to the bc list
				_neumann_bc->push_back(bc);
			}
		}
		if(s->name.compare("constraint") == 0)	
		{
			constraint c;
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("set") == 0)
					istringstream(s->values[j]) >> c.setname;
				if(s->keywords[j].compare("x") == 0)
					istringstream(s->values[j]) >> c.x;
				if(s->keywords[j].compare("y") == 0)
					istringstream(s->values[j]) >> c.y;
				if(s->keywords[j].compare("z") == 0)
					istringstream(s->values[j]) >> c.z;
			}
			// add to the bc list
			_constraints->push_back(c);
		}
	}

	return 0;
}
// ------------------------------------------------------------------------- //