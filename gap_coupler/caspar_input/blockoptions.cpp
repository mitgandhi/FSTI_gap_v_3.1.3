# include "inputfile.h"
# include <iostream>

using namespace std;

const double pi = 3.14159265358979323846;

// ------------------------------------------------------------------------- //
const char* optionsblock_file::general_kwds[] =
{
	"step_angle", "EHD_CB", "IM_CB", "EHD_VP", "IM_VP", 
	"Thermal_CB", "Thermal_VP", "StartWithTH", "macro_CB", "macro_VP", 
	"EHD_loop_debug", "dense_vtk_out_light", "dense_vtk_out"
};
// ------------------------------------------------------------------------- //
const char* optionsblock_file::position_kwds[] = 
{
	"hB1", "hB2", "hB3"
};
// ------------------------------------------------------------------------- //
const char* optionsblock_file::numeric_kwds[] = 
{
	"hmin", "contact", "epsilonB", "delta_v", "FEM_maxIters", "FEM_tolerance",
	"relax_CB", "relax_VP", "q_min_limit", "q_max_limit", "use_fsi_fb", 
	"use_sqz_hs", "use_sqz_hd"
};
// ------------------------------------------------------------------------- //
const char* optionsblock_file::fluidgrid_kwds[] = 
{
	"N", "M", "Q", "stl_cb", "stl_vp"
};
// ------------------------------------------------------------------------- //
int optionsblock_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of optionsblock.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 4)
	{
		inputlog << logger::error()
				 <<  "\n\nExpected 4 sections named:\n\n"
				 << "\t\"general\",\n"
				 << "\t\"numeric\",\n"
				 << "\t\"position\",\n"
				 << "\t\"fluidgrid\",\n"
				 << endl << "in " << "optionsblock.txt input file."
				 << endl << endl;
	}

	for(int i=0; i<4; i++)
	{
		if
		(
			sections[i]->name.compare("general") != 0 && 
			sections[i]->name.compare("numeric") != 0 && 
			sections[i]->name.compare("position") != 0 && 
			sections[i]->name.compare("fluidgrid") != 0 
		)
		{
			inputlog << logger::error()
					 << "\n\nSection name " << sections[i]->name << " is invalid."
					 << "\nExpected 4 sections named:\n\n"
					 << "\t\"general\",\n"
					 << "\t\"numeric\",\n"
					 << "\t\"position\",\n"
					 << "\t\"fluidgrid\",\n"
					 << endl << "in " << "optionsblock.txt input file."
					 << endl << endl;
		}	
	}

	// ------------- check if all geometric fields are defined --------------- //
	
	int sz = 0;
	for(unsigned int i=0; i<sections.size(); i++)
	{
		if(sections[i]->name.compare("general") == 0)
		{
			sz = sizeof(general_kwds)/sizeof(general_kwds[0]);
			sections[i]->check_keywords(general_kwds, sz, "optionsblock.txt");
		}
		else if(sections[i]->name.compare("numeric") == 0)
		{
			sz = sizeof(numeric_kwds)/sizeof(numeric_kwds[0]);
			sections[i]->check_keywords(numeric_kwds, sz, "optionsblock.txt");
		}
		else if(sections[i]->name.compare("position") == 0)
		{
			sz = sizeof(position_kwds)/sizeof(position_kwds[0]);
			sections[i]->check_keywords(position_kwds, sz, "optionsblock.txt");
		}
		else if(sections[i]->name.compare("fluidgrid") == 0)
		{
			sz = sizeof(fluidgrid_kwds)/sizeof(fluidgrid_kwds[0]);
			sections[i]->check_keywords(fluidgrid_kwds, sz, "optionsblock.txt");
		}
	}

	// ------------------ copy the geometric information --------------------- //

	for(unsigned int i=0; i<sections.size(); i++)
	{
		section* s = sections[i];
		// general
		if(s->name.compare("general") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("step_angle") == 0)
				{
					istringstream(s->values[j]) >> data.options_block.general.step_angle;
					data.options_block.general.step_angle *= pi/180.0; // convert to rad
				}
				if(s->keywords[j].compare("EHD_CB") == 0)
					istringstream(s->values[j]) >> data.options_block.general.EHD_CB;
				if(s->keywords[j].compare("IM_CB") == 0)
					istringstream(s->values[j]) >> data.options_block.general.IM_CB;
				if(s->keywords[j].compare("EHD_VP") == 0)
					istringstream(s->values[j]) >> data.options_block.general.EHD_VP;
				if(s->keywords[j].compare("IM_VP") == 0)
					istringstream(s->values[j]) >> data.options_block.general.IM_VP;
				if(s->keywords[j].compare("Thermal_CB") == 0)
					istringstream(s->values[j]) >> data.options_block.general.Thermal_CB;
				if(s->keywords[j].compare("Thermal_VP") == 0)
					istringstream(s->values[j]) >> data.options_block.general.Thermal_VP;
				if(s->keywords[j].compare("StartWithTH") == 0)
					istringstream(s->values[j]) >> data.options_block.general.StartWithTH;
				if(s->keywords[j].compare("macro_CB") == 0)
					istringstream(s->values[j]) >> data.options_block.general.macro_CB;
				if(s->keywords[j].compare("macro_VP") == 0)
					istringstream(s->values[j]) >> data.options_block.general.macro_VP;
				if(s->keywords[j].compare("EHD_loop_debug") == 0)
					istringstream(s->values[j]) >> data.options_block.general.EHD_loop_debug;
				if(s->keywords[j].compare("dense_vtk_out_light") == 0)
					istringstream(s->values[j]) >> data.options_block.general.dense_vtk_out_light;
				if(s->keywords[j].compare("dense_vtk_out") == 0)
					istringstream(s->values[j]) >> data.options_block.general.dense_vtk_out;
			}
		}
		// numeric
		if(s->name.compare("numeric") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("hmin") == 0)
				{
					istringstream(s->values[j]) >> data.options_block.numeric.hmin;
					data.options_block.numeric.hmin *= 1e-6;
				}
				if(s->keywords[j].compare("q_min_limit") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.q_min_limit;
				if(s->keywords[j].compare("q_max_limit") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.q_max_limit;
				if(s->keywords[j].compare("contact") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.contact;
				if(s->keywords[j].compare("delta_v") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.delta_v;
				if(s->keywords[j].compare("epsilonB") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.epsilonB;
				if(s->keywords[j].compare("FEM_maxIters") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.FEM_maxIters;
				if(s->keywords[j].compare("FEM_tolerance") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.FEM_tolerance;
				if(s->keywords[j].compare("relax_CB") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.relax_CB;
				if(s->keywords[j].compare("relax_VP") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.relax_VP;
				if(s->keywords[j].compare("use_fsi_fb") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.use_fsi_fb;
				if(s->keywords[j].compare("use_sqz_hs") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.use_sqz_hs;
				if(s->keywords[j].compare("use_sqz_hd") == 0)
					istringstream(s->values[j]) >> data.options_block.numeric.use_sqz_hd;

			}
		}
		// position
		if(s->name.compare("position") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("hB1") == 0)
				{
					istringstream(s->values[j]) >> data.options_block.position.hB1;
					data.options_block.position.hB1 *= 1e-6;
				}
				if(s->keywords[j].compare("hB2") == 0)
				{
					istringstream(s->values[j]) >> data.options_block.position.hB2;
					data.options_block.position.hB2 *= 1e-6;
				}
				if(s->keywords[j].compare("hB3") == 0)
				{
					istringstream(s->values[j]) >> data.options_block.position.hB3;
					data.options_block.position.hB3 *= 1e-6;
				}
			}
		}
		// fluid grid
		if(s->name.compare("fluidgrid") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("N") == 0)
					istringstream(s->values[j]) >> data.options_block.fluid_grid.N;
				if(s->keywords[j].compare("M") == 0)
					istringstream(s->values[j]) >> data.options_block.fluid_grid.M;
				if(s->keywords[j].compare("Q") == 0)
					istringstream(s->values[j]) >> data.options_block.fluid_grid.Q;
				if(s->keywords[j].compare("stl_cb") == 0)
					istringstream(s->values[j]) >> data.options_block.fluid_grid.stl_cb;
				if(s->keywords[j].compare("stl_vp") == 0)
					istringstream(s->values[j]) >> data.options_block.fluid_grid.stl_vp;
			}
		}

	}

	inputlog << "done!" << endl;

	return 0;
}