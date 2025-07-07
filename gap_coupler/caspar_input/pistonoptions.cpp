# include "inputfile.h"
# include <iostream>
# include <map>

using namespace std;

// ------------------------------------------------------------------------- //
const char* optionspiston_file::general_kwds[] =
{
	"ReadpFile", "ReynoldsMultiGrid", "EnergyEquation", "HeatTransfer", 
	"PressureDeformation", "PressureDeformationOMP", "ThermalDeformation", 
	"IM_piston_path", "IM_bushing_path",
	"McrK", "McrB", "McrK_file", "McrB_file", "zeta", "Friction_less", "Inclined_piston"
};
// ------------------------------------------------------------------------- //
const char* optionspiston_file::position_kwds[] =
{
	"xA", "xB", "yA", "yB"
};
// ------------------------------------------------------------------------- //
const char* optionspiston_file::numeric_kwds[] =
{
	"epsilonK", "jmax", "kmax", "delta_v", "AlphaP", "AlphaDef", "AlphaMu", 
	"AlphaTh", "Rmin_p", "Rmin_h", "nmax", "Rmin_R", "Rmin_E", "Simalphastep", "Simalphaplot",
	"EmodK", "vK", "EmodB", "vB", "hmin", "Tmax", "AlphaDC", "AlphaCase", "cgv", "pgv", "wgv"
};
// ------------------------------------------------------------------------- //
const char* optionspiston_file::GS_kwds[] =
{
	"N", "M", "Q"
};
// ------------------------------------------------------------------------- //
const char* optionspiston_file::GMG_kwds[] =
{
	"nL", "MG_M", "MG_N", "Q", "VW", "MGInt", "v1", "v2"
};
// ------------------------------------------------------------------------- //
int optionspiston_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of options_piston.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 5)
	{
		inputlog << logger::error()
				 << "\n\nExpected 5 sections named:\n\n"
				 << "\t\"general\",\n"
				 << "\t\"numeric\",\n"
				 << "\t\"position\",\n"
				 << "\t\"fluidgrid_GS\",\n"
				 << "\t\"fluidgrid_GMG\",\n"
				 << endl << "in " << "options_piston.txt input file."
				 << endl << endl;
	}

	for(int i=0; i<5; i++)
	{
		if
		(
			sections[i]->name.compare("general") != 0 && 
			sections[i]->name.compare("numeric") != 0 && 
			sections[i]->name.compare("position") != 0 && 
			sections[i]->name.compare("fluidgrid_GS") != 0 &&
			sections[i]->name.compare("fluidgrid_GMG") != 0
		)
		{
			inputlog << logger::error()
					 << "\n\nSection name " << sections[i]->name << " is invalid."
					 << "\nExpected 5 sections named:\n\n"
					 << "\t\"general\",\n"
					 << "\t\"numeric\",\n"
					 << "\t\"position\",\n"
					 << "\t\"fluidgrid_GS\",\n"
					 << "\t\"fluidgrid_GMG\",\n"
					 << endl << "in " << "options_piston.txt input file."
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
			sections[i]->check_keywords(general_kwds, sz, "options_piston.txt");
		}
		else if(sections[i]->name.compare("numeric") == 0)
		{
			sz = sizeof(numeric_kwds)/sizeof(numeric_kwds[0]);
			sections[i]->check_keywords(numeric_kwds, sz, "options_piston.txt");
		}
		else if(sections[i]->name.compare("position") == 0)
		{
			sz = sizeof(position_kwds)/sizeof(position_kwds[0]);
			sections[i]->check_keywords(position_kwds, sz, "options_piston.txt");
		}
		else if(sections[i]->name.compare("fluidgridGS") == 0)
		{
			sz = sizeof(GS_kwds)/sizeof(GS_kwds[0]);
			sections[i]->check_keywords(GS_kwds, sz, "options_piston.txt");
		}
		else if(sections[i]->name.compare("fluidgridGMG") == 0)
		{
			sz = sizeof(GMG_kwds)/sizeof(GMG_kwds[0]);
			sections[i]->check_keywords(GMG_kwds, sz, "options_piston.txt");
		}
	}
			
	// ------------------ copy the information --------------------- //

	for(unsigned int i=0; i<sections.size(); i++)
	{
		section* s = sections[i];
		// general
		if(s->name.compare("general") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("ReadpFile") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.ReadpFile;
				if(s->keywords[j].compare("ReynoldsMultiGrid") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.ReynoldsMultiGrid;
				if(s->keywords[j].compare("EnergyEquation") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.EnergyEquation;
				if(s->keywords[j].compare("HeatTransfer") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.HeatTransfer;
				if(s->keywords[j].compare("PressureDeformation") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.PressureDeformation;
				if(s->keywords[j].compare("PressureDeformationOMP") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.PressureDeformationOMP;
				if(s->keywords[j].compare("ThermalDeformation") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.ThermalDeformation;
				if(s->keywords[j].compare("EHDTestRig") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.EHDTestRig;
				if(s->keywords[j].compare("TriboTestRig") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.TriboTestRig;
				if(s->keywords[j].compare("IM_piston_path") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.IM_piston_path;
				if(s->keywords[j].compare("IM_bushing_path") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.IM_bushing_path;
				if(s->keywords[j].compare("McrK") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.McrK;
				if(s->keywords[j].compare("McrB") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.McrB;
				if(s->keywords[j].compare("McrK_file") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.McrK_file;
				if(s->keywords[j].compare("McrB_file") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.McrB_file;
				if(s->keywords[j].compare("Friction_less") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.friction_less;
				if(s->keywords[j].compare("Inclined_piston") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.inclined_piston;
				if(s->keywords[j].compare("zeta") == 0)
					istringstream(s->values[j]) >> data.options_piston.general.zeta;
				
			}
		}
		// numeric
		if(s->name.compare("numeric") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("epsilonK") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.epsilonK;
				if(s->keywords[j].compare("jmax") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.jmax;
				if(s->keywords[j].compare("kmax") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.kmax;
				if(s->keywords[j].compare("delta_v") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.delta_v;
				if(s->keywords[j].compare("AlphaP") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.AlphaP;
				if(s->keywords[j].compare("AlphaDef") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.AlphaDef;
				if(s->keywords[j].compare("AlphaMu") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.AlphaMu;
				if(s->keywords[j].compare("AlphaTh") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.AlphaTh;
				if(s->keywords[j].compare("Rmin_p") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Rmin_p;
				if(s->keywords[j].compare("Rmin_h") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Rmin_h;
				if(s->keywords[j].compare("nmax") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.nmax;
				if(s->keywords[j].compare("Rmin_R") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Rmin_R;
				if(s->keywords[j].compare("Rmin_E") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Rmin_E;
				if(s->keywords[j].compare("Simalphastep") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Simalphastep;
				if(s->keywords[j].compare("Simalphaplot") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Simalphaplot;
				if(s->keywords[j].compare("hmin") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.hmin, data.options_piston.numeric.hmin *= 1e-6;
				if(s->keywords[j].compare("EmodK") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.EmodK;
				if(s->keywords[j].compare("vK") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.vK;
				if(s->keywords[j].compare("EmodB") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.EmodB;
				if(s->keywords[j].compare("vB") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.vB;
				if(s->keywords[j].compare("Tmax") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.Tmax;
				if(s->keywords[j].compare("AlphaDC") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.AlphaDC;
				if(s->keywords[j].compare("AlphaCase") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.AlphaCase;
				if(s->keywords[j].compare("cgv") == 0)
				{
					istringstream iss(s->values[j]);
					while(true)
					{
						string tmp;
						iss >> tmp;
						if(tmp.size() > 0)
						{
							int val;
							istringstream(tmp) >> val;
							if(val > 0)
								data.options_piston.numeric.cgv.push_back(val);
						}
						else
							break;
					}
				}
				if(s->keywords[j].compare("pgv") == 0)
				{
					istringstream iss(s->values[j]);
					while(true)
					{
						string tmp;
						iss >> tmp;
						if(tmp.size() > 0)
						{
							double val;
							istringstream(tmp) >> val;
							data.options_piston.numeric.pgv.push_back(val);
						}
						else
							break;
					}
				}
				if(s->keywords[j].compare("wgv") == 0)
				{
					istringstream iss(s->values[j]);
					while(true)
					{
						string tmp;
						iss >> tmp;
						if(tmp.size() > 0)
						{
							double val;
							istringstream(tmp) >> val;
							data.options_piston.numeric.wgv.push_back(val);
						}
						else
							break;
					}
				}
				/*if(s->keywords[j].compare("pobgv") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.pobgv;
				if(s->keywords[j].compare("stgv") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.stgv;
				if(s->keywords[j].compare("wgv") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.wgv;
				if(s->keywords[j].compare("ngv") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.ngv;
				if(s->keywords[j].compare("spgv") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.spgv;
				if(s->keywords[j].compare("stploc") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.stploc;
				if(s->keywords[j].compare("stpdep") == 0)
					istringstream(s->values[j]) >> data.options_piston.numeric.stpdep;*/
			}
		}
		// position
		if(s->name.compare("position") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("xA") == 0)
				{
					istringstream(s->values[j]) >> data.options_piston.position.xA;
					data.options_piston.position.xA *= 1e-6;
				}
				if(s->keywords[j].compare("xB") == 0)
				{
					istringstream(s->values[j]) >> data.options_piston.position.xB;
					data.options_piston.position.xB *= 1e-6;

				}
				if(s->keywords[j].compare("yA") == 0)
				{
					istringstream(s->values[j]) >> data.options_piston.position.yA;
					data.options_piston.position.yA *= 1e-6;
				}
				if(s->keywords[j].compare("yB") == 0)
				{
					istringstream(s->values[j]) >> data.options_piston.position.yB;
					data.options_piston.position.yB *= 1e-6;
				}
			}
		}
		// GS fluid grid
		if(s->name.compare("fluidgrid_GS") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("N") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.GS.N;
				if(s->keywords[j].compare("M") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.GS.M;
				if(s->keywords[j].compare("Q") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.GS.Q;
			}
		}
		// GMG fluid grid
		if(s->name.compare("fluidgrid_GMG") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("nL") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.MG.nL;
				if(s->keywords[j].compare("MG_M") == 0)
				{
					istringstream iss(s->values[j]);
					while(true)
					{
						string tmp;
						iss >> tmp;
						if(tmp.size() > 0)
						{
							int val;
							istringstream(tmp) >> val;
							data.options_piston.fluid_grid.MG.MG_M.push_back(val);
						}
						else
							break;
					}
				}
				if(s->keywords[j].compare("MG_N") == 0)
				{
					istringstream iss(s->values[j]);
					while(true)
					{
						string tmp;
						iss >> tmp;
						if(tmp.size() > 0)
						{
							int val;
							istringstream(tmp) >> val;
							data.options_piston.fluid_grid.MG.MG_N.push_back(val);
						}
						else
							break;
					}
				}
				if(s->keywords[j].compare("Q") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.MG.Q;
				if(s->keywords[j].compare("VW") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.MG.VW;
				if(s->keywords[j].compare("MGInt") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.MG.MGInt;
				if(s->keywords[j].compare("v1") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.MG.v1;
				if(s->keywords[j].compare("v2") == 0)
					istringstream(s->values[j]) >> data.options_piston.fluid_grid.MG.v2;
			}
		}
	}

	// check integrity of the the levels
	if((data.options_piston.fluid_grid.MG.MG_N.size() != data.options_piston.fluid_grid.MG.nL)&&(data.options_piston.fluid_grid.MG.nL!=0))
	{

		inputlog << logger::error()
				 << "\nIn file options_piston.txt, the volumes in fluid circumference are " 
				 << data.options_piston.fluid_grid.MG.MG_N.size()
				 << " which differs from nL = " << data.options_piston.fluid_grid.MG.nL
				 << endl << endl;
		
	}
	if((data.options_piston.fluid_grid.MG.MG_M.size() != data.options_piston.fluid_grid.MG.nL)&&(data.options_piston.fluid_grid.MG.nL!=0))
	{
		inputlog << logger::error()
				 << "\nIn file options_piston.txt, the volumes in fluid length are " 
				 << data.options_piston.fluid_grid.MG.MG_M.size()
				 << " which differs from nL = " << data.options_piston.fluid_grid.MG.nL
				 << endl << endl;
	}

	inputlog << "done!" << endl;

	return 0;
}
