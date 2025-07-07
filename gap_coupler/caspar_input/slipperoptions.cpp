# include "inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
const char* optionsslipper_file::general_kwds[] =
{
	"reynolds_mu", "flowtype", "alphaD_KG", "EnableSlipperMacro", "SlipperMacroFile", "EnableRoughness", 
	"RoughnessRq", "CalcEnergy", "SlipperPressureDeformation", "SlipperThermoElastic",
	"SwashplatePressureDeformation", "SwashplateThermoElastic", "EHDsqueeze", "Explicit", "HybridForceBalance",
	"preFSIforceBalance", "ComplexPicard", "DebugMode", "DenseMode", "IMpath"
};
// ------------------------------------------------------------------------- //
const char* optionsslipper_file::position_kwds[] = 
{
	"hG1", "hG2", "hG3"
};


// ------------------------------------------------------------------------- //
const char* optionsslipper_file::numeric_kwds[] = 
{
	"AlphaReynolds", "AlphaEnergy", "AlphapG", "AlphaTEHD",
	"AlphaContact", "contact_def_lim", "contact_p_lim", 
	"Newton1", "Newton2", "Newton3", "FSIresidTol", "phistep_deg"
};
// ------------------------------------------------------------------------- //
const char* optionsslipper_file::fluidgrid_kwds[] = 
{
	"N", "M", "Q", "SealingLand", "Groove1Location", "Groove1r", "Groove1dr", "Groove2Location", "Groove2r", "Groove2dr"
};
// ------------------------------------------------------------------------- //
const char* optionsslipper_file::friction_kwds[] = 
{
	"ball_joint_friction", "C_joint", "mu_coef_dyn", "mu_coef_static", "mu_coef_pressure", "mu_coef_speed"
};

const char* optionsslipper_file::springforce_kwds[]=
{
	"hmaxG", "SlipHD", "FSlipperHD", "Fslipperspring", "Fspringconstant","Spring_length"
};
// ------------------------------------------------------------------------- //
int optionsslipper_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}
	
	inputlog << "Checking integrity of optionsslipper.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 6)
	{
		inputlog << logger::error()
				 << "\n\nExpected 6 sections named:\n\n"
				 << "\t\"general\",\n"
				 << "\t\"numeric\",\n"
				 << "\t\"position\",\n"
				 << "\t\"fluidgrid\",\n"
				 << "\t\"friction\",\n"
				 << "\t\"springforce\",\n"
				 << endl << "in " << "optionsslipper.txt input file."
				 << endl << endl;
	}

	for(int i=0; i<6; i++)
	{
		if
		(
			sections[i]->name.compare("general") != 0 && 
			sections[i]->name.compare("numeric") != 0 && 
			sections[i]->name.compare("position") != 0 && 
			sections[i]->name.compare("fluidgrid") != 0 &&
			sections[i]->name.compare("friction") != 0 &&
			sections[i]->name.compare("springforce") != 0 
		)
		{
			inputlog << logger::error()
					 << "\n\nSection name " << sections[i]->name << " is invalid."
					 << "\nExpected 5 sections named:\n\n"
					 << "\t\"general\",\n"
					 << "\t\"numeric\",\n"
					 << "\t\"position\",\n"
					 << "\t\"fluidgrid\",\n"
					 << "\t\"friction\",\n"
					 << "\t\"springforce\",\n"
					 << endl << "in " << "optionsslipper.txt input file."
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
			sections[i]->check_keywords(general_kwds, sz, "optionsslipper.txt");
		}
		else if(sections[i]->name.compare("numeric") == 0)
		{
			sz = sizeof(numeric_kwds)/sizeof(numeric_kwds[0]);
			sections[i]->check_keywords(numeric_kwds, sz, "optionsslipper.txt");
		}
		else if(sections[i]->name.compare("position") == 0)
		{
			sz = sizeof(position_kwds)/sizeof(position_kwds[0]);
			sections[i]->check_keywords(position_kwds, sz, "optionsslipper.txt");
		}
		else if(sections[i]->name.compare("fluidgrid") == 0)
		{
			sz = sizeof(fluidgrid_kwds)/sizeof(fluidgrid_kwds[0]);
			sections[i]->check_keywords(fluidgrid_kwds, sz, "optionsslipper.txt");
		}
		else if(sections[i]->name.compare("friction") == 0)
		{
			sz = sizeof(friction_kwds)/sizeof(friction_kwds[0]);
			sections[i]->check_keywords(friction_kwds, sz, "optionsslipper.txt");
		}
		else if(sections[i]->name.compare("springforce") == 0)
		{
			sz = sizeof(springforce_kwds)/sizeof(springforce_kwds[0]);
			sections[i]->check_keywords(springforce_kwds, sz, "optionsslipper.txt");
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
				if(s->keywords[j].compare("reynolds_mu") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.reynolds_mu;
				if(s->keywords[j].compare("flowtype") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.flowtype;
				if(s->keywords[j].compare("alphaD_KG") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.alphaD_KG;
				if(s->keywords[j].compare("EnableSlipperMacro") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.EnableSlipperMacro;
				if(s->keywords[j].compare("SlipperMacroFile") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.SlipperMacroFile;
				if(s->keywords[j].compare("EnableRoughness") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.EnableRoughness;
				if(s->keywords[j].compare("RoughnessRq") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.RoughnessRq;
				if(s->keywords[j].compare("CalcEnergy") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.CalcEnergy;
				if(s->keywords[j].compare("SlipperPressureDeformation") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.SlipperPressureDeformation;
				if(s->keywords[j].compare("SlipperThermoElastic") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.SlipperThermoElastic;
				if(s->keywords[j].compare("SwashplateThermoElastic") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.SwashplateThermoElastic;
				if(s->keywords[j].compare("SwashplatePressureDeformation") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.SwashplatePressureDeformation;
				if(s->keywords[j].compare("EHDsqueeze") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.EHDsqueeze;
				if(s->keywords[j].compare("HybridForceBalance") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.HybridForceBalance;
				if(s->keywords[j].compare("preFSIforceBalance") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.preFSIforceBalance;
				if(s->keywords[j].compare("Explicit") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.Explicit;
				if(s->keywords[j].compare("ComplexPicard") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.ComplexPicard;
				if(s->keywords[j].compare("DebugMode") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.DebugMode;
				if(s->keywords[j].compare("DenseMode") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.DenseMode;
				if(s->keywords[j].compare("IMpath") == 0)
					istringstream(s->values[j]) >> data.options_slipper.general.IMpath;
			}
		}
		// numeric
		if(s->name.compare("numeric") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("AlphaReynolds") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.AlphaReynolds;
				if(s->keywords[j].compare("AlphaEnergy") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.AlphaEnergy;
				if(s->keywords[j].compare("AlphapG") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.AlphapG;
				if(s->keywords[j].compare("AlphaTEHD") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.AlphaTEHD;
				if(s->keywords[j].compare("AlphaContact") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.AlphaContact;
				if(s->keywords[j].compare("contact_def_lim") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.contact_def_lim, data.options_slipper.numeric.contact_def_lim *= 1e-6;
				if(s->keywords[j].compare("contact_p_lim") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.contact_p_lim, data.options_slipper.numeric.contact_p_lim *= 1e5;
				if(s->keywords[j].compare("Newton1") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.Newton1;
				if(s->keywords[j].compare("Newton2") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.Newton2;
				if(s->keywords[j].compare("Newton3") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.Newton3;
				if(s->keywords[j].compare("FSIresidTol") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.FSIresidTol;
				if(s->keywords[j].compare("phistep_deg") == 0)
					istringstream(s->values[j]) >> data.options_slipper.numeric.phistep_deg;
			}
		}
		// position
		if(s->name.compare("position") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("hG1") == 0)
				{
					istringstream(s->values[j]) >> data.options_slipper.position.hG1;
					data.options_slipper.position.hG1 *= 1e-6;
				}
				if(s->keywords[j].compare("hG2") == 0)
				{
					istringstream(s->values[j]) >> data.options_slipper.position.hG2;
					data.options_slipper.position.hG2 *= 1e-6;
				}
				if(s->keywords[j].compare("hG3") == 0)
				{
					istringstream(s->values[j]) >> data.options_slipper.position.hG3;
					data.options_slipper.position.hG3 *= 1e-6;
				}
			}
		}
		// fluid grid
		if(s->name.compare("fluidgrid") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("N") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.N;
				if(s->keywords[j].compare("M") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.M;
				if(s->keywords[j].compare("Q") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Q;
				if(s->keywords[j].compare("SealingLand") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.SealingLand;
				if(s->keywords[j].compare("Groove1Location") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Groove1Location;
				if(s->keywords[j].compare("Groove1r") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Groove1r, data.options_slipper.fluid_grid.Groove1r *= 1e-3;
				if(s->keywords[j].compare("Groove1dr") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Groove1dr, data.options_slipper.fluid_grid.Groove1dr *= 1e-3;
				if(s->keywords[j].compare("Groove2Location") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Groove2Location;
				if(s->keywords[j].compare("Groove2r") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Groove2r, data.options_slipper.fluid_grid.Groove2r *= 1e-3;
				if(s->keywords[j].compare("Groove2dr") == 0)
					istringstream(s->values[j]) >> data.options_slipper.fluid_grid.Groove2dr, data.options_slipper.fluid_grid.Groove2dr *= 1e-3;

			}
		}
		// friction
		if(s->name.compare("friction") == 0)	
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("ball_joint_friction") == 0)
					istringstream(s->values[j]) >> data.options_slipper.friction.ball_joint_friction;
				if(s->keywords[j].compare("C_joint") == 0)
					istringstream(s->values[j]) >> data.options_slipper.friction.C_joint;
				if(s->keywords[j].compare("mu_coef_dyn") == 0)
					istringstream(s->values[j]) >> data.options_slipper.friction.mu_coef_dyn;
				if(s->keywords[j].compare("mu_coef_static") == 0)
					istringstream(s->values[j]) >> data.options_slipper.friction.mu_coef_static;
				if(s->keywords[j].compare("mu_coef_pressure") == 0)
					istringstream(s->values[j]) >> data.options_slipper.friction.mu_coef_pressure;
				if(s->keywords[j].compare("mu_coef_speed") == 0)
					istringstream(s->values[j]) >> data.options_slipper.friction.mu_coef_speed;
			}
			/*
			cout << endl << "slipper_options.cpp";
			cout << endl << "mu_coef_speed";
			cout << endl << data.options_slipper.friction.mu_coef_speed;
			cout << endl << "mu_coef_dyn";
			cout << endl << data.options_slipper.friction.mu_coef_dyn;
			cin.get();
			*/
		}
		//Spring-force-controller
		if(s->name.compare("springforce")==0)
		{
			for(unsigned int j=0; j<s->keywords.size(); j++)
			{
				if(s->keywords[j].compare("hmaxG") == 0)
					istringstream(s->values[j]) >> data.options_slipper.springforce.hmaxG, data.options_slipper.springforce.hmaxG *= 1e-6;
				if(s->keywords[j].compare("Fslipperspring") == 0)
					istringstream(s->values[j]) >> data.options_slipper.springforce.Fslipper;
				if(s->keywords[j].compare("SlipHD") == 0)
						istringstream(s->values[j]) >> data.options_slipper.springforce.SlipHD;
				if(s->keywords[j].compare("FSlipperHD") == 0)
						istringstream(s->values[j]) >> data.options_slipper.springforce.HoldDownStiffness;
				if(s->keywords[j].compare("Fspringconstant") == 0)
					istringstream(s->values[j]) >> data.options_slipper.springforce.Fspringconstant;
				if(s->keywords[j].compare("Spring_length") == 0)// from BDC bottom dead centre.
					istringstream(s->values[j]) >> data.options_slipper.springforce.LspringBDC;

			}
		}
	}

	inputlog << "done!" << endl;

	return 0;
}