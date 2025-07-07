# include "inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
const char* pmodule_file::kwds[] =
{
	"AreaFile", "nrevolutions", "Vdead", "alphaD_LP", "alphaD_HP", "V_LP", 
	"V_HP", "AD_LP", "AD_HP", "P1", "P2", "leakageoption"
};
// ------------------------------------------------------------------------- //
int pmodule_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of p_module.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 1)
	{
		inputlog << logger::error()
				 << "\n\nExpected just one section named \"p_module\" "
				 << " in p_module.txt input file."
				 << endl << endl;
	}

	if(sections[0]->name.compare("p_module") != 0)
	{
		inputlog << logger::error()
				 << "\n\nExpected section named \"p_module\" in "
				 << " p_module.txt input file."
				 << endl << endl;
	}

	int sz = sizeof(kwds)/sizeof(kwds[0]);
	sections[0]->check_keywords(kwds, sz, "p_module.txt");
			
	// ---------------- copy the input information ------------------- //

	section* s = sections[0];

	for(unsigned int j=0; j<s->keywords.size(); j++)
	{
		if(s->keywords[j].compare("AreaFile") == 0)
			istringstream(s->values[j]) >> data.p_module.AreaFile;
		if(s->keywords[j].compare("nrevolutions") == 0)
			istringstream(s->values[j]) >> data.p_module.nrevolutions;
		if(s->keywords[j].compare("Vdead") == 0)
			istringstream(s->values[j]) >> data.p_module.Vdead;
		if(s->keywords[j].compare("alphaD_LP") == 0)
			istringstream(s->values[j]) >> data.p_module.alphaD_LP;
		if(s->keywords[j].compare("alphaD_HP") == 0)
			istringstream(s->values[j]) >> data.p_module.alphaD_HP;
		if(s->keywords[j].compare("V_LP") == 0)
			istringstream(s->values[j]) >> data.p_module.V_LP;
		if(s->keywords[j].compare("V_HP") == 0)
			istringstream(s->values[j]) >> data.p_module.V_HP;
		if(s->keywords[j].compare("AD_LP") == 0)
			istringstream(s->values[j]) >> data.p_module.AD_LP;
		if(s->keywords[j].compare("AD_HP") == 0)
			istringstream(s->values[j]) >> data.p_module.AD_HP;
		if(s->keywords[j].compare("P1") == 0)
		{
			istringstream(s->values[j]) >> data.p_module.P1;
			data.p_module.P1 *= 1e5;	// bar -> Pa
		}
		if(s->keywords[j].compare("P2") == 0)
		{
			istringstream(s->values[j]) >> data.p_module.P2;
			data.p_module.P2 *= 1e5;	// bar -> Pa
		}
		if(s->keywords[j].compare("leakageoption") == 0)
			istringstream(s->values[j]) >> data.p_module.leakageoption;
		if(s->keywords[j].compare("Q_Leak") == 0)
		{
			istringstream(s->values[j]) >> data.p_module.Q_Leak;
			data.p_module.Q_Leak /= 60000;  // l/min -> m3/s
		}
		if(s->keywords[j].compare("Auto_Area") == 0)      
			istringstream(s->values[j]) >> data.p_module.Auto_Area;
		if(s->keywords[j].compare("HPvp_area") == 0)      
			istringstream(s->values[j]) >> data.p_module.HPvp_area;
		if(s->keywords[j].compare("AD_HPtoLine") == 0)    
			istringstream(s->values[j]) >> data.p_module.AD_HPtoLine;
		if(s->keywords[j].compare("alphaD_HPtoLine") == 0)
			istringstream(s->values[j]) >> data.p_module.alphaD_HPtoLine;
		if(s->keywords[j].compare("V_line") == 0)         
			istringstream(s->values[j]) >> data.p_module.V_line;
		if(s->keywords[j].compare("phi0") == 0)           
			istringstream(s->values[j]) >> data.p_module.phi0;
		if(s->keywords[j].compare("phi2") == 0)           
			istringstream(s->values[j]) >> data.p_module.phi2;
		if(s->keywords[j].compare("phi3") == 0)          
			istringstream(s->values[j]) >> data.p_module.phi3;
		if(s->keywords[j].compare("phi4") == 0)          
			istringstream(s->values[j]) >> data.p_module.phi4;
		if(s->keywords[j].compare("vpv") == 0)            
			istringstream(s->values[j]) >> data.p_module.vpv;
		if(s->keywords[j].compare("vdv") == 0)            
			istringstream(s->values[j]) >> data.p_module.vdv;
		if(s->keywords[j].compare("Momentum") == 0)      
			istringstream(s->values[j]) >> data.p_module.Momentum;
		if(s->keywords[j].compare("IntegralHPfile") == 0)
			istringstream(s->values[j]) >> data.p_module.IntegralHPfile;
		if(s->keywords[j].compare("IntegralLPfile") == 0)
			istringstream(s->values[j]) >> data.p_module.IntegralLPfile;
        if(s->keywords[j].compare("FV") == 0)             
			istringstream(s->values[j]) >> data.p_module.FV;
		if(s->keywords[j].compare("FVAreaFile") == 0)
			istringstream(s->values[j]) >> data.p_module.FVAreaFile;
		if(s->keywords[j].compare("Air") == 0)            
			istringstream(s->values[j]) >> data.p_module.Air;
		if(s->keywords[j].compare("AirAreaFile") == 0)
			istringstream(s->values[j]) >> data.p_module.AirAreaFile;
		if(s->keywords[j].compare("Solver") == 0)         
			istringstream(s->values[j]) >> data.p_module.Solver;
	}

	inputlog << "done!" << endl;

	return 0;
}
// ------------------------------------------------------------------------- //
