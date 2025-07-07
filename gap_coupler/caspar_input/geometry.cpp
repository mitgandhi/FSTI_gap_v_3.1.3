# include "inputfile.h"
# include <iostream>

using namespace std;

const double pi = 3.14159265358979323846;

// ------------------------------------------------------------------------- //
const char* geometry_file::kwds[] =
{
	"dB", "dK", "dZ", "lK", "lF", "le", "lZ0", "dDK", "lDK", "lKG", "lch", "rK_red",
	"lK_hs", "speedK", "doutG", "dinG", "dDG", "lDG", "aSock", "r_socket", "lSK", 
	"lSG", "lG", "mK", "vPocket", "mG", "d_gap_in", "d_ope_in", 
	"d_ope_out", "d_groove_in", "d_groove_out", "d_gap_out", "d_spherical", "dBa", 
	"lengthB", "lengthcanalB", "delta_z0", "ADC", "DC_mesh", "mB", "IMB", "Fblock"
};
// ------------------------------------------------------------------------- //
int geometry_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of geometry.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() != 1)
	{
		inputlog << logger::error()
				 << "\n\nExpected just one section named \"geometry\" in geometry.txt "
				 << "input file."
				 << endl << endl;
	}

	if(sections[0]->name.compare("geometry") != 0)
	{
		inputlog << logger::error()
				 << "\n\nExpected section named \"geometry\" in geometry.txt "
				 << " input file"
				 << endl << endl;
	}

	// ------------- check if all geometric fields are defined --------------- //
	
	int sz = sizeof(kwds)/sizeof(kwds[0]);
	sections[0]->check_keywords(kwds, sz, "geometry.txt");
			
	// ---------------- copy the geometric information ------------------- //

	section* s = sections[0];

	for(unsigned int j=0; j<s->keywords.size(); j++)
	{
		if(s->keywords[j].compare("dB") == 0)
			istringstream(s->values[j]) >> data.geometry.dB, data.geometry.dB *= 1e-3;
		if(s->keywords[j].compare("dK") == 0)
			istringstream(s->values[j]) >> data.geometry.dK, data.geometry.dK *= 1e-3;
		if(s->keywords[j].compare("dZ") == 0)
			istringstream(s->values[j]) >> data.geometry.dZ, data.geometry.dZ *= 1e-3;
		if(s->keywords[j].compare("lK") == 0)
			istringstream(s->values[j]) >> data.geometry.lK, data.geometry.lK *= 1e-3;
		if(s->keywords[j].compare("lF") == 0)
			istringstream(s->values[j]) >> data.geometry.lF, data.geometry.lF *= 1e-3;
		if(s->keywords[j].compare("le") == 0)
			istringstream(s->values[j]) >> data.geometry.le, data.geometry.le *= 1e-3;
		if(s->keywords[j].compare("lZ0") == 0)
			istringstream(s->values[j]) >> data.geometry.lZ0, data.geometry.lZ0 *= 1e-3;
		if(s->keywords[j].compare("dDK") == 0)
			istringstream(s->values[j]) >> data.geometry.dDK, data.geometry.dDK *= 1e-3;
		if(s->keywords[j].compare("lDK") == 0)
			istringstream(s->values[j]) >> data.geometry.lDK, data.geometry.lDK *= 1e-3;
		if(s->keywords[j].compare("lKG") == 0)
			istringstream(s->values[j]) >> data.geometry.lKG, data.geometry.lKG *= 1e-3;
		if(s->keywords[j].compare("lch") == 0)
			istringstream(s->values[j]) >> data.geometry.lch, data.geometry.lch *= 1e-3;
		if(s->keywords[j].compare("rK_red") == 0)
			istringstream(s->values[j]) >> data.geometry.rK_red, data.geometry.rK_red *= 1e-6;
		if(s->keywords[j].compare("lK_hs") == 0)
			istringstream(s->values[j]) >> data.geometry.lK_hs, data.geometry.lK_hs *= 1e-3;
		if(s->keywords[j].compare("speedK") == 0)
			istringstream(s->values[j]) >> data.geometry.speedK;
		if(s->keywords[j].compare("doutG") == 0)
			istringstream(s->values[j]) >> data.geometry.doutG, data.geometry.doutG *= 1e-3;
		if(s->keywords[j].compare("dinG") == 0)
			istringstream(s->values[j]) >> data.geometry.dinG, data.geometry.dinG *= 1e-3;
		if(s->keywords[j].compare("dDG") == 0)
			istringstream(s->values[j]) >> data.geometry.dDG, data.geometry.dDG *= 1e-3;
		if(s->keywords[j].compare("lDG") == 0)
			istringstream(s->values[j]) >> data.geometry.lDG, data.geometry.lDG *= 1e-3;
		if(s->keywords[j].compare("aSock") == 0)
			istringstream(s->values[j]) >> data.geometry.aSock, data.geometry.aSock *= 1e-6;
		if(s->keywords[j].compare("r_socket") == 0)
			istringstream(s->values[j]) >> data.geometry.r_socket, data.geometry.r_socket *= 1e-3;
		if(s->keywords[j].compare("lSK") == 0)
			istringstream(s->values[j]) >> data.geometry.lSK, data.geometry.lSK *= 1e-3;
		if(s->keywords[j].compare("lSG") == 0)
			istringstream(s->values[j]) >> data.geometry.lSG, data.geometry.lSG *= 1e-3;
		if(s->keywords[j].compare("lG") == 0)
			istringstream(s->values[j]) >> data.geometry.lG, data.geometry.lG *= 1e-3;
		if(s->keywords[j].compare("mK") == 0)
			istringstream(s->values[j]) >> data.geometry.mK, data.geometry.mK *= 1e-3;
		if(s->keywords[j].compare("vPocket") == 0)
			istringstream(s->values[j]) >> data.geometry.vPocket, data.geometry.vPocket *= 1e-9;
		if(s->keywords[j].compare("mG") == 0)
			istringstream(s->values[j]) >> data.geometry.mG, data.geometry.mG *= 1e-3;
		if(s->keywords[j].compare("d_gap_in") == 0)
			istringstream(s->values[j]) >> data.geometry.d_gap_in, data.geometry.d_gap_in *= 1e-3;
		if(s->keywords[j].compare("d_ope_in") == 0)
			istringstream(s->values[j]) >> data.geometry.d_ope_in, data.geometry.d_ope_in *= 1e-3;
		if(s->keywords[j].compare("d_ope_out") == 0)
			istringstream(s->values[j]) >> data.geometry.d_ope_out, data.geometry.d_ope_out *= 1e-3;
		if(s->keywords[j].compare("d_groove_in") == 0)
			istringstream(s->values[j]) >> data.geometry.d_groove_in, data.geometry.d_groove_in *= 1e-3;
		if(s->keywords[j].compare("d_groove_out") == 0)
			istringstream(s->values[j]) >> data.geometry.d_groove_out, data.geometry.d_groove_out *= 1e-3;
		if(s->keywords[j].compare("d_gap_out") == 0)
			istringstream(s->values[j]) >> data.geometry.d_gap_out, data.geometry.d_gap_out *= 1e-3;
		if(s->keywords[j].compare("d_spherical") == 0)
			istringstream(s->values[j]) >> data.geometry.d_spherical, data.geometry.d_spherical *= 1e-3;
		if(s->keywords[j].compare("dBa") == 0)
			istringstream(s->values[j]) >> data.geometry.dBa, data.geometry.dBa *= 1e-3;
		if(s->keywords[j].compare("lengthB") == 0)
			istringstream(s->values[j]) >> data.geometry.lengthB, data.geometry.lengthB *= 1e-3;
		if(s->keywords[j].compare("lengthcanalB") == 0)
			istringstream(s->values[j]) >> data.geometry.lengthcanalB, data.geometry.lengthcanalB *= 1e-3;
		if(s->keywords[j].compare("delta_z0") == 0)
			istringstream(s->values[j]) >> data.geometry.delta_z0, data.geometry.delta_z0 *= 1e-3;
		if(s->keywords[j].compare("ADC") == 0)
			istringstream(s->values[j]) >> data.geometry.ADC, data.geometry.ADC *= 1e-6; // mm2 ->m2
		if(s->keywords[j].compare("DC_mesh") == 0)
			istringstream(s->values[j]) >> data.geometry.DC_mesh;
		if(s->keywords[j].compare("mB") == 0)
			istringstream(s->values[j]) >> data.geometry.mB;
		if(s->keywords[j].compare("IMB") == 0)
			istringstream(s->values[j]) >> data.geometry.IMB;
		if(s->keywords[j].compare("Fblock") == 0)
			istringstream(s->values[j]) >> data.geometry.Fblock;
		if(s->keywords[j].compare("gamma") == 0)
		{
			istringstream(s->values[j]) >> data.geometry.gamma;
			data.geometry.gamma *= (pi/180.0); // [deg] -> [rad]
		}
		if(s->keywords[j].compare("offset_J") == 0)
		{
			istringstream(s->values[j]) >> data.geometry.offset_J;
			data.geometry.offset_J *= 1e-3; // [mm] -> [m]
		}
		if(s->keywords[j].compare("offset_K") == 0)
		{
			istringstream(s->values[j]) >> data.geometry.offset_K;
			data.geometry.offset_K *= 1e-3; // [mm] -> [m]
		}
	}

	inputlog << "done!" << endl;

	return 0;
}
// ------------------------------------------------------------------------- //