# include "inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
const char* oil_file::general_kwds[] =
{
	"oiltype", "oillambda", "oilC"
};
// ------------------------------------------------------------------------- //
const char* oil_file::constant_properties_kwds[] =
{
	"oildensity", "oilviscosity", "oilK"
};
// ------------------------------------------------------------------------- //
const char* oil_file::user_defined_kwds[] =
{
	"w", "nupf", "nup1", "nup2", "nuT1", "nuT2", "rho0", "rhop1", "rhop2", 
	"rhoT", "rhopT", "pmin", "pmax", "Tmin", "Tmax"
};
// ------------------------------------------------------------------------- //
const char* oil_file::user_defined2_kwds[] =
{
	"rho0", "rhop1", "rhop2", "rhoT", "rhopT", "rhop2T" ,
	"pmin", "pmax", "Tmin", "Tmax", "c_2","d_2","g_0","s_0"
};
// ------------------------------------------------------------------------- //
int oil_file::read(input_data& data)
{
	if(!is_open())
	{
		inputlog << logger::error() << "\nUnable to open " << file_name << endl;
	}

	inputlog << "Checking integrity of oil.txt input file ... ";

	// ----------------------- check the file integrity ---------------------- //
	if(sections.size() < 1)
	{
		inputlog << logger::error()
				 << "\n\nExpected at least the general section in in oil.txt input file."
				 << endl << endl;
	}

	if(sections[0]->name.compare("general") != 0)
	{
		inputlog << logger::error()
				 << "\n\nExpected first section named \"general\" in oil.txt "
				 << " input file"
				 << endl << endl;
	}

	int general_sz = sizeof(general_kwds)/sizeof(general_kwds[0]);
	sections[0]->check_keywords(general_kwds, general_sz, "oil.txt");

	for(unsigned int j=0; j<sections[0]->keywords.size(); j++)
	{
		if(sections[0]->keywords[j].compare("oiltype") == 0)
			istringstream(sections[0]->values[j]) >> data.oil.general.oiltype;
		if(sections[0]->keywords[j].compare("oillambda") == 0)
			istringstream(sections[0]->values[j]) >> data.oil.general.oillambda;
		if(sections[0]->keywords[j].compare("oilC") == 0)
			istringstream(sections[0]->values[j]) >> data.oil.general.oilC;
	}

	// user_defined oil
	if(data.oil.general.oiltype == 0)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("constant_properties") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(constant_properties_kwds)/sizeof(constant_properties_kwds[0]);
				sections[i]->check_keywords(constant_properties_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("oildensity") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.constant_properties.oildensity;
					if(sections[i]->keywords[j].compare("oilviscosity") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.constant_properties.oilviscosity;
					if(sections[i]->keywords[j].compare("oilK") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.constant_properties.oilK;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"constant_properties\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}

	// user_defined oil
	if(data.oil.general.oiltype == 1)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined_kwds)/sizeof(user_defined_kwds[0]);
				sections[i]->check_keywords(user_defined_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;
					if(sections[i]->keywords[j].compare("w") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.w;
					if(sections[i]->keywords[j].compare("nupf") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.nupf;
					if(sections[i]->keywords[j].compare("nup1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.nup1;
					if(sections[i]->keywords[j].compare("nup2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.nup2;
					if(sections[i]->keywords[j].compare("nuT1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.nuT1;
					if(sections[i]->keywords[j].compare("nuT2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.nuT2;
					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	// user_defined oil2
	if(data.oil.general.oiltype == 3)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}

	// red_oil
	if(data.oil.general.oiltype == 4)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	// HLP32_oil
	if(data.oil.general.oiltype == 2)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	// SAE10W
	if(data.oil.general.oiltype == 5)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	// milh5606_oil
	if(data.oil.general.oiltype == 6)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	// ExxonDTE10Excel32
	if(data.oil.general.oiltype == 7)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	// ISO46_oil
	if(data.oil.general.oiltype == 8)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}

	//Parker_skydrol_oil
		if(data.oil.general.oiltype == 9)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}

	//water_oil- 50
		if(data.oil.general.oiltype == 50)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}

	//MILPRF87257_oil- 10
		if(data.oil.general.oiltype == 10)
	{
		bool found = false;
		for(unsigned int i=0; i<sections.size(); i++)
		{
			if(sections[i]->name.compare("user_defined") == 0)
			{
				found = true;

				// check keywords
				int sz = sizeof(user_defined2_kwds)/sizeof(user_defined2_kwds[0]);
				sections[i]->check_keywords(user_defined2_kwds, sz, "oil.txt");

				for(unsigned int j=0; j<sections[i]->keywords.size(); j++)
				{
					if(sections[i]->keywords[j].compare("pmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmin;
					if(sections[i]->keywords[j].compare("pmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.pmax;
					if(sections[i]->keywords[j].compare("Tmin") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmin;
					if(sections[i]->keywords[j].compare("Tmax") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.Tmax;

					if(sections[i]->keywords[j].compare("rho0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rho0;
					if(sections[i]->keywords[j].compare("rhop1") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop1;
					if(sections[i]->keywords[j].compare("rhop2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2;
					if(sections[i]->keywords[j].compare("rhoT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhoT;
					if(sections[i]->keywords[j].compare("rhopT") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhopT;
					if(sections[i]->keywords[j].compare("rhop2T") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.rhop2T;

					if(sections[i]->keywords[j].compare("g_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.g_0;
					if(sections[i]->keywords[j].compare("s_0") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.s_0;
					if(sections[i]->keywords[j].compare("c_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.c_2;
					if(sections[i]->keywords[j].compare("d_2") == 0)
						istringstream(sections[i]->values[j]) >> data.oil.user_defined.d_2;
				}
			}
		}
		if(!found)
		{
			inputlog << logger::error()
				 << "\n\nExpected a section named \"user_defined\" in oil.txt "
				 << " input file"
				 << endl << endl;
		}
	}
	inputlog << "done!" << endl;

	return 0;
}