# ifndef __inputfile__
# define __inputfile__

# include <vector>
# include <sstream>
# include <fstream>
# include "./section.h"
# include "./inputdata.h"

struct input_data;

// ------------------------------------------------------------------------- //
class input_file
{

protected:

	bool isopen;

	std::ifstream file;

	std::vector<section*> sections;
	
	input_file(const char* filename);
	~input_file();
	
	// read the file
	bool read(const char* filename);
	
	// check the file sections
	virtual int check_sections();

public:

	//what is the filename
	std::string file_name;

	// return file status
	bool is_open() const { return isopen; }

};
// ------------------------------------------------------------------------- //
class geometry_file : public input_file
{
	static const char* kwds[];

public:
	
	geometry_file() : input_file("./input/geometry.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class operatingconditions_file : public input_file
{
	static const char* kwds[];

public:
	
	operatingconditions_file() : input_file("./input/operatingconditions.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class lubricationmodule_file : public input_file
{
	static const char* kwds[];

public:
	
	lubricationmodule_file() : input_file("./input/lubrication_module.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class pmodule_file : public input_file
{
	static const char* kwds[];

public:
	
	pmodule_file() : input_file("./input/p_module.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class oil_file : public input_file
{
	static const char* general_kwds[];
	static const char* constant_properties_kwds[];
	static const char* user_defined_kwds[];
	static const char* user_defined2_kwds[];

public:
	
	oil_file() : input_file("./input/oil.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class optionspiston_file : public input_file
{
	
	static const char* general_kwds[];
	static const char* position_kwds[];
	static const char* numeric_kwds[];
	static const char* GS_kwds[];
	static const char* GMG_kwds[];
	
public:
	
	optionspiston_file() : input_file("./input/options_piston.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class optionsblock_file : public input_file
{
	
	static const char* general_kwds[];
	static const char* position_kwds[];
	static const char* numeric_kwds[];
	static const char* fluidgrid_kwds[];
	
	
public:
	
	optionsblock_file() : input_file("./input/options_block.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class optionsslipper_file : public input_file
{
	
	static const char* general_kwds[];
	static const char* position_kwds[];
	static const char* numeric_kwds[];
	static const char* fluidgrid_kwds[];
	static const char* friction_kwds[];
	static const char* springforce_kwds[];

public:
	
	optionsslipper_file() : input_file("./input/options_slipper.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class materials_file : public input_file
{

	static const char* material_kwds[];

public:

	materials_file() : input_file("./input/materials.txt") {}
	int read(input_data& data);
};
// ------------------------------------------------------------------------- //
class thermal_file : public input_file
{

	const char* filename;
	static const char* general_kwds[];	
	static const char* boundary_kwds[];	
	static const char* dirichlet_kwds[];
	static const char* mixed_kwds[];
	static const char* neumann_kwds[];
	static const char* constraints_kwds[];	

public:

	thermal_file(const char* name) : input_file(name), filename(name) {}
	int read(input_data& data);

};
# endif