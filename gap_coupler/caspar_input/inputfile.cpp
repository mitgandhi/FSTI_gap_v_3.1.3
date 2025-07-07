# include "inputfile.h"
# include <iostream>

using namespace std;

// ------------------------------------------------------------------------- //
input_file::input_file(const char* filename)
{

	isopen = read(filename);

}
// ------------------------------------------------------------------------- //
input_file::~input_file()
{
	for(unsigned int i=0; i<sections.size(); i++)
		delete sections[i];
}
// ------------------------------------------------------------------------- //
int input_file::check_sections()
{
	return 0;
}
// ------------------------------------------------------------------------- //
bool input_file::read(const char* filename)
{
	file_name = filename;
	file.open(filename);

	if (!file.is_open()) 	
	{
		return false;
	}

	inputlog << "Reading " << filename << " input file ... ";

	string tmp;
	string line;
	double val = 0;
	int linenumber = 1;

	
	while(getline(file,line)) 
	{
		if(line.size() > 0)
		{

			istringstream iss (line, istringstream::in);
			iss >> tmp;
			size_t comment  = tmp.find("//"); // check if there is a comment
			
			// if is not a comment read
			if (comment == string::npos) 
			{
				// new section
				if(tmp.compare("section") == 0)
				{
					section* s = new section;
					sections.push_back(s);
					iss >> tmp;	
					if(tmp.compare("section") != 0)
					{
						s -> name = tmp;
						// get the keywords
						do
						{
							getline(file,line);
							if(line.find("endsection") == string::npos)
								s->push_back(istringstream(line));
							else
								break;
						} while(true);
					}
					else
					{
						inputlog << logger::error() << "\nIn file " << filename << " section name is missing!" 
								 << endl;
					}
				}
			}
		}
	}

	inputlog << "done!" << endl;

	return true;
}
