# ifndef __section__
# define __section__

# include <vector>
# include <sstream>
# include "logger.h"

//the global input log
extern logger inputlog;

struct section
{

	std::string name;
	std::vector<std::string> keywords;
	std::vector<std::string> values;
	
	// member functions
	section() {}
	void push_back(std::istringstream& field);
	int find_keyword(char* word);
	void check_keywords(const char* list[], unsigned int sz, const char* file);
};

# endif