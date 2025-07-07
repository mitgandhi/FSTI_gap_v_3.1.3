# include "./section.h"
# include <iostream>
# include <string>
# include <vector>
# include <sstream>

using namespace std;


//we need all of this because the C++ standards library sucks when dealing with strings and can't include it
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}


// ------------------------------------------------------------------------- //
void section::push_back(std::istringstream& field)
{
	string s;
	string v;
	field >> s;
	if(s.size() > 0 && s.find("//") == string::npos)
	{
		//instead of using the methodology to remove extra whitespace between values now, 
		//let's just shove the rest of the line into the value field
		//splitting will be done by the insertion operators in the other read() methods

		/*
		string tmp;
		do
		{
			tmp.clear();
			field >> tmp;
			if(tmp.size() > 0)
			{
				v.append(tmp);
				v.append(string("\t"));
			}
		} while(tmp.size() > 0);
		*/

		getline(field, v);
		//we can trim a string, that should be fine
		v = trim(v);

		//i think we should still push the keyword back, even if it is empty.
		//this is especially important so that we can have empty file names
		//thus i'm commenting out the if() statement - andrew
		
		// check if the value is not empty
		//if(v.size() > 0)
		{
			keywords.push_back(s);
			values.push_back(v);
		}
	}
}
// ------------------------------------------------------------------------- //
int section::find_keyword(char* word)
{
	int idx = -1;
	for(int i = 0; i < keywords.size(); i++)
	{
		if(keywords[i].compare(word) == 0)
			idx = i;
	}
	return idx;
}
// ------------------------------------------------------------------------- //
void section::check_keywords(const char* list[], unsigned int sz, const char* file)
{
	
	vector<string> missing(0);

	for(unsigned int i = 0; i < sz; i++)
	{
		bool found = false;
		for(unsigned int j = 0; j < keywords.size(); j++)
		{
			if(string(list[i]).compare(keywords[j]) == 0)
			{
				found = true; 
				break;
			}
		}	
		if(!found)
		{
			inputlog << logger::error()
					 << "\n\nIn file " << file << " (section " << name << ")\n" 
					 << "keyword \"" << list[i]  
					 << "\" was not found, or an invalid value "
					 << "is associated with it." << endl << endl;

		}
	}	
}