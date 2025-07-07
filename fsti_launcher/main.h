#include <iostream>
#include <fstream>
#include <string>
#include <windows.h>
#include <time.h>

//utility function declerations
bool dirExists(const std::string& dirName_in);
std::string getinfostring(void);
bool FileExists(const std::string szPath);


//this is the function decleration for the gap_coupler main hook
int fsti_main(int argc, char* argv[]);
