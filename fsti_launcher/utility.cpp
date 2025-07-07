#include "main.h"

using namespace std;

bool dirExists(const std::string& dirName_in)
{
  DWORD ftyp = GetFileAttributesA(dirName_in.c_str());
  if (ftyp == INVALID_FILE_ATTRIBUTES)
    return false;  //something is wrong with your path!

  if (ftyp & FILE_ATTRIBUTE_DIRECTORY)
    return true;   // this is a directory!

  return false;    // this is not a directory!
}

string getinfostring(void)
{
	time_t rawtime;
	struct tm timeinfo;
	char buffer [80];
	time ( &rawtime );
	localtime_s (&timeinfo, &rawtime);
	strftime (buffer,80,"%m%d%y_%H%M%S",&timeinfo);
	return buffer;
};

bool FileExists(const std::string szPath)
{
  DWORD dwAttrib = GetFileAttributes(szPath.c_str());

  return (dwAttrib != INVALID_FILE_ATTRIBUTES && 
         !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
};

