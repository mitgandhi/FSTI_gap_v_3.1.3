#include "logger.h"
#include "time.h"
#include <Windows.h>
#include <stdio.h>
#include <direct.h>
#include <iomanip>  // For setw
using namespace std;

// Include dynamically defined version
#include "version_define.h"
#ifndef PROJECT_VERSION
    #define PROJECT_VERSION "Unknown"
#endif

logger::logger(const string logFile, const bool Silent) : silent(Silent)
{
    logstream.open(logFile);
    if(!logstream.is_open())
    {
        cout << "ERROR: Unable to open " << logFile << " !" << endl;
        exit(1);
    }
    
    // Calculate padding needed to format version number nicely
    string version = PROJECT_VERSION;
    int padding = 51 - version.length();  // 51 is column width minus asterisks
    
    logstream << "*************************************************************************" << endl;
    logstream << "*                                                                       *" << endl;
    logstream << "* FSTI Gap Design                                                       *" << endl;
    logstream << "* Version: " << version << setw(padding) << "*" << endl;
    logstream << "* Smart Hydraulic Solutions                                             *" << endl;
    logstream << "* SHS                                                                   *" << endl;
    logstream << "* https://www.smart-hydraulic-solutions.com/casparfsti                  *" << endl;
    logstream << "* FSTI Version: " << version << setw(padding - 12) << "*" << endl;  // -12 to account for "FSTI Version: "
    logstream << "*                                                                       *" << endl;
    logstream << "*************************************************************************" << endl;
    logstream << endl;
    logstream << "Time: " << gettime() << endl;
    logstream << "Computer: " << getcomputername() << endl;
    logstream << "Working Directory: " << workingdir() << endl;
    logstream << endl;
    logstream << "Log Opened ..." << endl;
    logstream << endl;
}
logger::~logger(void) 
{
	if(logstream.is_open())
	{
		logstream << endl;
		logstream << "Time: " << gettime() << endl;
		logstream << "... Log Closed." << endl;
		logstream.close();
	}
}

//this handles endl, scientific, setprecision, etc..
logger& logger::operator<<( ostream&(*f) (ostream&))
{
	//should we send the output to the console?
	if(!silent)
	{
		//send the output to console
		cout << f;
	}

	//if the log file is open, also to the log file
	if(logstream.is_open())
	{
		logstream << f;
	}

	return *this;
}
//get the current date / time
string logger::gettime(void)
{
	time_t rawtime;
	struct tm timeinfo;
	char buffer [80];
	time ( &rawtime );
	localtime_s (&timeinfo, &rawtime);
	strftime (buffer,80,"%x %X",&timeinfo);
	return buffer;
};
//get the working program dir
string logger::workingdir(void)
{
	char cdir[FILENAME_MAX] = {0};
	_getcwd(cdir, sizeof(cdir));

	return cdir;
};
//get the windows system computer name
string logger::getcomputername(void)
{
	unsigned long computerNameLength = 256;
	char name[256];
	GetComputerNameA(name, &computerNameLength);
	return string(name);
};
