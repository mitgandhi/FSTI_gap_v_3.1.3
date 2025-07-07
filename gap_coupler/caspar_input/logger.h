# ifndef __logger__
# define __logger__

# include <fstream>
# include <string>
# include <iostream>

using namespace std;

class logger 
{
	std::ofstream logstream;
	std::string gettime(void);
	std::string workingdir(void);
	string logger::getcomputername(void);

public:

	//A logger stream 'manipulator' that indicates an error
	class error
	{
		const int err_code;
	public:
		error(const int error_code = 1) : err_code(error_code)
		{
		};

		//This will cause the process to exit when the object is distructed
		~error()
		{
			exit(err_code);
		}
	};

	//logger constructor / destructor
	logger (const string logFile, const bool Silent = true);
	~logger ();

	//if true, logging will only be directed to the input_log.txt
	//if false, logging will also be copied to inputlog
	bool silent;

	//insertion operators

	//handles ostream
	logger& operator<<(ostream&(*f) (ostream&));
	
	//handles the custom error manipulator
	logger& operator<< (const error& e)
	{
		*this << endl << endl;
		*this << "ERROR: ";
		return *this;
	};


	//this handles everything else but errors
	template<typename T> logger& operator<<(const T& str)
	{
		//should we send the output to the console?
		if(!silent)
		{ 
			//send the output to console
			std::cout << str;
		}

		//if the log file is open, also to the log file
		if(logstream.is_open())
		{
			logstream << str;
		}

		return *this;
	};
		
	//flush the underlying streams
	void flush()
	{
		logstream.flush();
		cout.flush();
	};

};

# endif
