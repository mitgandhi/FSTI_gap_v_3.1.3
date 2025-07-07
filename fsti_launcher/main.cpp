#include "main.h"

using namespace std;

int main(int argc, char* argv[])
{
	//THIS IS THE MAIN FSTI_GAP.EXE ENTRY POINT
	//WE NEED A SEPERATE EXE BECAUSE ALL THE LOG OBJECTS ARE USED AS EXTERN C OBJECTS
	//THUS WE ARE UNABLE TO MOVE THEM BEFORE THEY ARE CREATED UNLESS A SEPERATE (DELAYED LOADED) DLL
	//IS USED TO LOAD THEM

	//THIS FUNCTION SIMPLY MOVES THE OLD OUPUTS / LOGS IF THEY ALREADY EXIST
	//AND THEN LAUNCHES THE FSTI_MAIN DLL FUNCTION

	//COMPILE WITH THE FOLLOWING FLAGS ON THE LINKER INPUT PROPERTY PAGE
	//Delay Loaded Dlls: fsti_gap.dll
	//Additional Dependencies: Delayimp.lib

	if(FileExists("./output/info.txt"))
	{
		//we need to move this dir and log files
		cout << "Moving old outputs directory ... " << endl;

		//let's see if an info.txt file exists for the new dir name
		ifstream info("./output/info.txt");
		
		string newoutput = "";
		if(info.is_open())
		{
			info >> newoutput;
			info.close();
		}

		if(newoutput.length() == 0)
		{
			//set to the current time
			newoutput = getinfostring();
		}

		//try to move the output directory
		newoutput = "./Outputs from " + newoutput;
		MoveFileEx("./output", newoutput.c_str(), NULL);

		//test to see if the dir has moved
		if(dirExists(newoutput))
		{
			//move the rest of the log files in there
			cout << "Moving old log files ... " << endl;

			//this file will have been renamed if it exists
			MoveFileEx("./gap_log.txt", (newoutput+"/gap_log.txt").c_str(), NULL);
			MoveFileEx("./input_log.txt", (newoutput+"/input_log.txt").c_str(), NULL);
			MoveFileEx("./slipper_log.txt", (newoutput+"/slipper_log.txt").c_str(), NULL);
			MoveFileEx("./block_log.txt", (newoutput+"/block_log.txt").c_str(), NULL);
			MoveFileEx("./piston_log.txt", (newoutput+"/piston_log.txt").c_str(), NULL);
		}

	}

	cout << "Creating output directory ... " << endl;
	CreateDirectory("./output", NULL);
	ofstream info("./output/info.txt");
	if(info.is_open())
	{
		//write the date / time
		info << getinfostring() << endl;
		info.close();
	}

	//Launch the main fsti dll
	return fsti_main(argc, argv);
	
}