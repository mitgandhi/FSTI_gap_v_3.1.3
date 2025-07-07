# include "caspar_input.h"
# include <iostream>

using namespace std;

const double pi = 3.14159265358979323846;

void caspar_input::readpfile(void)
{

	char str[256];
	char cbuffer;
	double ttemp;
	double pDCtemp;
	double pHPtemp;
	double pLPtemp;

	ifstream fpFile(operating_conditions.pModuleFile);

	if (!fpFile.is_open())
	{
		inputlog << endl << "ERROR: CGapInput::readpFile -> Unable to open pFile ( " << operating_conditions.pModuleFile << " ) input file!" << endl;
		exit(1);
	}
	
	while(!fpFile.eof()) 
	{
		if( fpFile.get(cbuffer) && (cbuffer =='/'))
		{
			fpFile.getline(str,sizeof(str));continue;
		}
		else
		{
			fpFile.putback(cbuffer);
		}
			
		fpFile >> ttemp >> pDCtemp >> pHPtemp >> pLPtemp;
		pfile.time.push_back(ttemp);
		pfile.pDC.push_back(pDCtemp);
		pfile.pHP.push_back(pHPtemp);
		pfile.pLP.push_back(pLPtemp);
	}
	fpFile.close();
	
	//GapLog.message("CGapInput::readpFile -> pFile input file read successfully!");
}
caspar_input::caspar_input(const input & I) : input_data(I.data)
{
	//set the parent input
	parent_input = &I;

	//read the pressure file
	readpfile();

	//set the revolution period
	common.rev_period = (2.0*pi)/operating_conditions.speed;

	//The solver timestep
	common.phistep_deg = options_slipper.numeric.phistep_deg; 
	//common.phistep_deg = 1.0; //force a 1 deg timestep for now, should be input
	common.timestep = common.rev_period*(common.phistep_deg/360.0);
	
	//Check that our phi step is okay
	{
		double drev_steps = 360.0/common.phistep_deg;
		if(floor(drev_steps) != drev_steps)
		{
			inputlog << "ERROR: 360 is not evenly divisible by the timestep! Change the simulation timestep." << endl;
			exit(1);
		}
	
		common.rev_steps = int(floor(drev_steps));
	}

	//Set the phi deg tol
	common.phi_deg_tol = 0.01;

};
