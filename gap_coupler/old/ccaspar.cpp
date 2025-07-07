#include "ccaspar.h"

//Function that reads revolution from input file
int CheckInputs()
{
	int rev;
	PistonInput PInput;

	// open the input file
	ifstream pgeneral("./inputs_piston/General.gen");

	cout << "[Coupling]\t\tCheck Input files.................................." << endl;
	
	if (!pgeneral.is_open()) 	{
		cout << "[Coupling]\t\tUnable to open General.gen! - Generate General.gen and restart..." << endl;
		system("PAUSE");
		exit(1);
	}

	string tmp;
	string line;

	while(pgeneral) {
		getline(pgeneral,line);
		istringstream iss (line,istringstream::in);
		while(iss) {
			iss >> tmp;
			
			// check if there is a comment
			size_t comment  = tmp.find("//");
			if (comment!=string::npos)
			{
				break;
			}
			if (tmp == "mode"){
				iss >> PInput.mode;
				break;
			}
			if (tmp == "npistons"){
				iss >> PInput.npistons;
				break;
			}
			if (tmp == "nrevolutions"){
				iss >> PInput.nrevolutions;
				break;
			}
			if (tmp == "speed"){
				iss >> PInput.speed;
				break;
			}
			if (tmp == "dB"){
				iss >> PInput.dB;
				break;
			}
			if (tmp == "dK"){
				iss >> PInput.dK;
				break;
			}
			if (tmp == "beta"){
				iss >> PInput.beta;
				break;
			}
			if (tmp == "gamma"){
				iss >> PInput.gamma;
				break;
			}
			if (tmp == "betamax"){
				iss >> PInput.betamax;
				break;
			}
			if (tmp == "pHP"){
				iss >> PInput.pHP;
				break;
			}
			if (tmp == "pLP"){
				iss >> PInput.pLP;
				break;
			}
			if (tmp == "pCase"){
				iss >> PInput.pCase;
				break;
			}
			if (tmp == "THP"){
				iss >> PInput.THP;
				break;
			}
			if (tmp == "TLP"){
				iss >> PInput.TLP;
				break;
			}
			if (tmp == "TCase"){
				iss >> PInput.TCase;
				break;
			}
			if (tmp == "oiltype"){
				iss >> PInput.oiltype;
				break;
			}
			if (tmp == "oildensity"){
				iss >> PInput.oildensity;
				break;
			}
			if (tmp == "oilbetaP"){
				iss >> PInput.oilbetaP;
				break;
			}
			if (tmp == "oilbetaT"){
				iss >> PInput.oilbetaT;
				break;
			}
			if (tmp == "oilK"){
				iss >> PInput.oilK;
				break;
			}
			if (tmp == "oilbetaKP"){
				iss >> PInput.oilbetaKP;
				break;
			}
			if (tmp == "oilbetaKT"){
				iss >> PInput.oilbetaKT;
				break;
			}
			if (tmp == "oilviscosity"){
				iss >> PInput.oilviscosity;
				break;
			}
			if (tmp == "oilW"){
				iss >> PInput.oilW;
				break;
			}
			if (tmp == "oilTc1"){
				iss >> PInput.oilTc1;
				break;
			}
			if (tmp == "oilPc1"){
				iss >> PInput.oilPc1;
				break;
			}
			if (tmp == "oilTc2"){
				iss >> PInput.oilTc2;
				break;
			}
			if (tmp == "oilPc2"){
				iss >> PInput.oilPc2;
				break;
			}
			if (tmp == "oillambda"){
				iss >> PInput.oillambda;
				break;
			}
			if (tmp == "oilC"){
				iss >> PInput.oilC;
				break;
			}
			if (tmp == "alpha1"){
				iss >> PInput.alpha1;
				break;
			}
			if (tmp == "alpha2"){
				iss >> PInput.alpha2;
				break;
			}
			if (tmp == "alpha3"){
				iss >> PInput.alpha3;
				break;
			}
		}
	}

	pgeneral.close();


	rev = PInput.nrevolutions;

	return rev;
}

void All()
{
	cout << "[Coupling]\t\tStart Coupling " << endl;

	//int rev=CheckInputs(), mnum=1;		//Revolution, Number of modules -> Incase reading revolution from input file
	int rev=10, mnum=3;		//Revolution, Number of modules
	CurrentScheduler::Create(SchedulerPolicy(1, MaxConcurrency, 3));
	
	task_group tasks;

	//the local inteface events
	event * local_rev = new event [mnum];	
	//the global sync events
	event * global_rev = new event [mnum];

	//Run Gap Modules
	//Piston/Cylinder
	tasks.run([local_rev, global_rev](){piston_GUI(local_rev[0], global_rev[0]);});
	//CylinderBlock/Valveplate
	tasks.run([local_rev, global_rev](){run_block_gap(local_rev[1], global_rev[1]);});
	//Slipper/Swashplate
	tasks.run([local_rev, global_rev](){run_slipper_gap(local_rev[2], global_rev[2]);});

	for(int i = 0; i < rev; i++)
	{
		for(int j = 0; j < mnum; j++)
		{
			local_rev[j].wait();
			local_rev[j].reset();
		}

		//Trasfer forces
		cout << "[Coupling]\t\t------------------------------------- " << endl;
		cout << "[Coupling]\t\tRevolution #: " << i+1 << endl;
		
		for(int j = 0; j < mnum; j++)
		{
			global_rev[j].set();
		}
	}

	tasks.wait();

	cout << "[Coupling]\t\tFinished" << endl;

	delete [] local_rev;
	delete [] global_rev;
}
