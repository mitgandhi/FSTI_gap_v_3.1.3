#include "ccaspar.h"


void main()
{
	//CASPAR intro banner
	cout << "********************************************************************************";
	cout << "*" << setw(78) << " " << "*";
	cout << "*" << "\t" << left << setw(71) << "CASPAR Sunchronization" << "*";
	cout << "*" << "\t" << left << setw(71) << "June 2012 by Dongjune Kim" << "*";
	cout << "*" << setw(78) << " " << "*";
	cout << "*" << "\t" << left << setw(71) << "Purdue University" << "*";
	cout << "*" << "\t" << left << setw(71) << "Maha Fluid Power Research Center" << "*";
	cout << "*" << "\t" << left << setw(71) << "https://engineering.purdue.edu/Maha" << "*";
	cout << "********************************************************************************";
	cout << endl;

	bool quit = false;
	while(!quit)
	{
		cout << endl;
		cout << "Please select an option:" << endl;
		cout << endl;
		cout << "\t" << "1   Piston / Cylinder" << endl;
		cout << "\t" << "2   Cylinder Block / Valve plate" << endl;
		cout << "\t" << "3   Slipper / Swash plate" << endl;
		cout << "\t" << "4   All the above" << endl;
		cout << "\t" << "5   Quit" << endl;
		cout << endl;
		cout << "Enter Option Number: ";
		string input;
		int option;
		cin >> input;
		option = atoi(input.c_str());

		system("cls");

		cout << endl << "********************************************************************************" << endl;
		switch(option)
		{
		case 1:
			//Piston();
			break;
		case 2:
			//CylinderBlock();
			break;
		case 3:
			//slipper_standalone(int argc, char *argv[]);
			break;
		case 4:
			//Coupled interface simulation
			All();
			break;
		case 5:
			cout << "Quitting CASPAR..." << endl;
			quit = true;
			break;
		default:
			cout << "ERROR: Invalid Selection!" << endl;
			break;
		}
		cout << endl << "********************************************************************************" << endl;
	}
	cout << endl;

	cout << " Byebye " << endl;
}