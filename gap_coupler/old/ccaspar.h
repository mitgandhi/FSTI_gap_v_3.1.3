#include <iostream>
#include <ppl.h>
#include <algorithm>
#include <array>
#include <vector>
#include <concrt.h>
#include <concrtrm.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>

//Read in CASPAR Gap Modules through the form of dll
//Piston
__declspec(dllexport) int piston_standalone(void);
__declspec(dllexport) int piston_GUI(Concurrency::event & local_event, Concurrency::event & global_event);
//Block
__declspec(dllexport) int run_block_standalone(int argc, char *argv[]);
__declspec(dllexport) void run_block_gap(Concurrency::event & local_event, Concurrency::event & global_event);
//Slipper
//#include "caspar_slipper_dll.h"
__declspec(dllexport) void run_slipper_gap(Concurrency::event & local_event, Concurrency::event & global_event);
//This function will return the current version and date/time stamp of the caspar_slipper.dll
__declspec(dllexport) const char* slipper_version();
//This function can be called to run a standalone Caspar Slipper simulation executable
__declspec(dllexport) int slipper_standalone(int argc, char *argv[]);

using namespace std;
using namespace Concurrency;

void All();

//For reading inputs
struct PistonInput
{
	int mode, npistons, nrevolutions, oiltype;
	double speed, dB, dK, beta, gamma, betamax, pHP, pLP, pCase, THP, TLP, TCase, oildensity, oilbetaP, oilbetaT, oilK;
	double oilbetaKP, oilbetaKT, oilviscosity, oilW, oilTc1, oilPc1, oilTc2, oilPc2, oillambda, oilC, alpha1, alpha2, alpha3; 
};