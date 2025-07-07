#include <iostream>
#include <string>
#include "caspar_input\input.h"
#include <windows.h>
#include <time.h>
#include <ppl.h>

//these are the main lubrication dll methods

///the slipper
void run_slipper_gap(Concurrency::event & local_event, Concurrency::event & global_event, const input & I);

//the piston
int piston_GUI(Concurrency::event & local_event, Concurrency::event & global_event);

//the block
void run_block_gap(Concurrency::event & local_event, Concurrency::event & global_event, const input& in);

//the block "standalone" functions used for pre-check
int run_block_standalone(int argc, char* argv[]);

#ifndef DATA_PROCESSOR_H
#define DATA_PROCESSOR_H

#include <vector>
#include <string>

// Function to read temperature data from a file and store it in a vector
//std::vector<double> readTemperatureData(const std::string& filePath);
std::vector<double>* readTemperaturesFromFile(const std::string& filePath, std::string& errorMessage);


// Function to create a summary file with mean, min, and max temperatures
void createSummaryFile(const std::string& outputFilePath, const std::vector<std::vector<double>>& data);

#endif

bool FileExists(const std::string szPath)
{
  DWORD dwAttrib = GetFileAttributes(szPath.c_str());

  return (dwAttrib != INVALID_FILE_ATTRIBUTES && 
         !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
};

bool fileExists(const std::string& filePath) {
    std::ifstream file(filePath);
    return file.good();
}