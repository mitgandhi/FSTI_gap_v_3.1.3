
#include "main.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <numeric>
#include <regex>
#include <iterator>
#include <direct.h> // For directory functions
#include <sys/stat.h> // For stat structure
#include <string>
#include <limits>
#include <cctype>
#include <cmath> // For fabs() function
#include <vector>
using namespace std;
using namespace Concurrency;

logger logfile("gap_log.txt", false);
bool ininputs = false;
bool shouldexit = false;
const double TEMPERATURE_CONVERGENCE_THRESHOLD = 1.0; // Maximum temperature difference (Celsius) for convergence
const double pi = 3.14159265358979323846;

// Progress tracking structures and variables
struct InterfaceProgress {
    int current_rev;
    double current_phi;
    double progress_percentage;
    
    // New fields for real-time estimation
    clock_t calculation_start_time;
    bool calculation_in_progress;
    double estimated_progress;
    double avg_step_time;     // Average time per phi step
    double last_phi;          // Last recorded phi value
    clock_t last_update_time; // Time of last progress update
    
    InterfaceProgress() : current_rev(0), current_phi(0), progress_percentage(0),
                         calculation_start_time(0), calculation_in_progress(false),
                         estimated_progress(0), avg_step_time(1.0), last_phi(0),
                         last_update_time(0) {}
};


// New timing variables for periodic updates
clock_t last_progress_display = 0;
const double PROGRESS_UPDATE_INTERVAL = 10.0; // 10 seconds in seconds

// Array to store progress for all interfaces
InterfaceProgress interface_progress[3];

// Store previous temperature data for each interface
vector<double> prev_max_temp(3, 0.0);

// Function to check if directory exists
bool directoryExists(const string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
}



string gt() {
    time_t rawtime;
    struct tm timeinfo;
    char buffer[80];
    time(&rawtime);
    localtime_s(&timeinfo, &rawtime);
    strftime(buffer, 80, "%x %X", &timeinfo);
    return string("[") + buffer + "] ";
}

//------------------------------------------
// 6. Updated: getProgressBar() to show estimated
std::string getProgressBar(double percentage, int width = 20, bool estimated = false) {
    int completedWidth = static_cast<int>(percentage * width / 100.0);
    std::ostringstream oss;
    oss << "[";

    for (int i = 0; i < width; i++) {
        if (i < completedWidth) {
            oss << "=";
        } else if (i == completedWidth) {
            oss << ">";
        } else {
            oss << " ";
        }
    }

    oss << "] ";
    oss << std::setw(3) << static_cast<int>(percentage) << "%";
    if (estimated) oss << " (est)";
    return oss.str();
}


// Calculate continuous progress (0-100%)
double calculateContinuousProgress(int currentRev, double currentPhi, int totalRev) {
    if (totalRev <= 0) return 0.0;
    double fractionComplete = currentPhi / 360.0;
    double totalProgress = ((double)currentRev + fractionComplete) / (double)totalRev * 100.0;
    return min(totalProgress, 100.0);
}


// Function to mark the start of a calculation for an interface
void startInterfaceCalculation(int interface_idx) {
    interface_progress[interface_idx].calculation_start_time = clock();
    interface_progress[interface_idx].calculation_in_progress = true;
    interface_progress[interface_idx].last_update_time = clock();
}

// Function to mark the end of a calculation
void endInterfaceCalculation(int interface_idx) {
    interface_progress[interface_idx].calculation_in_progress = false;
}

// Function to update the estimates for calculation time (call this when actual progress is recorded)
void updateCalculationTimeEstimates(int interface_idx, double new_phi) {
    if (!interface_progress[interface_idx].calculation_in_progress)
        return;

    double phi_change = new_phi - interface_progress[interface_idx].last_phi;

    // Handle wrap-around safely
    if (phi_change <= 0) {
        if (new_phi < interface_progress[interface_idx].last_phi && 
            fabs(new_phi - interface_progress[interface_idx].last_phi) > 180.0) {
            phi_change = (360.0 - interface_progress[interface_idx].last_phi) + new_phi;
        } else {
            return; // ignore regression
        }
    }

    clock_t current_time = clock();
    double elapsed_seconds = (current_time - interface_progress[interface_idx].last_update_time) / (double)CLOCKS_PER_SEC;

    if (elapsed_seconds > 0.1) {
        double time_per_degree = elapsed_seconds / phi_change;

        if (interface_progress[interface_idx].avg_step_time > 0) {
            interface_progress[interface_idx].avg_step_time =
                0.7 * interface_progress[interface_idx].avg_step_time +
                0.3 * time_per_degree;
        } else {
            interface_progress[interface_idx].avg_step_time = time_per_degree;
        }

        interface_progress[interface_idx].last_phi = new_phi;
        interface_progress[interface_idx].last_update_time = current_time;
    }
}
// Function to estimate current progress based on time elapsed
void updateEstimatedProgress(int interface_idx, int total_rev) {
    // Only update if calculation is in progress
    if (!interface_progress[interface_idx].calculation_in_progress)
        return;
        
    // Get elapsed time since last known position update
    clock_t current_time = clock();
    double elapsed_seconds = (current_time - interface_progress[interface_idx].last_update_time) 
                           / (double)CLOCKS_PER_SEC;
    
    // Estimate progress based on average calculation time
    double avg_time_per_degree = interface_progress[interface_idx].avg_step_time;
    
    // Estimate phi progress based on elapsed time
    double estimated_phi_progress = (elapsed_seconds / avg_time_per_degree);
    
    // Ensure estimated progress doesn't exceed a full revolution
    double max_estimated_phi = 360.0 - interface_progress[interface_idx].last_phi;
    estimated_phi_progress = min(estimated_phi_progress, max_estimated_phi);
    
    // Calculate total estimated phi
    double estimated_phi = interface_progress[interface_idx].last_phi + estimated_phi_progress;
    
    // Make sure we don't exceed 360 degrees
    if (estimated_phi > 360.0) {
        estimated_phi = 360.0;
    }
    
    // Calculate estimated progress percentage
    interface_progress[interface_idx].estimated_progress = 
        calculateContinuousProgress(interface_progress[interface_idx].current_rev, 
                                   estimated_phi, total_rev);
}

// Function to update the progress information for a specific interface
void updateInterfaceProgress(int interface_idx, int rev_value, double phi_value, int total_rev) {
    // Store previous phi to calculate step time
    double prev_phi = interface_progress[interface_idx].current_phi;
    
    // Update with actual progress values
    interface_progress[interface_idx].current_rev = rev_value;
    interface_progress[interface_idx].current_phi = phi_value;
    interface_progress[interface_idx].progress_percentage = 
        calculateContinuousProgress(rev_value, phi_value, total_rev);
        
    // Update calculation time estimates
    updateCalculationTimeEstimates(interface_idx, phi_value);
}
//CLEAR THE CONSOLE
void clearConsoleLine() {
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hStdout == INVALID_HANDLE_VALUE) return;
    
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    if (!GetConsoleScreenBufferInfo(hStdout, &csbi)) return;
    
    COORD coords = { 0, csbi.dwCursorPosition.Y };
    DWORD written;
    
    // Fill the line with spaces
    FillConsoleOutputCharacter(hStdout, ' ', csbi.dwSize.X, coords, &written);
    // Reset cursor position to start of line
    SetConsoleCursorPosition(hStdout, coords);
}

// Add this function to get the base folder name
std::string getBaseFolderName() {
    char currentPath[MAX_PATH];
    if (GetCurrentDirectory(MAX_PATH, currentPath) == 0) {
        return "Simulation"; // Default if we can't get the path
    }
    
    // Extract just the folder name from the full path
    std::string fullPath(currentPath);
    size_t lastSlash = fullPath.find_last_of("\\/");
    if (lastSlash != std::string::npos) {
        return fullPath.substr(lastSlash + 1);
    }
    return fullPath;
}

// Then modify the updateConsoleTitle function
void updateConsoleTitle(bool piston_enabled, bool block_enabled, bool slipper_enabled,
                       bool piston_disabled, bool block_disabled, bool slipper_disabled,
                       bool use_estimated = false) {
    static char titleBuffer[256];
    
    // Get the base folder name instead of using "Simulation Progress"
    std::string folderName = getBaseFolderName();
    sprintf_s(titleBuffer, sizeof(titleBuffer), "%s - ", folderName.c_str());
    
    char* pos = titleBuffer + strlen(titleBuffer);
    
    if (piston_enabled && !piston_disabled) {
        double progress = use_estimated && interface_progress[0].calculation_in_progress ? 
                         interface_progress[0].estimated_progress : 
                         interface_progress[0].progress_percentage;
        pos += sprintf(pos, "Piston: %.1f%% ", progress);
    }
    
    if (block_enabled && !block_disabled) {
        double progress = use_estimated && interface_progress[1].calculation_in_progress ? 
                         interface_progress[1].estimated_progress : 
                         interface_progress[1].progress_percentage;
        pos += sprintf(pos, "Block: %.1f%% ", progress);
    }
    
    if (slipper_enabled && !slipper_disabled) {
        double progress = use_estimated && interface_progress[2].calculation_in_progress ? 
                         interface_progress[2].estimated_progress : 
                         interface_progress[2].progress_percentage;
        pos += sprintf(pos, "Slipper: %.1f%% ", progress);
    }
    
    SetConsoleTitle(titleBuffer);
}
// 4. Updated: printProgressOnce()
void printProgressOnce(const std::vector<std::string>& progressSegments, int totalWidth = 160) {
    static std::string lastPrintedLine;

    std::ostringstream oss;
    for (size_t i = 0; i < progressSegments.size(); ++i) {
        oss << progressSegments[i] << "  ";
    }

    std::string currentLine = oss.str();
    if (currentLine.length() < static_cast<size_t>(totalWidth)) {
        currentLine += std::string(totalWidth - currentLine.length(), ' ');
    } else {
        currentLine = currentLine.substr(0, totalWidth);
    }

    if (currentLine == lastPrintedLine) return;

#ifdef _WIN32
    // std::cout << "\r" << currentLine << std::flush; // Commented to suppress output
#else
    // if (!lastPrintedLine.empty()) {
    //     std::cout << "\033[F";  // move cursor up
    // }
    // std::cout << "\r" << currentLine << std::flush; // Commented to suppress output
#endif

    lastPrintedLine = currentLine;
}

void displayEstimatedProgress(int total_rev, logger& logfile, 
                              const std::vector<bool>& interface_enabled, 
                              const std::vector<bool>& interface_disabled) {
    static const int LINE_WIDTH = 160;
    static clock_t last_log_time = 0;
    std::vector<std::string> segments;

    for (int i = 0; i < 3; ++i) {
        if (!interface_enabled[i] || interface_disabled[i]) {
            continue;
        }

        std::string interface_name;
        if (i == 0) interface_name = "Piston";
        else if (i == 1) interface_name = "Block";
        else interface_name = "Slipper";

        double progress = interface_progress[i].calculation_in_progress ?
                          interface_progress[i].estimated_progress :
                          interface_progress[i].progress_percentage;

        std::ostringstream seg;
        seg << interface_name << " " << getProgressBar(progress)
            << " Rev: " << interface_progress[i].current_rev
            << "+" << std::fixed << std::setprecision(0)
            << (interface_progress[i].current_phi / 3.6) << "% of " << total_rev;

        if (interface_progress[i].calculation_in_progress) {
            seg << " (calculating...)";
        }

        segments.push_back(seg.str());
    }

    // Print to console in-place
    printProgressOnce(segments, LINE_WIDTH);

    // Log to file every minute (optional)
    clock_t now = clock();
    if ((now - last_log_time) > CLOCKS_PER_SEC * 10) { // update every 10 seconds
        std::ostringstream flat;
        for (size_t i = 0; i < segments.size(); ++i) {
            flat << gt() << segments[i] << std::endl;
        }
        //std::cout << flat.str();
        last_log_time = now;
    }
}



// Function to read the last line from a file
bool readLastLineFromFile(const string& filepath, string& lastLine) {
    ifstream file(filepath);
    if (!file.is_open()) return false;

    string currentLine;
    while (getline(file, currentLine)) {
        if (currentLine.empty() || currentLine.find_first_not_of(" \t\r\n") == string::npos)
            continue;
        lastLine = currentLine;
    }

    file.close();
    return !lastLine.empty();
}
// Function to update progress information from output files
void updateProgressFromFiles(bool piston_enabled, bool block_enabled, bool slipper_enabled,
                           bool piston_disabled, bool block_disabled, bool slipper_disabled,
                           int total_rev, bool updateConsole) {
    // Process piston progress
    if (piston_enabled && !piston_disabled) {
        string lastLine;
        if (readLastLineFromFile("output/piston/piston.txt", lastLine)) {
            istringstream iss(lastLine);
            vector<double> vals;
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(vals));
            
            if (vals.size() >= 3) {
                // Get revolution and phi from the data
                int rev_value = static_cast<int>(vals[1]);
                double phi_value = vals[2];
                
                // Update progress information
                updateInterfaceProgress(0, rev_value, phi_value, total_rev);
            }
        }
    }
    
    // Process block progress (same implementation as before)
    if (block_enabled && !block_disabled) {
        string lastLine;
        if (readLastLineFromFile("output/block/block.txt", lastLine)) {
            istringstream iss(lastLine);
            vector<double> vals;
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(vals));
            
            if (vals.size() >= 3) {
                // Get revolution and phi from the data
                int rev_value = static_cast<int>(vals[1]);
                double phi_value = vals[2];
                
                // Update progress information
                updateInterfaceProgress(1, rev_value, phi_value, total_rev);
            }
        }
    }
    
    // Process slipper progress (same implementation as before)
    if (slipper_enabled && !slipper_disabled) {
        string lastLine;
        if (readLastLineFromFile("output/slipper/slipper.txt", lastLine)) {
            istringstream iss(lastLine);
            vector<double> vals;
            copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(vals));
            
            if (vals.size() >= 3) {
                // Get revolution and phi from the data
                int rev_value = static_cast<int>(vals[1]);
                double phi_value = vals[2];
                
                // Update progress information
                updateInterfaceProgress(2, rev_value, phi_value, total_rev);
            }
        }
    }
    
    // Update console title if requested
    if (updateConsole) {
        updateConsoleTitle(piston_enabled, block_enabled, slipper_enabled,
                         piston_disabled, block_disabled, slipper_disabled);
    }
}

//progress bar ends..........//


// Function to calculate efficiency metrics for the piston interface
void calculatePistonEfficiency(int revolutionNumber, logger& logfile) {
    // Get operating conditions
    input gapinput;
    
    // Open piston.txt
    std::ifstream pistonFile("output/piston/piston.txt");
    if (!pistonFile.is_open()) {
        logfile << gt() << "Warning: Could not open piston.txt for efficiency calculation." << std::endl;
        return;
    }
    
    // Read the file and find revolution data
    std::vector<std::vector<double>> revData;
    std::string line;
    
    // Read all lines to find data for the requested revolution
    while (std::getline(pistonFile, line)) {
        std::istringstream iss(line);
        std::vector<double> values;
        double value;
        
        while (iss >> value) {
            values.push_back(value);
        }
        
        // Check if this line has enough values and if it belongs to our revolution
        if (values.size() >= 62) {
            double lineRev = values[1]; // revolution column
            
            // If we're in the target revolution range, store the data
            if (floor(lineRev) == revolutionNumber) {
                revData.push_back(values);
            }
        }
    }
    
    // Skip calculation if no data found for this revolution
    if (revData.empty()) {
        logfile << gt() << "No data found for piston revolution " << revolutionNumber << std::endl;
        return;
    }
    
    // Calculate DP
    double V = 110.0;
    double speed = gapinput.data.operating_conditions.speed * 60/(2*pi);
    double dp = (gapinput.data.operating_conditions.HP - gapinput.data.operating_conditions.LP)/1e5;
    
    // Calculate efficiency metrics using formulas from the MATLAB code
    double betaNorm = gapinput.data.operating_conditions.beta / gapinput.data.operating_conditions.betamax;
    double flow_theo = speed * V * betaNorm / 1000.0;  // Theoretical flow in l/min
    
    // Calculate average leakage for this revolution
    double leakageSum = 0.0;
    for (size_t i = 0; i < revData.size(); ++i) {
        const std::vector<double>& row = revData[i];
        if (row.size() >= 60) {
            leakageSum += row[59]; // Total_Leakage
        }
    }

    double leakage = leakageSum / revData.size() * 60000.0;  // Convert to l/min
    
    // Calculate flow and volumetric efficiency
    double flow = flow_theo - leakage;
    double nu_vol = flow / flow_theo;
    
    // Calculate power losses
    double PlossF_sum = 0.0;
    double PlossL_sum = 0.0;
    
    for (size_t i = 0; i < revData.size(); ++i) {
        const std::vector<double>& row = revData[i];
        if (row.size() >= 63) {
            PlossF_sum += row[61]; // Mechanical power loss
            
            PlossL_sum += row[62]; // Volumetric power loss
            
        }
    }
    
    double PlossF = PlossF_sum / revData.size();
    double PlossL = PlossL_sum / revData.size();
    double Ploss = PlossF + PlossL;
    
    // Calculate power and overall efficiency
    double P_theo = dp * flow_theo / 600.0;
    double P_actual = P_theo - Ploss / 1000.0;
    double nu_ges = P_actual / P_theo;
    double nu_hm = nu_ges / nu_vol;
    
    // Output the results
    logfile << "\n--- PISTON EFFICIENCY ANALYSIS: Revolution " << revolutionNumber << " ---" << std::endl;
    logfile << "Speed       : " << speed << " rpm" << std::endl;
    logfile << "HP / LP     : " << std::fixed << std::setprecision(2) << gapinput.data.operating_conditions.HP/1e5 << " / " 
            << gapinput.data.operating_conditions.LP/1e5 << " bar (Δp = " << dp << " bar)" << std::endl;
    logfile << "Flow Theo   : " << std::fixed << std::setprecision(2) << flow_theo << " l/min" << std::endl;
    logfile << "Leakage     : " << std::fixed << std::setprecision(2) << leakage << " l/min" << std::endl;
    logfile << "Net Flow    : " << std::fixed << std::setprecision(2) << flow << " l/min" << std::endl;
    logfile << "Ploss Total : " << std::fixed << std::setprecision(2) << Ploss << " W" << std::endl;
    logfile << "Ploss Frict : " << std::fixed << std::setprecision(2) << PlossF << " W" << std::endl;
    logfile << "Ploss Leak  : " << std::fixed << std::setprecision(2) << PlossL << " W" << std::endl;
    logfile << "Power Theo  : " << std::fixed << std::setprecision(2) << P_theo << " kW" << std::endl;
    logfile << "Power Actual: " << std::fixed << std::setprecision(2) << P_actual << " kW" << std::endl;
    logfile << "Volumetric Efficiency : " << std::fixed << std::setprecision(2) << (nu_vol * 100) << " %" << std::endl;
    logfile << "Hydro-Mechanical Efficiency : " << std::fixed << std::setprecision(2) << (nu_hm * 100) << " %" << std::endl;
    logfile << "Total Efficiency : " << std::fixed << std::setprecision(2) << (nu_ges * 100) << " %" << std::endl;
}

// Function to calculate efficiency metrics for the block interface
void calculateBlockEfficiency(int revolutionNumber, logger& logfile) {
    // Get operating conditions
    input gapinput;
    
    // Open block.txt
    std::ifstream blockFile("output/block/block.txt");
    if (!blockFile.is_open()) {
        logfile << gt() << "Warning: Could not open block.txt for efficiency calculation." << std::endl;
        return;
    }
    
    // Read the file and find revolution data
    std::vector<std::vector<double>> revData;
    std::string line;
    
    // Read all lines to find data for the requested revolution
    while (std::getline(blockFile, line)) {
        // Skip header lines
        if (line.empty() || line[0] == '%' || line[0] == 't') {
            continue;
        }
        
        std::istringstream iss(line);
        std::vector<double> values;
        double value;
        
        while (iss >> value) {
            values.push_back(value);
        }
        
        // Check if this line has enough values and if it belongs to our revolution
        if (values.size() >= 25) {
            double lineRev = values[1]; // revolution column
            
            // If we're in the target revolution range, store the data
            if (floor(lineRev) == revolutionNumber) {
                revData.push_back(values);
            }
        }
    }
    
    // Skip calculation if no data found for this revolution
    if (revData.empty()) {
        logfile << gt() << "No data found for block revolution " << revolutionNumber << std::endl;
        return;
    }
    
    // Calculate DP
    double V = 110.0;
    double speed = gapinput.data.operating_conditions.speed * 60/(2*pi);
    double dp =( gapinput.data.operating_conditions.HP - gapinput.data.operating_conditions.LP)/1e5;
    
    // Calculate efficiency metrics
    double betaNorm = gapinput.data.operating_conditions.beta / gapinput.data.operating_conditions.betamax;
    double flow_theo = speed * V * betaNorm / 1000.0;  // Theoretical flow in l/min
    
    // Calculate average leakage and power loss for this revolution
    double leakageSum = 0.0;
    double powerLossSum = 0.0;
    double powerLossLSum = 0.0;
    double powerLossFSum = 0.0;
    
    for (size_t i = 0; i < revData.size(); ++i) {
        const std::vector<double>& row = revData[i];
        if (row.size() >= 25) {
            leakageSum += row[21];     // Leak column
            powerLossSum += row[22];    // Ploss column
            
            if (row.size() >= 24) {
                powerLossLSum += row[23]; // PlossL column
                powerLossFSum += row[24]; // PlossF column
            }
        }
    }

    double leakage = leakageSum / revData.size();
    double powerLoss = powerLossSum / revData.size();
    double powerLossL = powerLossLSum / revData.size();
    double powerLossF = powerLossFSum / revData.size();
    
    // Calculate flow and volumetric efficiency
    double flow = flow_theo - leakage;
    double nu_vol = flow_theo > 0 ? flow / flow_theo : 0;
    
    // Calculate power and overall efficiency
    double P_theo = dp * flow_theo / 600.0;
    double P_actual = P_theo - powerLoss / 1000.0;
    double nu_ges = P_theo > 0 ? P_actual / P_theo : 0;
    double nu_hm = nu_vol > 0 ? nu_ges / nu_vol : 0;
    
    // Output the results
    logfile << "\n--- BLOCK EFFICIENCY ANALYSIS: Revolution " << revolutionNumber << " ---" << std::endl;
    logfile << "Speed       : " << speed << " rpm" << std::endl;
    logfile << "HP / LP     : " << std::fixed << std::setprecision(2) << gapinput.data.operating_conditions.HP/1e5<< " / " 
            << gapinput.data.operating_conditions.LP/1e5 << " bar (Δp = " << dp << " bar)" << std::endl;
    logfile << "Flow Theo   : " << std::fixed << std::setprecision(2) << flow_theo << " l/min" << std::endl;
    logfile << "Leakage     : " << std::fixed << std::setprecision(2) << leakage << " l/min" << std::endl;
    logfile << "Net Flow    : " << std::fixed << std::setprecision(2) << flow << " l/min" << std::endl;
    logfile << "Ploss Total : " << std::fixed << std::setprecision(2) << powerLoss << " W" << std::endl;
    logfile << "Ploss Frict : " << std::fixed << std::setprecision(2) << powerLossF << " W" << std::endl;
    logfile << "Ploss Leak  : " << std::fixed << std::setprecision(2) << powerLossL << " W" << std::endl;
    logfile << "Power Theo  : " << std::fixed << std::setprecision(2) << P_theo << " kW" << std::endl;
    logfile << "Power Actual: " << std::fixed << std::setprecision(2) << P_actual << " kW" << std::endl;
    logfile << "Volumetric Efficiency : " << std::fixed << std::setprecision(2) << (nu_vol * 100) << " %" << std::endl;
    logfile << "Hydro-Mechanical Efficiency : " << std::fixed << std::setprecision(2) << (nu_hm * 100) << " %" << std::endl;
    logfile << "Total Efficiency : " << std::fixed << std::setprecision(2) << (nu_ges * 100) << " %" << std::endl;
 
}
// Function to calculate efficiency metrics for the slipper interface
void calculateSlipperEfficiency(int revolutionNumber, logger& logfile) {
    // Get operating conditions
    input gapinput;
    
    // Open slipper.txt
    std::ifstream slipperFile("output/slipper/slipper.txt");
    if (!slipperFile.is_open()) {
        logfile << gt() << "Warning: Could not open slipper.txt for efficiency calculation." << std::endl;
        return;
    }
    
    // Read the file and find revolution data
    std::vector<std::vector<double>> revData;
    std::string line;
    
    // Read all lines to find data for the requested revolution
    while (std::getline(slipperFile, line)) {
        // Skip header lines
        if (line.empty() || line[0] == '%') {
            continue;
        }
        
        std::istringstream iss(line);
        std::vector<double> values;
        double value;
        
        while (iss >> value) {
            values.push_back(value);
        }
        
        // Check if this line has enough values and if it belongs to our revolution
        if (values.size() >= 18) {
            double lineRev = values[1]; // revolution column
            
            // If we're in the target revolution range, store the data
            if (floor(lineRev) == revolutionNumber) {
                revData.push_back(values);
            }
        }
    }
    
    // Skip calculation if no data found for this revolution
    if (revData.empty()) {
        logfile << gt() << "No data found for slipper revolution " << revolutionNumber << std::endl;
        return;
    }
    
    // Calculate DP
    double V = 110.0;
    double speed = gapinput.data.operating_conditions.speed * 60/(2*pi);
    double dp = (gapinput.data.operating_conditions.HP - gapinput.data.operating_conditions.LP)/1e5;
    
    // Calculate efficiency metrics
    double betaNorm = gapinput.data.operating_conditions.beta / gapinput.data.operating_conditions.betamax;
    double flow_theo = speed * V * betaNorm / 1000.0;  // Theoretical flow in l/min
    
    // Calculate average leakage and power loss for this revolution
    double leakageSum = 0.0;
    double powerLossSum = 0.0;
    double torqueLossSum = 0.0;
    
    for (size_t i = 0; i < revData.size(); ++i) {
        const std::vector<double>& row = revData[i];
        if (row.size() >= 18) {
            leakageSum += row[12];      // QSG column
            torqueLossSum += row[13];   // TorqueLoss column
            powerLossSum += row[14];    // PowerLoss column
        }
    }

    double leakage = leakageSum / revData.size();
    double torqueLoss = torqueLossSum / revData.size();
    double powerLoss = powerLossSum / revData.size();
    
    // Calculate flow and volumetric efficiency
    double flow = flow_theo - leakage;
    double nu_vol = flow_theo > 0 ? flow / flow_theo : 0;
    
    // Calculate power and overall efficiency
    double P_theo = dp * flow_theo / 600.0;
    double P_actual = P_theo - powerLoss / 1000.0;
    double nu_ges = P_theo > 0 ? P_actual / P_theo : 0;
    double nu_hm = nu_vol > 0 ? nu_ges / nu_vol : 0;
    
    // Output the results
    logfile << "\n--- SLIPPER EFFICIENCY ANALYSIS: Revolution " << revolutionNumber << " ---" << std::endl;
    logfile << "Speed       : " << speed << " rpm" << std::endl;
    logfile << "HP / LP     : " << std::fixed << std::setprecision(2) << gapinput.data.operating_conditions.HP/1e5 << " / " 
            << gapinput.data.operating_conditions.LP/1e5 << " bar (Δp = " << dp << " bar)" << std::endl;
    logfile << "Flow Theo   : " << std::fixed << std::setprecision(2) << flow_theo << " l/min" << std::endl;
    logfile << "Leakage     : " << std::fixed << std::setprecision(2) << leakage << " l/min" << std::endl;
    logfile << "Net Flow    : " << std::fixed << std::setprecision(2) << flow << " l/min" << std::endl;
    logfile << "Torque Loss : " << std::fixed << std::setprecision(2) << torqueLoss << " Nm" << std::endl;
    logfile << "Power Loss  : " << std::fixed << std::setprecision(2) << powerLoss << " W" << std::endl;
    logfile << "Power Theo  : " << std::fixed << std::setprecision(2) << P_theo << " kW" << std::endl;
    logfile << "Power Actual: " << std::fixed << std::setprecision(2) << P_actual << " kW" << std::endl;
     logfile <<"Volumetric Efficiency : " << std::fixed << std::setprecision(2) << (nu_vol * 100) << " %" << std::endl;
    logfile << "Hydro-Mechanical Efficiency : " << std::fixed << std::setprecision(2) << (nu_hm * 100) << " %" << std::endl;
    logfile << "Total Efficiency : " << std::fixed << std::setprecision(2) << (nu_ges * 100) << " %" << std::endl;

}

// Main function to analyze efficiency for all active interfaces
void analyzeEfficiencyPerRevolution(int currentRevolution, logger& logfile, 
                                   bool piston_enabled, bool block_enabled, bool slipper_enabled) {
    logfile << gt() << "Beginning efficiency analysis for revolution " << currentRevolution << "..." << std::endl;
    
    // Calculate efficiency for each enabled interface
    if (piston_enabled) {
        calculatePistonEfficiency(currentRevolution, logfile);
    }
    
    if (block_enabled) {
        calculateBlockEfficiency(currentRevolution, logfile);
    }
    
    if (slipper_enabled) {
        calculateSlipperEfficiency(currentRevolution, logfile);
    }
    
    // Add summary information
    input gapinput;
    double speed = gapinput.data.operating_conditions.speed * 60/(2*pi);
    double dp = (gapinput.data.operating_conditions.HP - gapinput.data.operating_conditions.LP)/1e5;
    
    logfile << gt() << "Efficiency analysis complete for revolution " << currentRevolution << std::endl;
    logfile << gt() << "Operating conditions: Speed=" << speed 
            << " rpm, HP=" << gapinput.data.operating_conditions.HP/1e5 << " bar, LP=" << gapinput.data.operating_conditions.LP/1e5
            << " bar, Beta=" << gapinput.data.operating_conditions.beta*180/pi << ", MaxBeta=" << gapinput.data.operating_conditions.betamax*180/pi << std::endl;
    logfile.flush();
}



// Thermal ouputs -----------------------------------------------
struct ThermalStats {
    double T_min;
    double T_mean;
    double T_max;
    double D_min;
    double D_max;
    bool valid;

    ThermalStats()
        : T_min(0), T_mean(0), T_max(0),
          D_min(0), D_max(0), valid(false) {}
};

vector<double> extract_scalar_values(const string& filepath, const string& scalar_name) {
    ifstream vtkfile(filepath.c_str());
    vector<double> values;
    string line;

    while (getline(vtkfile, line)) {
        if (line.find("SCALARS " + scalar_name) != string::npos) {
            while (getline(vtkfile, line)) {
                if (line.find("LOOKUP_TABLE") != string::npos) break;
            }
            while (getline(vtkfile, line)) {
                if (!line.empty() && isalpha(line[0])) break;
                istringstream iss(line);
                double val;
                while (iss >> val) {
                    values.push_back(val);
                }
            }
            break;
        }
    }
    return values;
}


// Enhanced parse_vtk function with better debugging
ThermalStats parse_vtk(const string& filepath) {
    ThermalStats stats;
    
    // Check if file exists first
    if (!fileExists(filepath)) {
        return stats; // Return invalid stats
    }

    vector<double> temps;
    vector<double> defos;
    
    // Check if this is a slipper/swashplate file
    bool is_slipper_interface = (filepath.find("slipper") != string::npos || 
                                filepath.find("swashplate") != string::npos);
    
    if (is_slipper_interface) {
        // For slipper/swashplate: try multiple possible field names
        // C++98 compatible way - try each field name individually
        temps = extract_scalar_values(filepath, "temperature");
        if (temps.empty()) {
            temps = extract_scalar_values(filepath, "Temperature");
        }
        if (temps.empty()) {
            temps = extract_scalar_values(filepath, "T_[C]");
        }
        if (temps.empty()) {
            temps = extract_scalar_values(filepath, "T");
    } 
	}else {
        // For piston/block: use existing field names
        temps = extract_scalar_values(filepath, "T_[C]");
        defos = extract_scalar_values(filepath, "Total_Deformation_[um]");
    }

    if (!temps.empty()) {
        stats.T_min = *min_element(temps.begin(), temps.end());
        stats.T_max = *max_element(temps.begin(), temps.end());
        stats.T_mean = accumulate(temps.begin(), temps.end(), 0.0) / temps.size();
        
        if (!defos.empty()) {
            stats.D_min = *min_element(defos.begin(), defos.end());
            stats.D_max = *max_element(defos.begin(), defos.end());
        } else {
            // No deformation data (slipper files or missing data)
            stats.D_min = 0.0;
            stats.D_max = 0.0;
        }
        stats.valid = true;
    }
    
    return stats;
};
// Check if temperature convergence is achieved for an interface
bool check_temperature_converged(const ThermalStats& stats, int interface_idx, double prev_temp) {
    // Skip first revolution - we need data from two revolutions to compare
    if (prev_temp == 0.0) {
        return false;
    }
    
    // Calculate temperature difference between current and previous revolution
    double temp_diff = fabs(stats.T_max - prev_temp);
    return stats.valid && temp_diff <= TEMPERATURE_CONVERGENCE_THRESHOLD;
};

void log_thermal_data(const string& label, const string& filepath, int rev_num, ostream& out) {
    ThermalStats stats = parse_vtk(filepath);
    if (!stats.valid) return;

    out << label << " Revolution " << rev_num << " complete.\n";
    out << "Starting Thermal Analysis:\n";
    out << "T_" << label << "_Min:   " << fixed << setprecision(1) << stats.T_min << " C\n";
    out << "T_" << label << "_Mean:  " << stats.T_mean << " C\n";
    out << "T_" << label << "_Max:   " << stats.T_max << " C\n\n";
    out << "Thermal Deformation after Rev " << rev_num << ":\n";
    out << "DT_" << label << ":min:  " << stats.D_min << " mu\n";
    out << "DT_" << label << ":max:  " << stats.D_max << " mu\n\n";
}

// NEW TEMPERATURE FILE LOGGING FUNCTIONS
// Function to initialize temperature files
void initializeTemperatureFiles(bool piston_enabled, bool block_enabled, bool slipper_enabled) {
    // Create temperature output directory
    if (!directoryExists("output/temperature")) {
        if (_mkdir("output/temperature") != 0) {
            // Handle error if needed - could log to main logfile
        }
    }
    
    if (piston_enabled) {
        string piston_temp_file = "output/temperature/piston_temperatures.txt";
        ofstream tempFile(piston_temp_file, ios::trunc);
        if (tempFile.is_open()) {
            tempFile << "=== PISTON INTERFACE TEMPERATURE DATA ===" << endl;
            tempFile << "Generated: " << gt() << endl;
            tempFile << "Temperature data for Piston and Cylinder components" << endl;
            tempFile << "All temperatures in Celsius (°C)" << endl;
            tempFile << "=================================================" << endl;
            tempFile.close();
        }
    }
    
    if (block_enabled) {
        string block_temp_file = "output/temperature/block_temperatures.txt";
        ofstream tempFile(block_temp_file, ios::trunc);
        if (tempFile.is_open()) {
            tempFile << "=== BLOCK INTERFACE TEMPERATURE DATA ===" << endl;
            tempFile << "Generated: " << gt() << endl;
            tempFile << "Temperature data for Cylinder Block and Valve Plate components" << endl;
            tempFile << "All temperatures in Celsius (°C)" << endl;
            tempFile << "=================================================" << endl;
            tempFile.close();
        }
    }
    
    if (slipper_enabled) {
        string slipper_temp_file = "output/temperature/slipper_temperatures.txt";
        ofstream tempFile(slipper_temp_file, ios::trunc);
        if (tempFile.is_open()) {
            tempFile << "=== SLIPPER INTERFACE TEMPERATURE DATA ===" << endl;
            tempFile << "Generated: " << gt() << endl;
            tempFile << "Temperature data for Slipper and Swashplate components" << endl;
            tempFile << "All temperatures in Celsius (°C)" << endl;
            tempFile << "=================================================" << endl;
            tempFile.close();
        }
    }
}

// Enhanced thermal logging function with better file detection and debugging
void thermal_log_per_revolution_with_files(int rev, logger& logfile,
                                          bool piston_enabled, bool block_enabled, bool slipper_enabled) {
    ostringstream thermal_summary;
    
    // Flags to track if we need to write headers for each component
    static bool piston_header_written = false;
    static bool cylinder_header_written = false;
    static bool block_header_written = false;
    static bool valve_header_written = false;
    static bool slipper_header_written = false;
    static bool swash_header_written = false;

    if (piston_enabled) {
        ostringstream ss_piston; ss_piston << rev;
        string piston_vtk = "output/piston/vtk/piston_th_rev_" + ss_piston.str() + ".vtk";
        ostringstream ss_cylinder; ss_cylinder << rev;
        string cylinder_vtk = "output/piston/vtk/cylinder_th_rev_" + ss_cylinder.str() + ".vtk";
        
        // Log to console/logfile
        log_thermal_data("Piston", piston_vtk, rev, thermal_summary);
        log_thermal_data("Cylinder", cylinder_vtk, rev, thermal_summary);
        
        // Write to temperature files
        string piston_temp_file = "output/temperature/piston_temperatures.txt";
        ThermalStats piston_stats = parse_vtk(piston_vtk);
        ThermalStats cylinder_stats = parse_vtk(cylinder_vtk);
        
        if (piston_stats.valid) {
            ofstream tempFile(piston_temp_file, ios::app);
            if (tempFile.is_open()) {
                if (!piston_header_written) {
                    tempFile << endl << "PISTON:" << endl;
                    tempFile << "Revolution\tMax (°C)\tMin (°C)\tMean (°C)" << endl;
                    tempFile << "----------\t--------\t--------\t---------" << endl;
                    piston_header_written = true;
                }
                tempFile << fixed << setprecision(1);
                tempFile << rev + 1 << "\t\t" << piston_stats.T_max << "\t\t" 
                         << piston_stats.T_min << "\t\t" << piston_stats.T_mean << endl;
                tempFile.close();
            }
        }
        
        if (cylinder_stats.valid) {
            ofstream tempFile(piston_temp_file, ios::app);
            if (tempFile.is_open()) {
                if (!cylinder_header_written) {
                    tempFile << endl << "CYLINDER:" << endl;
                    tempFile << "Revolution\tMax (°C)\tMin (°C)\tMean (°C)" << endl;
                    tempFile << "----------\t--------\t--------\t---------" << endl;
                    cylinder_header_written = true;
                }
                tempFile << fixed << setprecision(1);
                tempFile << rev + 1 << "\t\t" << cylinder_stats.T_max << "\t\t" 
                         << cylinder_stats.T_min << "\t\t" << cylinder_stats.T_mean << endl;
                tempFile.close();
            }
        }
    }
    
    if (block_enabled) {
        ostringstream ss_block; ss_block << rev;
        string block_vtk = "output/block/rev."+ ss_block.str()+"/cylinderblock.vtk";
        ostringstream ss_valve; ss_valve << rev;
        string valve_vtk = "output/block/rev."+ ss_valve.str() +"/valveplate.vtk";
        
        // Log to console/logfile
        log_thermal_data("Block", block_vtk, rev, thermal_summary);
        log_thermal_data("Valve", valve_vtk, rev, thermal_summary);
        
        // Write to temperature files
        string block_temp_file = "output/temperature/block_temperatures.txt";
        ThermalStats block_stats = parse_vtk(block_vtk);
        ThermalStats valve_stats = parse_vtk(valve_vtk);
        
        if (block_stats.valid) {
            ofstream tempFile(block_temp_file, ios::app);
            if (tempFile.is_open()) {
                if (!block_header_written) {
                    tempFile << endl << "CYLINDER BLOCK:" << endl;
                    tempFile << "Revolution\tMax (°C)\tMin (°C)\tMean (°C)" << endl;
                    tempFile << "----------\t--------\t--------\t---------" << endl;
                    block_header_written = true;
                }
                tempFile << fixed << setprecision(1);
                tempFile << rev + 1 << "\t\t" << block_stats.T_max << "\t\t" 
                         << block_stats.T_min << "\t\t" << block_stats.T_mean << endl;
                tempFile.close();
            }
        }
        
        if (valve_stats.valid) {
            ofstream tempFile(block_temp_file, ios::app);
            if (tempFile.is_open()) {
                if (!valve_header_written) {
                    tempFile << endl << "VALVE PLATE:" << endl;
                    tempFile << "Revolution\tMax (°C)\tMin (°C)\tMean (°C)" << endl;
                    tempFile << "----------\t--------\t--------\t---------" << endl;
                    valve_header_written = true;
                }
                tempFile << fixed << setprecision(1);
                tempFile << rev + 1 << "\t\t" << valve_stats.T_max << "\t\t" 
                         << valve_stats.T_min << "\t\t" << valve_stats.T_mean << endl;
                tempFile.close();
            }
        }
    }
    
    if (slipper_enabled) {
        // Enhanced slipper file detection with multiple possible naming patterns
        vector<string> slipper_patterns;
        vector<string> swash_patterns;
        
        // Try different possible naming patterns
        int phi = 360 * rev;
        ostringstream ss_slipper; ss_slipper << phi;
        ostringstream ss_rev; ss_rev << rev;
        
        // Pattern 1: Correct format based on actual file names (slipper_thermal.0.vtk)
        //slipper_patterns.push_back("output/slipper/vtk/slipper_thermal." + ss_rev.str() + ".vtk");
        //swash_patterns.push_back("output/slipper/vtk/swashplate_thermal." + ss_rev.str() + ".vtk");
        
        // Pattern 2: Original pattern (phi-based naming)
        slipper_patterns.push_back("output/slipper/vtk/slipper_thermal." + ss_slipper.str() + ".vtk");
        swash_patterns.push_back("output/slipper/vtk/swashplate_thermal." + ss_slipper.str() + ".vtk");
        
        // Pattern 3: Rev-based naming with underscore
        //slipper_patterns.push_back("output/slipper/vtk/slipper_thermal_rev_" + ss_rev.str() + ".vtk");
        //swash_patterns.push_back("output/slipper/vtk/swashplate_thermal_rev_" + ss_rev.str() + ".vtk");
        
        // Pattern 4: Different naming convention
        //slipper_patterns.push_back("output/slipper/vtk/slipper_th_rev_" + ss_rev.str() + ".vtk");
        //swash_patterns.push_back("output/slipper/vtk/swashplate_th_rev_" + ss_rev.str() + ".vtk");
        
        // Pattern 5: Simple naming
        slipper_patterns.push_back("output/slipper/slipper_thermal.vtk");
        swash_patterns.push_back("output/slipper/swashplate_thermal.vtk");
        
        // Find which file exists
        string slipper_vtk = "";
        string swash_vtk = "";
        
        for (size_t i = 0; i < slipper_patterns.size(); i++) {
            if (fileExists(slipper_patterns[i])) {
                slipper_vtk = slipper_patterns[i];
                break;
            }
        }
        
        for (size_t i = 0; i < swash_patterns.size(); i++) {
            if (fileExists(swash_patterns[i])) {
                swash_vtk = swash_patterns[i];
                break;
            }
        }
        
        // Debug logging to help identify the issue
        if (slipper_vtk.empty()) {
            logfile << gt() << "DEBUG: No slipper thermal file found. Tried patterns:" << endl;
            for (size_t i = 0; i < slipper_patterns.size(); i++) {
                logfile << gt() << "  -> " << slipper_patterns[i] << " (exists: " << (fileExists(slipper_patterns[i]) ? "YES" : "NO") << ")" << endl;
            }
        }
        
        if (swash_vtk.empty()) {
            logfile << gt() << "DEBUG: No swashplate thermal file found. Tried patterns:" << endl;
            for (size_t i = 0; i < swash_patterns.size(); i++) {
                logfile << gt() << "  -> " << swash_patterns[i] << " (exists: " << (fileExists(swash_patterns[i]) ? "YES" : "NO") << ")" << endl;
            }
        }
        
        // Process files if found
        if (!slipper_vtk.empty()) {
            log_thermal_data("Slipper", slipper_vtk, rev, thermal_summary);
        }
        if (!swash_vtk.empty()) {
            log_thermal_data("Swash", swash_vtk, rev, thermal_summary);
        }
        
        // Write to temperature files
        string slipper_temp_file = "output/temperature/slipper_temperatures.txt";
        ThermalStats slipper_stats, swash_stats;
        
        if (!slipper_vtk.empty()) {
            slipper_stats = parse_vtk(slipper_vtk);
        }
        if (!swash_vtk.empty()) {
            swash_stats = parse_vtk(swash_vtk);
        }
        
        if (slipper_stats.valid) {
            ofstream tempFile(slipper_temp_file, ios::app);
            if (tempFile.is_open()) {
                if (!slipper_header_written) {
                    tempFile << endl << "SLIPPER:" << endl;
                    tempFile << "Revolution\tMax (°C)\tMin (°C)\tMean (°C)" << endl;
                    tempFile << "----------\t--------\t--------\t---------" << endl;
                    slipper_header_written = true;
                }
                tempFile << fixed << setprecision(1);
                tempFile << rev + 1 << "\t\t" << slipper_stats.T_max << "\t\t" 
                         << slipper_stats.T_min << "\t\t" << slipper_stats.T_mean << endl;
                tempFile.close();
            }
        } else if (!slipper_vtk.empty()) {
            // File exists but parsing failed
            logfile << gt() << "WARNING: Slipper VTK file exists but parsing failed: " << slipper_vtk << endl;
        }
        
        if (swash_stats.valid) {
            ofstream tempFile(slipper_temp_file, ios::app);
            if (tempFile.is_open()) {
                if (!swash_header_written) {
                    tempFile << endl << "SWASHPLATE:" << endl;
                    tempFile << "Revolution\tMax (°C)\tMin (°C)\tMean (°C)" << endl;
                    tempFile << "----------\t--------\t--------\t---------" << endl;
                    swash_header_written = true;
                }
                tempFile << fixed << setprecision(1);
                tempFile << rev + 1 << "\t\t" << swash_stats.T_max << "\t\t" 
                         << swash_stats.T_min << "\t\t" << swash_stats.T_mean << endl;
                tempFile.close();
            }
        } else if (!swash_vtk.empty()) {
            // File exists but parsing failed
            logfile << gt() << "WARNING: Swashplate VTK file exists but parsing failed: " << swash_vtk << endl;
        }
        
        // If no files were found at all, log this information
        if (slipper_vtk.empty() && swash_vtk.empty()) {
            logfile << gt() << "WARNING: No slipper thermal VTK files found for revolution " << rev << endl;
            logfile << gt() << "  -> Check if slipper thermal analysis is enabled and generating output files" << endl;
        }
    }

    logfile << thermal_summary.str();
    cout << thermal_summary.str();
}


void thermal_log_per_revolution(int rev, logger& logfile,
                            bool piston_enabled, bool block_enabled, bool slipper_enabled) {
    ostringstream thermal_summary;

    if (piston_enabled) {
        ostringstream ss_piston; ss_piston << rev;
        string piston_vtk = "output/piston/vtk/piston_th_rev_" + ss_piston.str() + ".vtk";
        ostringstream ss_cylinder; ss_cylinder << rev;
        string cylinder_vtk = "output/piston/vtk/cylinder_th_rev_" + ss_cylinder.str() + ".vtk";
        log_thermal_data("Piston", piston_vtk, rev, thermal_summary);
        log_thermal_data("Cylinder", cylinder_vtk, rev, thermal_summary);
    }
    if (block_enabled) {
        ostringstream ss_block; ss_block << rev;
        string block_vtk = "output/block/rev."+ ss_block.str()+"/cylinderblock.vtk";
        ostringstream ss_valve; ss_valve << rev;
        string valve_vtk = "output/block/rev."+ ss_valve.str() +"/valveplate.vtk";
        log_thermal_data("Block", block_vtk, rev, thermal_summary);
        log_thermal_data("Valve", valve_vtk, rev, thermal_summary);
    }
    if (slipper_enabled) {
        int phi = 360 * rev;
        ostringstream ss_slipper; ss_slipper << phi;
        string slipper_vtk = "output/slipper/vtk/slipper_thermal." + ss_slipper.str() + ".vtk";
        ostringstream ss_swash; ss_swash << phi;
        string swash_vtk = "output/slipper/vtk/swashplate_thermal." + ss_swash.str() + ".vtk";
        log_thermal_data("Slipper", slipper_vtk, rev, thermal_summary);
        log_thermal_data("Swash", swash_vtk, rev, thermal_summary);
    }

    logfile << thermal_summary.str();
    cout << thermal_summary.str();
}

// This function checks temperature convergence between consecutive revolutions
// Returns true if any interface was disabled due to convergence
bool check_temperature_convergence(int rev, vector<bool>& interface_enabled, vector<bool>& interface_disabled, 
                                int& active_interfaces, logger& logfile) {
    // Skip checking first revolution - need at least two revolutions to compare
    if (rev <= 0) return false;
    
    bool any_disabled = false;
    
    if (interface_enabled[0] && !interface_disabled[0]) {
        ostringstream ss_piston; ss_piston << rev;
        string piston_vtk = "output/piston/vtk/piston_th_rev_" + ss_piston.str() + ".vtk";
        ThermalStats piston_stats = parse_vtk(piston_vtk);
        
        if (piston_stats.valid) {
            // Check if this is not the first revolution and if converged
            if (prev_max_temp[0] > 0.0) {
                double temp_diff = fabs(piston_stats.T_max - prev_max_temp[0]);
                logfile << gt() << "Piston interface temperature difference: " 
                       << fixed << setprecision(2) << temp_diff << " C (threshold: " 
                       << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
                
                if (temp_diff <= TEMPERATURE_CONVERGENCE_THRESHOLD) {
                    logfile << gt() << "Piston interface temperature converged (diff: " 
                           << temp_diff << " C <= " << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
                    logfile << gt() << "Disabling Piston interface for future revolutions." << endl;
                    interface_disabled[0] = true;
                    active_interfaces--;
                    any_disabled = true;
                }
            }
            // Update previous temperature for next comparison
            prev_max_temp[0] = piston_stats.T_max;
        }
    }

    if (interface_enabled[1] && !interface_disabled[1]) {
        ostringstream ss_block; ss_block << rev;
        string block_vtk = "output/block/rev."+ ss_block.str()+"/cylinderblock.vtk";
        ThermalStats block_stats = parse_vtk(block_vtk);
        
        if (block_stats.valid) {
            // Check if this is not the first revolution and if converged
            if (prev_max_temp[1] > 0.0) {
                double temp_diff = fabs(block_stats.T_max - prev_max_temp[1]);
                logfile << gt() << "Block interface temperature difference: " 
                       << fixed << setprecision(2) << temp_diff << " C (threshold: " 
                       << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
                
                if (temp_diff <= TEMPERATURE_CONVERGENCE_THRESHOLD) {
                    logfile << gt() << "Block interface temperature converged (diff: " 
                           << temp_diff << " C <= " << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
                    logfile << gt() << "Disabling Block interface for future revolutions." << endl;
                    interface_disabled[1] = true;
                    active_interfaces--;
                    any_disabled = true;
                }
            }
            // Update previous temperature for next comparison
            prev_max_temp[1] = block_stats.T_max;
        }
    }

    if (interface_enabled[2] && !interface_disabled[2]) {
        int phi = 360 * rev;
        ostringstream ss_slipper; ss_slipper << phi;
        string slipper_vtk = "output/slipper/vtk/slipper_thermal." + ss_slipper.str() + ".vtk";
        ThermalStats slipper_stats = parse_vtk(slipper_vtk);
        
        if (slipper_stats.valid) {
            // Check if this is not the first revolution and if converged
            if (prev_max_temp[2] > 0.0) {
                double temp_diff = fabs(slipper_stats.T_max - prev_max_temp[2]);
                logfile << gt() << "Slipper interface temperature difference: " 
                       << fixed << setprecision(2) << temp_diff << " C (threshold: " 
                       << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
                
                if (temp_diff <= TEMPERATURE_CONVERGENCE_THRESHOLD) {
                    logfile << gt() << "Slipper interface temperature converged (diff: " 
                           << temp_diff << " C <= " << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
                    logfile << gt() << "Disabling Slipper interface for future revolutions." << endl;
                    interface_disabled[2] = true;
                    active_interfaces--;
                    any_disabled = true;
                }
            }
            // Update previous temperature for next comparison
            prev_max_temp[2] = slipper_stats.T_max;
        }
    }
    
    return any_disabled;
}

// Add this new function before the main function
void generatePumpSummaryReport(int finalRevolution, logger& logfile,
                              bool piston_enabled, bool block_enabled, bool slipper_enabled,
                              bool piston_disabled, bool block_disabled, bool slipper_disabled) {
    
    logfile << gt() << "Generating comprehensive pump summary report..." << std::endl;
    
    // Create pump.txt file
    std::ofstream pumpFile("output/pump.txt", std::ios::trunc);
    if (!pumpFile.is_open()) {
        logfile << gt() << "ERROR: Could not create pump.txt file." << std::endl;
        return;
    }
    
    // Get operating conditions
    input gapinput;
    
    // File header
    pumpFile << "========================================================================" << std::endl;
    pumpFile << "                    COMPREHENSIVE PUMP PERFORMANCE REPORT" << std::endl;
    pumpFile << "========================================================================" << std::endl;
    pumpFile << "Generated: " << gt() << std::endl;
    pumpFile << "Final Revolution: " << (finalRevolution + 1) << std::endl;
    pumpFile << "========================================================================" << std::endl;
    pumpFile << std::endl;
    
    // Operating Conditions Summary
    double speed = gapinput.data.operating_conditions.speed * 60/(2*pi);
    double dp = (gapinput.data.operating_conditions.HP - gapinput.data.operating_conditions.LP)/1e5;
    double betaNorm = gapinput.data.operating_conditions.beta / gapinput.data.operating_conditions.betamax;
    double V = 110.0; // Displacement volume
    double flow_theo = speed * V * betaNorm / 1000.0;  // Theoretical flow in l/min
    
    pumpFile << "OPERATING CONDITIONS:" << std::endl;
    pumpFile << "--------------------" << std::endl;
    pumpFile << std::fixed << std::setprecision(2);
    pumpFile << "Speed            : " << speed << " rpm" << std::endl;
    pumpFile << "High Pressure    : " << gapinput.data.operating_conditions.HP/1e5 << " bar" << std::endl;
    pumpFile << "Low Pressure     : " << gapinput.data.operating_conditions.LP/1e5 << " bar" << std::endl;
    pumpFile << "Pressure Drop    : " << dp << " bar" << std::endl;
    pumpFile << "Beta Angle       : " << gapinput.data.operating_conditions.beta*180/pi << " degrees" << std::endl;
    pumpFile << "Max Beta Angle   : " << gapinput.data.operating_conditions.betamax*180/pi << " degrees" << std::endl;
    pumpFile << "Theoretical Flow : " << flow_theo << " l/min" << std::endl;
    pumpFile << std::endl;
    
    // Interface Status
    pumpFile << "INTERFACE STATUS:" << std::endl;
    pumpFile << "-----------------" << std::endl;
    pumpFile << "Piston Interface : " << (piston_enabled ? (piston_disabled ? "Converged" : "Active") : "Disabled") << std::endl;
    pumpFile << "Block Interface  : " << (block_enabled ? (block_disabled ? "Converged" : "Active") : "Disabled") << std::endl;
    pumpFile << "Slipper Interface: " << (slipper_enabled ? (slipper_disabled ? "Converged" : "Active") : "Disabled") << std::endl;
    pumpFile << std::endl;
    
    // Variables to store summary data
    double total_leakage = 0.0;
    double total_mechanical_loss = 0.0;
    double total_volumetric_loss = 0.0;
    double weighted_gap_height = 0.0;
    double max_temperature = 0.0;
    double mean_temperature = 0.0;
    int active_interface_count = 0;
    
    // ===========================================
    // PISTON INTERFACE ANALYSIS
    // ===========================================
    if (piston_enabled) {
        pumpFile << "========================================================================" << std::endl;
        pumpFile << "                           PISTON INTERFACE" << std::endl;
        pumpFile << "========================================================================" << std::endl;
        
        std::ifstream pistonFile("output/piston/piston.txt");
        if (pistonFile.is_open()) {
            std::vector<std::vector<double>> revData;
            std::string line;
            
            // Read data for the final revolution
            while (std::getline(pistonFile, line)) {
                std::istringstream iss(line);
                std::vector<double> values;
                double value;
                
                while (iss >> value) {
                    values.push_back(value);
                }
                
                if (values.size() >= 62) {
                    double lineRev = values[1];
                    if (floor(lineRev) == finalRevolution) {
                        revData.push_back(values);
                    }
                }
            }
            pistonFile.close();
            
            if (!revData.empty()) {
                active_interface_count++;
                
                // Calculate averages
                double avg_gap = 0.0, avg_leakage = 0.0, avg_mech_loss = 0.0, avg_vol_loss = 0.0;
                for (size_t idx = 0; idx < revData.size(); ++idx) {
					const std::vector<double>& row = revData[idx];
					avg_gap += row[11];        // Gap height
					avg_leakage += row[59];    // Total leakage
					avg_mech_loss += row[61];  // Mechanical power loss
					avg_vol_loss += row[62];   // Volumetric power loss
				}
                avg_gap /= revData.size();
                avg_leakage /= revData.size();
                avg_mech_loss /= revData.size();
                avg_vol_loss /= revData.size();
                
                // Convert leakage to l/min
                avg_leakage *= 60000.0;
                
                // Calculate efficiency metrics
                double leakage = avg_leakage;
                double PlossF = avg_mech_loss;
                double PlossL = avg_vol_loss;
                double Ploss = PlossF + PlossL;
                double flow = flow_theo - leakage;
                double nu_vol = flow_theo > 0 ? flow / flow_theo : 0;
                double P_theo = dp * flow_theo / 600.0;
                double P_actual = P_theo - Ploss / 1000.0;
                double nu_ges = P_theo > 0 ? P_actual / P_theo : 0;
                double nu_hm = nu_vol > 0 ? nu_ges / nu_vol : 0;
                
                // Get thermal data
                std::ostringstream ss_piston; ss_piston << finalRevolution;
                std::string piston_vtk = "output/piston/vtk/piston_th_rev_" + ss_piston.str() + ".vtk";
                ThermalStats piston_thermal = parse_vtk(piston_vtk);
                
                // Output piston data
                pumpFile << "Performance Metrics:" << std::endl;
                pumpFile << "-------------------" << std::endl;
                pumpFile << std::fixed << std::setprecision(3);
                pumpFile << "Mean Gap Height         : " << avg_gap << " µm" << std::endl;
                pumpFile << "Leakage Flow            : " << leakage << " l/min" << std::endl;
                pumpFile << "Net Flow                : " << flow << " l/min" << std::endl;
                pumpFile << "Mechanical Power Loss   : " << PlossF << " W" << std::endl;
                pumpFile << "Volumetric Power Loss   : " << PlossL << " W" << std::endl;
                pumpFile << "Total Power Loss        : " << Ploss << " W" << std::endl;
                pumpFile << std::fixed << std::setprecision(2);
                pumpFile << "Volumetric Efficiency   : " << (nu_vol * 100) << " %" << std::endl;
                pumpFile << "Hydro-Mechanical Eff.   : " << (nu_hm * 100) << " %" << std::endl;
                pumpFile << "Total Efficiency        : " << (nu_ges * 100) << " %" << std::endl;
                
                if (piston_thermal.valid) {
                    pumpFile << std::endl << "Thermal Analysis:" << std::endl;
                    pumpFile << "----------------" << std::endl;
                    pumpFile << std::fixed << std::setprecision(1);
                    pumpFile << "Min Temperature         : " << piston_thermal.T_min << " °C" << std::endl;
                    pumpFile << "Mean Temperature        : " << piston_thermal.T_mean << " °C" << std::endl;
                    pumpFile << "Max Temperature         : " << piston_thermal.T_max << " °C" << std::endl;
                    pumpFile << "Min Deformation         : " << piston_thermal.D_min << " mu" << std::endl;
                    pumpFile << "Max Deformation         : " << piston_thermal.D_max << " mu" << std::endl;
                    
                  // NEW CODE (C++98 compatible):
					if (piston_thermal.T_max > max_temperature) {
						max_temperature = piston_thermal.T_max;
					}
                    mean_temperature += piston_thermal.T_mean;
                }
                
                pumpFile << std::endl;
                
                // Update totals
                total_leakage += leakage;
                total_mechanical_loss += PlossF;
                total_volumetric_loss += PlossL;
                weighted_gap_height += avg_gap;
            }
        } else {
            pumpFile << "ERROR: Could not read piston data file." << std::endl;
        }
        pumpFile << std::endl;
    }
    
    // ===========================================
    // BLOCK INTERFACE ANALYSIS
    // ===========================================
    if (block_enabled) {
        pumpFile << "========================================================================" << std::endl;
        pumpFile << "                            BLOCK INTERFACE" << std::endl;
        pumpFile << "========================================================================" << std::endl;
        
        std::ifstream blockFile("output/block/block.txt");
        if (blockFile.is_open()) {
            std::vector<std::vector<double>> revData;
            std::string line;
            
            // Read data for the final revolution
            while (std::getline(blockFile, line)) {
                if (line.empty() || line[0] == '%' || line[0] == 't') continue;
                
                std::istringstream iss(line);
                std::vector<double> values;
                double value;
                
                while (iss >> value) {
                    values.push_back(value);
                }
                
                if (values.size() >= 25) {
                    double lineRev = values[1];
                    if (floor(lineRev) == finalRevolution) {
                        revData.push_back(values);
                    }
                }
            }
            blockFile.close();
            
            if (!revData.empty()) {
                active_interface_count++;
                
                // Calculate averages
                double avg_hmin = 0.0, avg_hmean = 0.0, avg_hmax = 0.0;
                double avg_leakage = 0.0, avg_mech_loss = 0.0, avg_vol_loss = 0.0;
                
                for (size_t idx = 0; idx < revData.size(); ++idx) {
					const std::vector<double>& row = revData[idx];
					avg_hmin += row[12];      // hmin
					avg_hmean += row[13];     // hmean  
					avg_hmax += row[14];      // hmax
					avg_leakage += row[21];   // Leak
					avg_mech_loss += row[24]; // PlossF
					avg_vol_loss += row[23];  // PlossL
				}
                avg_hmin /= revData.size();
                avg_hmean /= revData.size();
                avg_hmax /= revData.size();
                avg_leakage /= revData.size();
                avg_mech_loss /= revData.size();
                avg_vol_loss /= revData.size();
                
                // Calculate efficiency metrics
                double leakage = avg_leakage;
                double PlossF = avg_mech_loss;
                double PlossL = avg_vol_loss;
                double Ploss = PlossF + PlossL;
                double flow = flow_theo - leakage;
                double nu_vol = flow_theo > 0 ? flow / flow_theo : 0;
                double P_theo = dp * flow_theo / 600.0;
                double P_actual = P_theo - Ploss / 1000.0;
                double nu_ges = P_theo > 0 ? P_actual / P_theo : 0;
                double nu_hm = nu_vol > 0 ? nu_ges / nu_vol : 0;
                
                // Get thermal data
                std::ostringstream ss_block; ss_block << finalRevolution;
                std::string block_vtk = "output/block/rev." + ss_block.str() + "/cylinderblock.vtk";
                ThermalStats block_thermal = parse_vtk(block_vtk);
                
                // Output block data
                pumpFile << "Performance Metrics:" << std::endl;
                pumpFile << "-------------------" << std::endl;
                pumpFile << std::fixed << std::setprecision(3);
                pumpFile << "Min Gap Height          : " << avg_hmin << " µm" << std::endl;
                pumpFile << "Mean Gap Height         : " << avg_hmean << " µm" << std::endl;
                pumpFile << "Max Gap Height          : " << avg_hmax << " µm" << std::endl;
                pumpFile << "Leakage Flow            : " << leakage << " l/min" << std::endl;
                pumpFile << "Net Flow                : " << flow << " l/min" << std::endl;
                pumpFile << "Mechanical Power Loss   : " << PlossF << " W" << std::endl;
                pumpFile << "Volumetric Power Loss   : " << PlossL << " W" << std::endl;
                pumpFile << "Total Power Loss        : " << Ploss << " W" << std::endl;
                pumpFile << std::fixed << std::setprecision(2);
                pumpFile << "Volumetric Efficiency   : " << (nu_vol * 100) << " %" << std::endl;
                pumpFile << "Hydro-Mechanical Eff.   : " << (nu_hm * 100) << " %" << std::endl;
                pumpFile << "Total Efficiency        : " << (nu_ges * 100) << " %" << std::endl;
                
                if (block_thermal.valid) {
                    pumpFile << std::endl << "Thermal Analysis:" << std::endl;
                    pumpFile << "----------------" << std::endl;
                    pumpFile << std::fixed << std::setprecision(1);
                    pumpFile << "Min Temperature         : " << block_thermal.T_min << " °C" << std::endl;
                    pumpFile << "Mean Temperature        : " << block_thermal.T_mean << " °C" << std::endl;
                    pumpFile << "Max Temperature         : " << block_thermal.T_max << " °C" << std::endl;
                    pumpFile << "Min Deformation         : " << block_thermal.D_min << " mu" << std::endl;
                    pumpFile << "Max Deformation         : " << block_thermal.D_max << " mu" << std::endl;
                    
                    // Update summary variables
                    if (block_thermal.T_max > max_temperature) {
						max_temperature = block_thermal.T_max;
					}
                    mean_temperature += block_thermal.T_mean;
                }
                
                pumpFile << std::endl;
                
                // Update totals
                total_leakage += leakage;
                total_mechanical_loss += PlossF;
                total_volumetric_loss += PlossL;
                weighted_gap_height += avg_hmean; // Use mean gap height for block
            }
        } else {
            pumpFile << "ERROR: Could not read block data file." << std::endl;
        }
        pumpFile << std::endl;
    }
    
    // ===========================================
    // SLIPPER INTERFACE ANALYSIS
    // ===========================================
    if (slipper_enabled) {
        pumpFile << "========================================================================" << std::endl;
        pumpFile << "                          SLIPPER INTERFACE" << std::endl;
        pumpFile << "========================================================================" << std::endl;
        
        std::ifstream slipperFile("output/slipper/slipper.txt");
        if (slipperFile.is_open()) {
            std::vector<std::vector<double>> revData;
            std::string line;
            
            // Read data for the final revolution
            while (std::getline(slipperFile, line)) {
                if (line.empty() || line[0] == '%') continue;
                
                std::istringstream iss(line);
                std::vector<double> values;
                double value;
                
                while (iss >> value) {
                    values.push_back(value);
                }
                
                if (values.size() >= 18) {
                    double lineRev = values[1];
                    if (floor(lineRev) == finalRevolution) {
                        revData.push_back(values);
                    }
                }
            }
            slipperFile.close();
            
            if (!revData.empty()) {
                active_interface_count++;
                
                // Calculate averages
                double avg_hmin = 0.0, avg_hmean = 0.0, avg_hmax = 0.0;
                double avg_leakage = 0.0, avg_torque_loss = 0.0, avg_power_loss = 0.0;
                
                for (size_t idx = 0; idx < revData.size(); ++idx) {
					const std::vector<double>& row = revData[idx];
					avg_hmin += row[9];       // hmin
					avg_hmean += row[10];     // hmean
					avg_hmax += row[11];      // hmax
					avg_leakage += row[12];   // QSG (leakage)
					avg_torque_loss += row[13]; // TorqueLoss
					avg_power_loss += row[14];  // PowerLoss
				}
                avg_hmin /= revData.size();
                avg_hmean /= revData.size();
                avg_hmax /= revData.size();
                avg_leakage /= revData.size();
                avg_torque_loss /= revData.size();
                avg_power_loss /= revData.size();
                
                // Calculate efficiency metrics
                double leakage = avg_leakage;
                double torqueLoss = avg_torque_loss;
                double powerLoss = avg_power_loss;
                double flow = flow_theo - leakage;
                double nu_vol = flow_theo > 0 ? flow / flow_theo : 0;
                double P_theo = dp * flow_theo / 600.0;
                double P_actual = P_theo - powerLoss / 1000.0;
                double nu_ges = P_theo > 0 ? P_actual / P_theo : 0;
                double nu_hm = nu_vol > 0 ? nu_ges / nu_vol : 0;
                
                // Get thermal data
                int phi = 360 * finalRevolution;
                std::ostringstream ss_slipper; ss_slipper << phi;
                std::string slipper_vtk = "output/slipper/vtk/slipper_thermal." + ss_slipper.str() + ".vtk";
                ThermalStats slipper_thermal = parse_vtk(slipper_vtk);
                
                // Output slipper data
                pumpFile << "Performance Metrics:" << std::endl;
                pumpFile << "-------------------" << std::endl;
                pumpFile << std::fixed << std::setprecision(3);
                pumpFile << "Min Gap Height          : " << avg_hmin << " µm" << std::endl;
                pumpFile << "Mean Gap Height         : " << avg_hmean << " µm" << std::endl;
                pumpFile << "Max Gap Height          : " << avg_hmax << " µm" << std::endl;
                pumpFile << "Leakage Flow            : " << leakage << " l/min" << std::endl;
                pumpFile << "Net Flow                : " << flow << " l/min" << std::endl;
                pumpFile << "Torque Loss             : " << torqueLoss << " Nm" << std::endl;
                pumpFile << "Power Loss              : " << powerLoss << " W" << std::endl;
                pumpFile << std::fixed << std::setprecision(2);
                pumpFile << "Volumetric Efficiency   : " << (nu_vol * 100) << " %" << std::endl;
                pumpFile << "Hydro-Mechanical Eff.   : " << (nu_hm * 100) << " %" << std::endl;
                pumpFile << "Total Efficiency        : " << (nu_ges * 100) << " %" << std::endl;
                
                if (slipper_thermal.valid) {
                    pumpFile << std::endl << "Thermal Analysis:" << std::endl;
                    pumpFile << "----------------" << std::endl;
                    pumpFile << std::fixed << std::setprecision(1);
                    pumpFile << "Min Temperature         : " << slipper_thermal.T_min << " °C" << std::endl;
                    pumpFile << "Mean Temperature        : " << slipper_thermal.T_mean << " °C" << std::endl;
                    pumpFile << "Max Temperature         : " << slipper_thermal.T_max << " °C" << std::endl;
                    pumpFile << "Min Deformation         : " << slipper_thermal.D_min << " mu" << std::endl;
                    pumpFile << "Max Deformation         : " << slipper_thermal.D_max << " mu" << std::endl;
                    
                    // Update summary variables
                    if (slipper_thermal.T_max > max_temperature) {
						max_temperature = slipper_thermal.T_max;
					}
                    mean_temperature += slipper_thermal.T_mean;
                }
                
                pumpFile << std::endl;
                
                // Update totals
                total_leakage += leakage;
                total_mechanical_loss += powerLoss; // For slipper, power loss is the mechanical loss
                weighted_gap_height += avg_hmean; // Use mean gap height for slipper
            }
        } else {
            pumpFile << "ERROR: Could not read slipper data file." << std::endl;
        }
        pumpFile << std::endl;
    }
    
    // ===========================================
    // OVERALL PUMP SUMMARY
    // ===========================================
    if (active_interface_count > 0) {
        pumpFile << "========================================================================" << std::endl;
        pumpFile << "                          OVERALL PUMP SUMMARY" << std::endl;
        pumpFile << "========================================================================" << std::endl;
        
        // Calculate overall metrics
        double net_flow = flow_theo - total_leakage;
        double total_power_loss = total_mechanical_loss + total_volumetric_loss;
        double overall_vol_eff = flow_theo > 0 ? net_flow / flow_theo : 0;
        double P_theo_total = dp * flow_theo / 600.0;
        double P_actual_total = P_theo_total - total_power_loss / 1000.0;
        double overall_eff = P_theo_total > 0 ? P_actual_total / P_theo_total : 0;
        double overall_hm_eff = overall_vol_eff > 0 ? overall_eff / overall_vol_eff : 0;
        
        // Average temperatures and gap heights
        if (active_interface_count > 0) {
            mean_temperature /= active_interface_count;
            weighted_gap_height /= active_interface_count;
        }
        
        pumpFile << "Combined Performance:" << std::endl;
        pumpFile << "--------------------" << std::endl;
        pumpFile << std::fixed << std::setprecision(3);
        pumpFile << "Active Interfaces       : " << active_interface_count << std::endl;
        pumpFile << "Theoretical Flow        : " << flow_theo << " l/min" << std::endl;
        pumpFile << "Total Leakage           : " << total_leakage << " l/min" << std::endl;
        pumpFile << "Net Flow                : " << net_flow << " l/min" << std::endl;
        pumpFile << "Average Gap Height      : " << weighted_gap_height << " mu" << std::endl;
        pumpFile << "Total Mechanical Loss   : " << total_mechanical_loss << " W" << std::endl;
        pumpFile << "Total Volumetric Loss   : " << total_volumetric_loss << " W" << std::endl;
        pumpFile << "Total Power Loss        : " << total_power_loss << " W" << std::endl;
        pumpFile << std::fixed << std::setprecision(2);
        pumpFile << "Overall Vol. Efficiency : " << (overall_vol_eff * 100) << " %" << std::endl;
        pumpFile << "Overall H-M Efficiency  : " << (overall_hm_eff * 100) << " %" << std::endl;
        pumpFile << "Overall Total Efficiency: " << (overall_eff * 100) << " %" << std::endl;
        
        pumpFile << std::endl << "Thermal Summary:" << std::endl;
        pumpFile << "---------------" << std::endl;
        pumpFile << std::fixed << std::setprecision(1);
        pumpFile << "Average Mean Temperature: " << mean_temperature << " °C" << std::endl;
        pumpFile << "Peak Temperature        : " << max_temperature << " °C" << std::endl;
    }
    
    // Footer
    pumpFile << std::endl;
    pumpFile << "========================================================================" << std::endl;
    pumpFile << "                              END OF REPORT" << std::endl;
    pumpFile << "========================================================================" << std::endl;
    
    pumpFile.close();
    
    logfile << gt() << "Pump summary report generated successfully: output/pump.txt" << std::endl;
}


void atexithandler(void) {
    if (ininputs) {
        logfile << endl << endl;
        logfile << "There was an error reading the input txt files. Below is a copy of the input_log.txt file:" << endl << endl;
        ifstream inputlog("./input_log.txt");
        if (!inputlog.is_open()) {
            logfile << "Unable to open input_log.txt for reading. Please refer directly to this file for debugging information." << endl;
            return;
        }
        string s;
        while (getline(inputlog, s)) {
            logfile << s << endl;
        }
        logfile << endl;
    } else if (!shouldexit) {
        logfile << endl;
        logfile << "//" << endl;
        logfile << "// WARNING: Main program thread caught an expected exit signal." << endl;
        logfile << gt() << "WARNING: Beginning main program exit ..." << endl;
        logfile << gt() << "WARNING: If this exit was not user activated, check all log files for errors." << endl;
        logfile << gt() << "WARNING: Standby for main simulation thread termination ..." << endl;
    }
}

string i2word(const int i) {
    switch (i) {
        case 0: return "Piston";
        case 1: return "Block";
        case 2: return "Slipper";
        default: return "Badi";
    }
}

__declspec(dllexport) int fsti_main(int argc, char* argv[]) {
    SetConsoleTitle("DEBUG: Simulation started – waiting for revolutions...");

    if (argc == 2) {
        if (string(argv[1]) == "-BlockPreCheck_checkGap") {
            logfile << gt() << "* Running the Block Pre-Check of the gap definition ..." << endl;
            char* a[2]; a[1] = "-checkGap";
            return run_block_standalone(2, a);
        }
        if (string(argv[1]) == "-BlockPreCheck_checkEHD") {
            logfile << gt() << "* Running the Block Pre-Check of the influence matrices EHD definition ..." << endl;
            char* a[2]; a[1] = "-checkEHD";
            return run_block_standalone(2, a);
        }
        if (string(argv[1]) == "-BlockPreCheck_checkThermal") {
            logfile << gt() << "* Running the Block Pre-Check of the thermal inputs ..." << endl;
            char* a[2]; a[1] = "-checkThermal";
            return run_block_standalone(2, a);
        }
        if (string(argv[1]) == "-BlockPreCheck_checkOil") {
            logfile << gt() << "* Running the Block Pre-Check of the thermal inputs ..." << endl;
            char* a[2]; a[1] = "-getOilProperties";
            return run_block_standalone(2, a);
        }
    }

    logfile << gt() << "* Starting the main lubrication simulation thread ...";
    logfile.flush();
    atexit(atexithandler);
    logfile << " done." << endl;

    logfile << gt() << "* Loading txt file inputs ...";
    logfile.flush();
    ininputs = true;
    input gapinput;
    ininputs = false;
    logfile << " done." << endl;

    logfile << gt() << "* Preforming other pre-start file and system checks ... done." << endl << endl;

    // Check if output directories exist or create them
    logfile << gt() << "* Checking output directories..." << endl;
    if (!directoryExists("output")) {
        logfile << gt() << "  -> Main output directory does not exist, creating..." << endl;
        if (_mkdir("output") != 0) {
            logfile << gt() << "  -> ERROR: Failed to create output directory!" << endl;
        } else {
            logfile << gt() << "  -> Successfully created output directory." << endl;
        }
    } else {
        logfile << gt() << "  -> Main output directory exists." << endl;
    }

    logfile << "//---------------------------------------------------------------------//" << endl;
    logfile << "// Main simulation settings:" << endl;
    int rev = gapinput.data.lubrication_module.n_lubrication_revolutions;
    
    // Check if using temperature-based convergence
	bool temperature_based_convergence = gapinput.data.lubrication_module.temperature_convergence;
    if (temperature_based_convergence) {
        logfile << "//   Temperature convergence enabled (difference threshold: " << TEMPERATURE_CONVERGENCE_THRESHOLD << " C)" << endl;
        logfile << "//   Maximum number of lubrication module revolutions = " << abs(rev) << endl;
        // Convert to positive for internal use
        rev = abs(rev);
    } else {
        logfile << "//   Number of lubrication module revolutions = " << rev << endl;
    }

    vector<bool> interface_enabled(3, false);
    vector<bool> interface_disabled(3, false);
    int active_interfaces = 0;
    bool any_interfaces_active = false;
    
    if (gapinput.data.lubrication_module.solve_piston == 1) { 
        interface_enabled[0] = true; 
        active_interfaces++;
        any_interfaces_active = true;
        logfile << "//   Solve piston interface = true" << endl; 
        
        // Check for piston directory
        if (!directoryExists("output/piston")) {
            logfile << gt() << "  -> Piston output directory does not exist, creating..." << endl;
            if (_mkdir("output/piston") != 0) {
                logfile << gt() << "  -> ERROR: Failed to create piston output directory!" << endl;
            } else {
                if (_mkdir("output/piston/vtk") != 0) {
                    logfile << gt() << "  -> ERROR: Failed to create piston/vtk output directory!" << endl;
                } else {
                    logfile << gt() << "  -> Successfully created piston output directories." << endl;
                }
            }
        } else {
            logfile << gt() << "  -> Piston output directory exists." << endl;
            if (!directoryExists("output/piston/vtk")) {
                if (_mkdir("output/piston/vtk") != 0) {
                    logfile << gt() << "  -> ERROR: Failed to create piston/vtk output directory!" << endl;
                }
            }
        }
    }
    else { logfile << "//   Solve piston interface = false" << endl; }
    
    if (gapinput.data.lubrication_module.solve_block == 1) { 
        interface_enabled[1] = true; 
        active_interfaces++;
        any_interfaces_active = true;
        logfile << "//   Solve block interface = true" << endl; 
        
        // Check for block directory
        if (!directoryExists("output/block")) {
            logfile << gt() << "  -> Block output directory does not exist, creating..." << endl;
            if (_mkdir("output/block") != 0) {
                logfile << gt() << "  -> ERROR: Failed to create block output directory!" << endl;
            } else {
                logfile << gt() << "  -> Successfully created block output directory." << endl;
            }
        } else {
            logfile << gt() << "  -> Block output directory exists." << endl;
        }
    }
    else { logfile << "//   Solve block interface = false" << endl; }
    
    if (gapinput.data.lubrication_module.solve_slipper == 1) { 
        interface_enabled[2] = true; 
        active_interfaces++;
        any_interfaces_active = true;
        logfile << "//   Solve slipper interface = true" << endl; 
        
        // Check for slipper directory
        if (!directoryExists("output/slipper")) {
            logfile << gt() << "  -> Slipper output directory does not exist, creating..." << endl;
            if (_mkdir("output/slipper") != 0) {
                logfile << gt() << "  -> ERROR: Failed to create slipper output directory!" << endl;
            } else {
                if (_mkdir("output/slipper/vtk") != 0) {
                    logfile << gt() << "  -> ERROR: Failed to create slipper/vtk output directory!" << endl;
                } else {
                    logfile << gt() << "  -> Successfully created slipper output directories." << endl;
                }
            }
        } else {
            logfile << gt() << "  -> Slipper output directory exists." << endl;
            if (!directoryExists("output/slipper/vtk")) {
                if (_mkdir("output/slipper/vtk") != 0) {
                    logfile << gt() << "  -> ERROR: Failed to create slipper/vtk output directory!" << endl;
                }
            }
        }
    }
    else { logfile << "//   Solve slipper interface = false" << endl; }
    
    // Initialize temperature files after directory creation
    initializeTemperatureFiles(interface_enabled[0], interface_enabled[1], interface_enabled[2]);
    
    logfile << "//---------------------------------------------------------------------//" << endl << endl;

    // Check if no interfaces are active
    if (!any_interfaces_active) {
        logfile << gt() << "ERROR: No interfaces enabled for simulation. At least one interface must be enabled." << endl;
        logfile << gt() << "Exiting simulation." << endl;
        return -1;
    }

    task_group tasks;
    event* local_rev = new event[3];
    event* global_rev = new event[3];

    logfile << gt() << "* Creating worker threads ... done." << endl;
    logfile << gt() << "* Setting thread work pool scheduler policy ... done." << endl;
    CurrentScheduler::Create(SchedulerPolicy(1, MaxConcurrency, 3));

    logfile << gt() << "* Standby for worker thread startup ... " << endl;
    if (interface_enabled[0]) { tasks.run([&]() { piston_GUI(local_rev[0], global_rev[0]); }); logfile << gt() << "   -> Launching piston lubrication worker ... done." << endl; }
    if (interface_enabled[1]) { tasks.run([&]() { run_block_gap(local_rev[1], global_rev[1], gapinput); }); logfile << gt() << "   -> Launching block lubrication worker ... done." << endl; }
    if (interface_enabled[2]) { tasks.run([&]() { run_slipper_gap(local_rev[2], global_rev[2], gapinput); }); logfile << gt() << "   -> Launching slipper lubrication worker ... done." << endl; }

    logfile << gt() << "* Worker thread startup complete." << endl;
    logfile << endl << gt() << "* Starting main program loop ..." << endl << endl;

    // Initialize calculation status for all enabled interfaces
    for (int i = 0; i < 3; i++) {
        if (interface_enabled[i]) {
            startInterfaceCalculation(i);
        }
    }
	
    // Variables for timing the console updates
    static clock_t last_title_update = 0;
    static clock_t last_progress_poll = 0;
    static clock_t last_progress_display = clock(); // For periodic 10-second updates

    // Run the main simulation loop
    for(int i = 0; true; i++)
    {
        vector<bool> should_wait(3, true);

        // Wait for the original mechanism to do its work first before trying to read files
        for(int j = 0; true; j++)
        {
            if(j >= 3) j = 0;
            if(!interface_enabled[j] || interface_disabled[j]) continue;
        
            if(should_wait[j])
            {
                // Check if it's time to update the progress display while waiting
                clock_t current_time = clock();
                
                // Poll for file updates every 100ms
                if (current_time - last_progress_poll > CLOCKS_PER_SEC / 10) {
                    last_progress_poll = current_time;
                    // Update progress from files
                    updateProgressFromFiles(interface_enabled[0], interface_enabled[1], interface_enabled[2],
                                          interface_disabled[0], interface_disabled[1], interface_disabled[2],
                                          rev, (current_time - last_title_update > CLOCKS_PER_SEC / 2));
                    
                    // If we updated the console title, update the timestamp
                    if (current_time - last_title_update > CLOCKS_PER_SEC / 2) {
                        last_title_update = current_time;
                    }
                }
                
                // Check if we should update estimated progress (every 10 seconds)
                if (current_time - last_progress_display > PROGRESS_UPDATE_INTERVAL * CLOCKS_PER_SEC) {
                    last_progress_display = current_time;
                    
                    // Update estimated progress for active interfaces
                    for (int k = 0; k < 3; k++) {
                        if (interface_enabled[k] && !interface_disabled[k]) {
                            updateEstimatedProgress(k, rev);
                        }
                    }
                    
					
                    // Display the estimated progress
                    displayEstimatedProgress(rev, logfile, interface_enabled, interface_disabled);
                    
                    // Update console title with estimated progress
                    updateConsoleTitle(interface_enabled[0], interface_enabled[1], interface_enabled[2],
                                     interface_disabled[0], interface_disabled[1], interface_disabled[2],
                                     true); // Use estimated progress
                }
            
                // Wait with timeout for the event
                if(local_rev[j].wait(100) != COOPERATIVE_WAIT_TIMEOUT) // Shorter timeout to allow more frequent progress updates
                {
                    // Interface completed a step
                    should_wait[j] = false;
                    endInterfaceCalculation(j); // Mark calculation as complete
                    logfile.flush();
                    local_rev[j].reset();
                    
                    // Start timing for the next calculation
                    startInterfaceCalculation(j);
                }
            }
        
            bool can_break = true;
            for(int k=0; k<3; k++)
            {
                if(interface_enabled[k] && !interface_disabled[k] && should_wait[k]) can_break = false;
            }
        
            if(can_break) break;
			
        }
		std::cout << std::endl;
        // ---------- Process output files and create summary ----------
        ostringstream summary;
        string timeStr = gt();

        //--- Process piston.txt with progress bar ---
		if (interface_enabled[0] && !interface_disabled[0]) {
			// Ensure directory exists before trying to read the file
			if (directoryExists("output/piston")) {
				string lastLine;
				if (readLastLineFromFile("output/piston/piston.txt", lastLine)) {
					istringstream iss(lastLine);
					vector<double> vals;
					copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(vals));

					if (vals.size() >= 62) {
						// Get the revolution and phi from the data
						int rev_value = static_cast<int>(vals[1]);
						double phi_value = vals[2];
                
						// Update progress information
						updateInterfaceProgress(0, rev_value, phi_value, rev);
                
						// Calculate clearance percentages
						double e1 = vals[3];   // e_1
						double e2 = vals[4];   // e_2  
						double e3 = vals[5];   // e_3
						double e4 = vals[6];   // e_4
                
						// Get geometry values from input
						double dBushing = gapinput.data.geometry.dZ;  // Bushing diameter
						double dPiston = gapinput.data.geometry.dK;   // Piston diameter
						double clearance_total = (dBushing - dPiston) / 2.0;
                
						double C1 = (sqrt(e1*e1 + e2*e2)) / clearance_total * 100.0;
						double C2 = (sqrt(e3*e3 + e4*e4)) / clearance_total * 100.0;
                
						double QL = vals[60];
						double ML = vals[61];
                
						// Get progress bar - use actual progress, not estimated
						string progressBar = getProgressBar(interface_progress[0].progress_percentage);
                
						// Updated output with clearance percentages
						summary << timeStr << "Piston " << progressBar << " Rev: " << rev_value 
								<< "+" << fixed << setprecision(0) << (phi_value / 3.6) << "% Rev " 
								<< fixed << setprecision(2) << (rev_value + phi_value/360.0) << " of " << rev 
								<< ", C1: " << fixed << setprecision(0) << C1 << "% C2: " 
								<< C2 << "%, QL: " << QL << " l/min, ML: " << ML << " Nm" << endl;
											}
				}
			}
		}

        //--- Process block.txt with progress bar ---
        if (interface_enabled[1] && !interface_disabled[1]) {
            // Check if directory exists
            if (directoryExists("output/block")) {
                string lastLine;
                if (readLastLineFromFile("output/block/block.txt", lastLine)) {
                    istringstream iss(lastLine);
                    vector<double> vals; 
                    copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(vals));

                    if (vals.size() >= 23) {
                        // Get the revolution and phi values
                        int rev_value = static_cast<int>(vals[1]);
                        double phi_value = vals[2];
                        
                        // Update progress information
                        updateInterfaceProgress(1, rev_value, phi_value, rev);
                        
                        double hmin = vals[12];
                        double hmax = vals[13];
                        double hmean = vals[14];
                        double QL = vals[21];
                        double ML = vals[22];
                        
                        // Get progress bar - use actual progress, not estimated
                        string progressBar = getProgressBar(interface_progress[1].progress_percentage);

                        summary << timeStr << "Block  " << progressBar << " " << fixed << setprecision(0) 
								<< (phi_value / 3.6) << "% Rev " << fixed << setprecision(2) 
								<< (rev_value + phi_value/360.0) << " of " << rev << ", h (" << hmin 
								<< " / " << hmean << " / " << hmax << " mu), QL: " << QL 
								<< " l/min, ML: " << ML << " Nm" << endl;
                    }
                }
            }
        }

        //--- Process slipper.txt with progress bar ---
        if (interface_enabled[2] && !interface_disabled[2]) {
            // Check if directory exists
            if (directoryExists("output/slipper")) {
                string lastLine;
                if (readLastLineFromFile("output/slipper/slipper.txt", lastLine)) {
                    istringstream iss(lastLine);
                    vector<double> vals; 
                    copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(vals));

                    if (vals.size() >= 14) {
                        // Get the revolution and phi values
                        int rev_value = static_cast<int>(vals[1]);
                        double phi_value = vals[2];
                        
                        // Update progress information
                        updateInterfaceProgress(2, rev_value, phi_value, rev);
                        
                        double hmin = vals[9];
                        double hmean = vals[10];
                        double hmax = vals[11];
                        double QL = vals[12];
                        double ML = vals[13];
                        
                        // Get progress bar - use actual progress, not estimated
                        string progressBar = getProgressBar(interface_progress[2].progress_percentage);

                        summary << timeStr << "Slipper " << progressBar << " " << fixed << setprecision(0) 
								<< (phi_value / 3.6) << "% Rev " << fixed << setprecision(2) 
								<< (rev_value + phi_value/360.0) << " of " << rev << ", h (" << hmin 
								<< " / " << hmean << " / " << hmax << " mu), QL: " << QL 
								<< " l/min, ML: " << ML << " Nm" << endl;
                    }
                }
            }
        }
        
        // Update console title with latest progress
        updateConsoleTitle(interface_enabled[0], interface_enabled[1], interface_enabled[2],
                         interface_disabled[0], interface_disabled[1], interface_disabled[2]);
        
        // Output the summary to log file
        logfile << summary.str();
        logfile.flush();
        
        // STEP 1: Log thermal data for active interfaces only (with temperature file logging)
        thermal_log_per_revolution_with_files(i, logfile, 
                                             interface_enabled[0] && !interface_disabled[0], 
                                             interface_enabled[1] && !interface_disabled[1], 
                                             interface_enabled[2] && !interface_disabled[2]);

        // STEP 2: Check temperature convergence between consecutive revolutions if enabled
        if (temperature_based_convergence && i > 0) {
            check_temperature_convergence(i, interface_enabled, interface_disabled, active_interfaces, logfile);
        }

        // STEP 3: Modified convergence check
        bool should_exit_loop = false;
        if (temperature_based_convergence) {
            if (active_interfaces == 0) {
                logfile << gt() << "* Temperature convergence achieved - all interfaces have converged." << endl;
                logfile << gt() << "* Simulation completed after " << (i+1) << " revolutions." << endl;
                logfile.flush();
                should_exit_loop = true;
            } 
            else if (i >= rev-1) {
                // Reached max revolutions but not all interfaces converged
                logfile << gt() << "* WARNING: Reached maximum revolutions (" << rev << ") without complete temperature convergence." << endl;
                for (int j = 0; j < 3; j++) {
                    if (interface_enabled[j] && !interface_disabled[j]) {
                        logfile << gt() << "* Interface " << i2word(j) << " did not achieve temperature convergence." << endl;
                    }
                }
                logfile << gt() << "* Simulation completed after " << (i+1) << " revolutions." << endl;
                logfile.flush();
                should_exit_loop = true;
            }
        }
        else if (!temperature_based_convergence && i >= rev-1) {
            // Non-temperature based convergence, just check revolutions
            logfile << gt() << "* Reached target revolutions (" << rev << "), exiting loop" << endl;
            logfile.flush();
            should_exit_loop = true;
        }

        // Release the threads for the next revolution, but only for active interfaces
        for(int j = 0; j < 3; j++)
        {
            if(!interface_enabled[j] || interface_disabled[j]) continue;
        
            logfile.flush();
            global_rev[j].set();
        }

        // Log a newline for formatting
        logfile << endl;
        
        // Exit the loop after all processing is complete
        if (should_exit_loop) {
			cout << endl;
            break;
        }
		
    }
    
    // Run final efficiency analysis
    if (rev > 0) {
        logfile << gt() << "Final efficiency analysis after last revolution..." << std::endl;
        analyzeEfficiencyPerRevolution(rev - 1, logfile, 
                                      interface_enabled[0] && !interface_disabled[0],
                                      interface_enabled[1] && !interface_disabled[1],
                                      interface_enabled[2] && !interface_disabled[2]);
    }

	// Generate comprehensive pump summary report
    generatePumpSummaryReport(rev - 1, logfile, 
                             interface_enabled[0], interface_enabled[1], interface_enabled[2],
                             interface_disabled[0], interface_disabled[1], interface_disabled[2]);

    logfile << gt() << "* Worker threads have completed all revolutions." << endl;
    logfile << gt() << "* Standby for worker thread shutdown ..."; logfile.flush();
    tasks.wait(); logfile << " done." << endl;
    logfile << gt() << "* Removing temporary files ..."; system("rd /S /Q .\\temp"); logfile << " done." << endl;
    logfile << gt() << "* Preforming post-shutdown file and system checks ..."; delete[] local_rev; delete[] global_rev; logfile << " done." << endl;
    logfile << endl << gt() << "* Standby for main simulation thread termination ... " << endl;
    shouldexit = true;

    return 0;
}


